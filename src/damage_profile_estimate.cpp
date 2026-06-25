// Core damage estimation: SamplePassState, run_sample_pass, estimate_damage.
// damage_profile.cpp
// estimate_damage() — deamination + QC full pipeline damage estimation.
// Called by both 'derep' (Pass 0) and 'extend' (Pass 0).

#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"

#include "taph/frame_selector_decl.hpp"
#include "taph/length_gc_joint_mixture.hpp"
#include "taph/length_stratified_profile.hpp"
#include "taph/library_interpretation.hpp"
#include "taph/log_length_gmm.hpp"
// taph/sample_damage_profile.hpp is included transitively via damage_profile.hpp

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

static constexpr double MIN_COV_DP = 100.0;

// Pass-state captured by the unified sample/profile pass. The QC path
// also reads hex3_terminal, hex3_count, and short_read_frac from here.
struct SamplePassState {
    taph::SampleDamageProfile deam_profile;
    int64_t reads_scanned = 0;
    int64_t record_pos    = 0;
    int64_t len_sum       = 0;
    int64_t short_reads   = 0;   // L<30 (not contributing to deam_profile but seen)
    std::array<uint32_t, 4096> hex3_terminal{};
    uint64_t                   hex3_count = 0;
    int64_t                    adapter_reads_scanned = 0;
    std::array<uint64_t, LSD_HIST_BINS> lsd_hist{};
};

// Per-record accumulation, shared by the serial and parallel paths. Touches only the count state
// (commutative), so summing thread-local states via merge_sample_profiles reproduces the serial
// result bit-for-bit. record_pos / max_reads gating is the caller's job. `buf` is a scratch reused
// across calls to avoid per-record allocation under clipping.
static inline void accumulate_record(SamplePassState& s, const std::string& seq,
                                     const std::vector<std::string>& clip_stubs5,
                                     const std::vector<std::string>& clip_stubs3,
                                     int64_t adapter_scan_reads, std::string& buf)
{
    const std::string* seq_p = &seq;
    if (!clip_stubs5.empty() || !clip_stubs3.empty()) {
        buf = seq;
        if (!clip_stubs5.empty() && (int)buf.size() >= 6) {
            for (const auto& stub : clip_stubs5)
                if (buf.compare(0, 6, stub) == 0) { buf.erase(0, 6); break; }
        }
        if (!clip_stubs3.empty()) {
            bool trimmed;
            do {
                trimmed = false;
                int L = static_cast<int>(buf.size());
                if (L < 12) break;
                for (const auto& stub : clip_stubs3) {
                    if (buf.compare(L - 6, 6, stub) == 0) {
                        buf.erase(L - 6, 6);
                        trimmed = true;
                        break;
                    }
                }
            } while (trimmed);
        }
        seq_p = &buf;
    }

    int L = static_cast<int>(seq_p->size());
    if (L < 30) {
        if (L > 0) s.short_reads++;
        return;
    }
    s.len_sum += L;
    ++s.lsd_hist[lsd_hist_bin(L)];
    taph::FrameSelector::update_sample_profile(s.deam_profile, *seq_p);
    s.reads_scanned++;

    if (adapter_scan_reads == 0 || s.adapter_reads_scanned < adapter_scan_reads) {
        if (L >= 12) {
            int code = taph::encode_hex_at(*seq_p, L - 6);
            if (code >= 0) { ++s.hex3_terminal[code]; ++s.hex3_count; }
        }
        ++s.adapter_reads_scanned;
    }
}

static void merge_pass_state(SamplePassState& dst, SamplePassState& src)
{
    taph::FrameSelector::merge_sample_profiles(dst.deam_profile, src.deam_profile);
    dst.reads_scanned         += src.reads_scanned;
    dst.short_reads           += src.short_reads;
    dst.len_sum               += src.len_sum;
    dst.adapter_reads_scanned += src.adapter_reads_scanned;
    dst.hex3_count            += src.hex3_count;
    for (size_t i = 0; i < dst.hex3_terminal.size(); ++i) dst.hex3_terminal[i] += src.hex3_terminal[i];
    for (size_t i = 0; i < dst.lsd_hist.size();      ++i) dst.lsd_hist[i]      += src.lsd_hist[i];
}

static void run_sample_pass(FastqReaderBase& reader,
                            taph::SampleDamageProfile::LibraryType forced_lib,
                            int64_t max_reads,
                            int64_t adapter_scan_reads,
                            const std::vector<std::string>& clip_stubs5,
                            const std::vector<std::string>& clip_stubs3,
                            SamplePassState& s,
                            size_t threads = 1)
{
    s.deam_profile.forced_library_type = forced_lib;

    if (threads <= 1) {
        FastqRecord rec;
        std::string buf;
        while (reader.read(rec)) {
            s.record_pos++;
            if (max_reads > 0 && s.record_pos > max_reads) break;
            accumulate_record(s, rec.seq, clip_stubs5, clip_stubs3, adapter_scan_reads, buf);
        }
        return;
    }

    // Parallel: one reader thread (serial gzip decode) feeds batches to N workers, each accumulating
    // into a thread-local state; merged afterwards. update_sample_profile sums commutative counts, so
    // the merged estimate is identical to the serial one. Mirrors the LSD-scan producer/consumer.
    const int n = static_cast<int>(threads);
    constexpr int BATCH = 4096;
    struct Queue {
        std::mutex mtx;
        std::condition_variable cv_fill, cv_drain;
        std::vector<std::vector<std::string>> q;
        bool done = false;
        int cap;
        explicit Queue(int c) : cap(c) {}
    } queue(n * 2);

    std::vector<SamplePassState> wstate(n);
    for (auto& w : wstate) w.deam_profile.forced_library_type = forced_lib;
    // Split the adapter-scan budget across workers so the aggregate stays ~adapter_scan_reads
    // (0 = scan every read, as in the serial path).
    const int64_t per_worker_adapter =
        adapter_scan_reads > 0 ? (adapter_scan_reads + n - 1) / n : 0;

    std::vector<std::thread> workers;
    workers.reserve(n);
    for (int t = 0; t < n; ++t) {
        workers.emplace_back([&, t] {
            std::vector<std::string> batch;
            std::string buf;
            for (;;) {
                {
                    std::unique_lock<std::mutex> lk(queue.mtx);
                    queue.cv_fill.wait(lk, [&] { return !queue.q.empty() || queue.done; });
                    if (queue.q.empty()) break;
                    batch = std::move(queue.q.back());
                    queue.q.pop_back();
                    queue.cv_drain.notify_one();
                }
                for (auto& seq : batch)
                    accumulate_record(wstate[t], seq, clip_stubs5, clip_stubs3,
                                      per_worker_adapter, buf);
            }
        });
    }

    int64_t total = 0;
    FastqRecord rec;
    std::vector<std::string> batch;
    batch.reserve(BATCH);
    while (reader.read(rec)) {
        total++;
        if (max_reads > 0 && total > max_reads) break;  // matches serial: cap-th+1 record dropped
        batch.push_back(std::move(rec.seq));
        if (static_cast<int>(batch.size()) == BATCH) {
            std::unique_lock<std::mutex> lk(queue.mtx);
            queue.cv_drain.wait(lk, [&] { return static_cast<int>(queue.q.size()) < queue.cap; });
            queue.q.push_back(std::move(batch));
            queue.cv_fill.notify_one();
            batch.clear();
            batch.reserve(BATCH);
        }
    }
    if (!batch.empty()) {
        std::unique_lock<std::mutex> lk(queue.mtx);
        queue.cv_drain.wait(lk, [&] { return static_cast<int>(queue.q.size()) < queue.cap; });
        queue.q.push_back(std::move(batch));
        queue.cv_fill.notify_one();
    }
    {
        std::unique_lock<std::mutex> lk(queue.mtx);
        queue.done = true;
        queue.cv_fill.notify_all();
    }
    for (auto& w : workers) w.join();

    s.record_pos = total;
    for (auto& w : wstate) merge_pass_state(s, w);
}

static DamageProfile estimate_damage_impl(
    FastqReaderBase& reader,
    const std::string& path,
    double mask_threshold,
    taph::SampleDamageProfile::LibraryType forced_lib,
    int64_t max_reads,
    SamplePassState* out_state = nullptr,
    int64_t adapter_scan_reads = 0,
    const std::vector<std::string>* clip_stubs5 = nullptr,
    const std::vector<std::string>* clip_stubs3 = nullptr,
    size_t threads = 1)
{
    std::unique_ptr<SamplePassState> local_owner;
    if (!out_state) local_owner = std::make_unique<SamplePassState>();
    SamplePassState& s = out_state ? *out_state : *local_owner;
    static const std::vector<std::string> kEmptyStubs;
    run_sample_pass(reader, forced_lib, max_reads, adapter_scan_reads,
                    clip_stubs5 ? *clip_stubs5 : kEmptyStubs,
                    clip_stubs3 ? *clip_stubs3 : kEmptyStubs,
                    s, threads);

    int      typical_len   = 0;
    int64_t& reads_scanned = s.reads_scanned;
    int64_t& record_pos    = s.record_pos;
    int64_t& len_sum       = s.len_sum;
    taph::SampleDamageProfile& deam_profile = s.deam_profile;

    if (reads_scanned > 0)
        typical_len = static_cast<int>(len_sum / reads_scanned);

    taph::FrameSelector::finalize_sample_profile(deam_profile);

    double d_max_5   = deam_profile.d_max_5prime;
    double d_max_3   = deam_profile.d_max_3prime;
    double lambda_5  = deam_profile.lambda_5prime;
    double lambda_3  = deam_profile.lambda_3prime;
    double bg_5      = deam_profile.fit_baseline_5prime;
    double bg_3      = deam_profile.fit_baseline_3prime;
    double background = (bg_5 + bg_3) / 2.0;
    double d_max_combined = deam_profile.d_max_combined;

    // P0-2 fix: UNKNOWN must NOT silently fall back to DS or to a posterior
    // pick — both produced biased masks. Map library_type → Damage3PrimeMode
    // unconditionally; UNKNOWN → neutral (3' damage prior dropped entirely).
    Damage3PrimeMode mode_3p;
    switch (deam_profile.library_type) {
        case taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED:
            mode_3p = Damage3PrimeMode::ds_ga; break;
        case taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED:
            mode_3p = Damage3PrimeMode::ss_ct; break;
        case taph::SampleDamageProfile::LibraryType::UNKNOWN:
        default:
            mode_3p = Damage3PrimeMode::neutral; break;
    }
    bool is_ss = (mode_3p == Damage3PrimeMode::ss_ct);

    DamageProfile profile;
    profile.d_max_5prime   = d_max_5;
    profile.d_max_3prime   = d_max_3;
    profile.lambda_5prime  = lambda_5;
    profile.lambda_3prime  = lambda_3;
    profile.background     = background;
    profile.mask_threshold = mask_threshold;
    profile.bg_5_tc        = deam_profile.fit_baseline_5prime;
    profile.bg_3_channel   = deam_profile.fit_baseline_3prime;
    profile.mixture_d_damaged = deam_profile.mixture_d_damaged;
    profile.mixture_pi_damaged = deam_profile.mixture_pi_damaged;
    profile.mixture_d_population_highgc = deam_profile.mixture_d_population_highgc;
    profile.mixture_n_components = deam_profile.mixture_n_components;
    profile.mixture_converged = deam_profile.mixture_converged;
    profile.mixture_identifiable = deam_profile.mixture_identifiable;
    profile.pi_point          = deam_profile.pi.point;  // SOLUTION §6.3 validated pi (shadow step 2)
    profile.pi_lo             = deam_profile.pi.lo;
    profile.pi_hi             = deam_profile.pi.hi;
    profile.pi_state          = deam_profile.pi.state;
    profile.pi_detected       = (deam_profile.pi.state == taph::DamageConfidence::DETECTED);
    profile.mode_3prime    = mode_3p;
    profile.ss_mode        = (mode_3p == Damage3PrimeMode::ss_ct);  // back-compat
    profile.enabled        = (d_max_combined > 0.02 || d_max_5 > 0.02 || d_max_3 > 0.02
                              || deam_profile.damage_validated);

    // Empirical per-position mask from the deamination model'''s normalized frequencies.
    // DS: max(5' C→T excess, 3' G→A excess) > threshold.
    // SS: max(5' C→T excess, 3' T/(T+C) excess) > threshold — both ends have C→T.
    double bg_tc = deam_profile.baseline_t_freq /
                   (deam_profile.baseline_t_freq + deam_profile.baseline_c_freq + 1e-9);
    for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (deam_profile.tc_total_5prime[p] >= MIN_COV_DP)
            excess_5 = deam_profile.t_freq_5prime[p] - bg_5;
        // P0-2: neutral mode skips the 3' channel entirely; mask comes from 5' only.
        if (mode_3p == Damage3PrimeMode::ss_ct) {
            if (deam_profile.tc_total_3prime[p] >= MIN_COV_DP) {
                double tc3 = deam_profile.tc_total_3prime[p];
                double ct3_freq = (tc3 > 0)
                    ? deam_profile.t_freq_3prime[p] / tc3
                    : bg_tc;
                excess_3 = ct3_freq - bg_tc;
            }
        } else if (mode_3p == Damage3PrimeMode::ds_ga) {
            if (deam_profile.ag_total_3prime[p] >= MIN_COV_DP)
                excess_3 = deam_profile.a_freq_3prime[p] - bg_3;
        }
        profile.mask_pos[p] = (excess_5 > mask_threshold) || (excess_3 > mask_threshold);
    }

    // Fill short gaps in the empirical damage mask (≤1 position) to handle
    // sampling noise creating non-contiguous masks.
    {
        constexpr int MAX_GAP = 1;
        int last_p = -1;
        for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p) {
            if (profile.mask_pos[p]) {
                if (last_p >= 0 && (p - last_p - 1) <= MAX_GAP) {
                    for (int g = last_p + 1; g < p; ++g)
                        profile.mask_pos[g] = true;
                }
                last_p = p;
            }
        }
    }

    std::string sample_note = (max_reads > 0 && record_pos > max_reads)
        ? " [sampled " + std::to_string(max_reads) + " reads]" : "";
    log_info("Pass 0: deamination damage estimation — " +
             std::to_string(reads_scanned) + " reads in " + path + sample_note);
    log_info("  Library type: " + std::string(deam_profile.library_type_str()) +
             (deam_profile.library_type_auto_detected ? " (auto-detected)" : " (forced)"));
    log_info("  5'-end: d_max=" + std::to_string(d_max_5) +
             " lambda=" + std::to_string(lambda_5) +
             " bg=" + std::to_string(bg_5));
    log_info("  3'-end: d_max=" + std::to_string(d_max_3) +
             " lambda=" + std::to_string(lambda_3) +
             " bg=" + std::to_string(bg_3));
    log_info("  d_max_combined=" + std::to_string(d_max_combined) +
             " (source=" + deam_profile.d_max_source_str() + ")" +
             (deam_profile.mixture_converged ? " [mixture]" : "") +
             (deam_profile.damage_validated  ? " [validated]" : "") +
             (deam_profile.damage_artifact   ? " [ARTIFACT]" : ""));
    if (profile.enabled && typical_len > 0) {
        profile.print_info(typical_len);
    } else {
        log_info("  Damage below threshold — standard exact hashing will be used");
    }
    return profile;
}

DamageProfile estimate_damage(
    const std::string& path,
    double mask_threshold,
    taph::SampleDamageProfile::LibraryType forced_lib,
    int64_t max_reads)
{
    DamageEstimateOptions o;
    o.mask_threshold = mask_threshold;
    o.forced_lib     = forced_lib;
    o.max_reads      = max_reads;
    o.qc_enabled     = false;
    o.clip_policy    = DamageClipPolicy::Off;
    return estimate_damage_with_qc(path, o).profile;
}

DamageEstimate estimate_damage_with_qc(const std::string& path,
                                       const DamageEstimateOptions& opts)
{
    DamageEstimate out;

    SamplePassState state;
    {
        auto reader = make_fastq_reader(path);
        out.profile = estimate_damage_impl(*reader, path,
                                           opts.mask_threshold,
                                           opts.forced_lib,
                                           opts.max_reads,
                                           &state,
                                           opts.qc_enabled ? opts.adapter_scan_reads : 0,
                                           nullptr, nullptr, opts.threads);
    }
    out.lsd_hist.assign(state.lsd_hist.begin(), state.lsd_hist.end());

    // Bug 3: small-input warning. deamination fit needs ~100k reads to converge to
    // d_max>0; below that the fit silently returns 0 and downstream code
    // sees "no damage". Warn when the user actually asked for the fit.
    if (!opts.skip_deam_fit) {
        const int64_t n_seen = state.reads_scanned + state.short_reads;
        if (n_seen > 0 && n_seen < 100000) {
            log_warn("WARNING: damage profile fit needs >=100k reads to converge; got " +
                     std::to_string(n_seen) + " — d_max may be 0");
        }
        const double dmax = std::max({
            static_cast<double>(out.profile.d_max_5prime),
            static_cast<double>(out.profile.d_max_3prime),
            static_cast<double>(state.deam_profile.d_max_combined)});
        if (n_seen >= 100000 && !(dmax > 0.0)) {
            log_warn("deamination fit returned d_max=0; QC adapter/hexamer scan still ran");
        }
    }

    // Bug 2: QC-only mode — caller wants adapter/hex stats + QC flags but
    // does NOT want deamination damage masking applied to dedup hashing. Replace
    // the fitted profile with a stub (enabled=false) and continue into the
    // QC blocks below.
    if (opts.skip_deam_fit) {
        DamageProfile stub;
        stub.mask_threshold = opts.mask_threshold;
        out.profile = stub;
    }

    if (!opts.qc_enabled) {
        out.qc.enabled = false;
        return out;
    }

    out.qc.enabled               = true;
    out.qc.adapter_reads_scanned = state.adapter_reads_scanned;

    // Adapter-stub detection from finalized dp + per-pass hex3 terminal table.
    out.qc.adapter = taph::detect_adapter_stubs(state.deam_profile,
                                                state.hex3_terminal.data(),
                                                state.hex3_count);

    // Optional refit pass: re-profile with detected stubs stripped.
    if (opts.clip_policy == DamageClipPolicy::Refit &&
        (out.qc.adapter.adapter_clipped || out.qc.adapter.adapter3_clipped))
    {
        SamplePassState clipped;
        auto reader2 = make_fastq_reader(path);
        out.profile = estimate_damage_impl(*reader2, path,
                                           opts.mask_threshold,
                                           opts.forced_lib,
                                           opts.max_reads,
                                           &clipped,
                                           opts.adapter_scan_reads,
                                           &out.qc.adapter.stubs5,
                                           &out.qc.adapter.stubs3,
                                           opts.threads);
        out.qc.profile_clipped = true;
        // Recompute top enriched hexamers from the post-clip dp.
        out.qc.adapter.top_enriched =
            taph::compute_hex_enriched_5prime(clipped.deam_profile);
        out.qc.adapter.flag_hex_artifact = !out.qc.adapter.top_enriched.empty()
            && out.qc.adapter.top_enriched[0].log2fc > 1.5
            && !out.qc.adapter.top_enriched[0].damage_consistent;
        // Hex stats / preservation use the post-clip dp.
        out.qc.hex_stats = taph::compute_hex_stats(clipped.deam_profile);

        const bool is_ss = (clipped.deam_profile.library_type ==
            taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
        const int64_t total = clipped.reads_scanned + clipped.short_reads;
        out.qc.short_read_frac = total > 0
            ? static_cast<double>(clipped.short_reads) / static_cast<double>(total)
            : 0.0;

        auto cpg = taph::compute_cpg_score(clipped.deam_profile);
        out.qc.flags = taph::compute_library_qc_flags(
            clipped.deam_profile, is_ss, out.qc.adapter.flag_hex_artifact,
            out.qc.hex_stats.jsd, out.qc.hex_stats.entropy_terminal,
            out.qc.short_read_frac);
        out.qc.preservation = taph::compute_preservation_summary(
            clipped.deam_profile, is_ss,
            out.qc.adapter.adapter_clipped, out.qc.adapter.flag_hex_artifact,
            cpg.z, /*oxog_score_z=*/0.0,
            std::numeric_limits<double>::quiet_NaN(),
            out.qc.hex_stats.shift_p);
    } else {
        // Report-only: stats from the original sample pass.
        out.qc.hex_stats = taph::compute_hex_stats(state.deam_profile);

        const bool is_ss = (state.deam_profile.library_type ==
            taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
        const int64_t total = state.reads_scanned + state.short_reads;
        out.qc.short_read_frac = total > 0
            ? static_cast<double>(state.short_reads) / static_cast<double>(total)
            : 0.0;

        auto cpg = taph::compute_cpg_score(state.deam_profile);
        out.qc.flags = taph::compute_library_qc_flags(
            state.deam_profile, is_ss, out.qc.adapter.flag_hex_artifact,
            out.qc.hex_stats.jsd, out.qc.hex_stats.entropy_terminal,
            out.qc.short_read_frac);
        out.qc.preservation = taph::compute_preservation_summary(
            state.deam_profile, is_ss,
            out.qc.adapter.adapter_clipped, out.qc.adapter.flag_hex_artifact,
            cpg.z, /*oxog_score_z=*/0.0,
            std::numeric_limits<double>::quiet_NaN(),
            out.qc.hex_stats.shift_p);
    }

    // String-ify enabled QC flags for downstream JSON consumers.
    auto& f = out.qc.flags;
    auto& fn = out.qc.flag_names;
    if (f.adapter_remnant_5prime)   fn.emplace_back("adapter_remnant_5prime");
    if (f.adapter_remnant_3prime)   fn.emplace_back("adapter_remnant_3prime");
    if (f.hexamer_composition_bias) fn.emplace_back("hexamer_composition_bias");
    if (f.hexamer_terminal_shift)   fn.emplace_back("hexamer_terminal_shift");
    if (f.short_read_spike)         fn.emplace_back("short_read_spike");
    if (f.depurination)             fn.emplace_back("depurination");
    if (f.ds_3prime_signal_absent)  fn.emplace_back("ds_3prime_signal_absent");
    if (f.ga3_inward_displaced)     fn.emplace_back("ga3_inward_displaced");
    if (f.hexamer_artifact_bias)    fn.emplace_back("hexamer_artifact_bias");

    return out;
}

