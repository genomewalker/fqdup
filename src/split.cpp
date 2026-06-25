// fqdup split — classify reads as damaged/undamaged without deduplication.
// Runs damage profile estimation (Pass 0), builds the LLR split model,
// then streams the input once routing each read to --out-damaged or --out-undamaged.
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <memory>
#include <string>
#include <thread>

#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"
#include "fqdup/logger.hpp"
#include "fqdup/version.hpp"

static void usage(const char* prog) {
    std::cerr <<
        "Usage: fqdup split -i INPUT [options]\n"
        "\nClassify reads as damaged or undamaged without deduplication.\n"
        "Estimates damage profile, builds per-read LLR split model, streams\n"
        "input once routing each read to --out-damaged / --out-undamaged.\n"
        "\nRequired:\n"
        "  -i FILE              Input FASTQ (raw or .gz); does NOT need to be sorted\n"
        "\nOutputs (at least one required):\n"
        "  --out-damaged  FILE  Write damaged reads here\n"
        "  --out-undamaged FILE Write undamaged reads here\n"
        "\nOptions:\n"
        "  --library-type STR   ds (default), ss, or auto (inferred from profile)\n"
        "  --split-model STR    auto (default), bulk, or empirical\n"
        "                         auto: empirical when d_max5>0.01, else bulk\n"
        "                         bulk: exponential model from Pass 0 estimate\n"
        "                         empirical: always run length-stratified LSD scan\n"
        "  --split-threshold F  LLR threshold for damaged call (default: 0.0)\n"
        "  --damage-deam-sample N  Max reads for Pass 0 damage scan (default: 5000000)\n"
        "  -t N                 Threads (default: all available, capped at 16 for I/O)\n"
        "  -h, --help           Show this help\n";
}

int split_main(int argc, char** argv) {
    std::string in_path, out_damaged_path, out_undamaged_path;
    float split_threshold = 0.0f;
    int64_t damage_deam_max_reads = 5'000'000;
    unsigned threads = std::max(1u, std::thread::hardware_concurrency());

    enum class SplitModelMode { Auto, Bulk, Empirical } split_model_mode = SplitModelMode::Auto;
    auto forced_library_type = taph::SampleDamageProfile::LibraryType::UNKNOWN;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if ((a == "-i" || a == "--input") && i+1 < argc)       { in_path = argv[++i]; }
        else if (a == "--out-damaged"  && i+1 < argc)           { out_damaged_path = argv[++i]; }
        else if (a == "--out-undamaged" && i+1 < argc)          { out_undamaged_path = argv[++i]; }
        else if (a == "--split-threshold" && i+1 < argc)        { split_threshold = std::stof(argv[++i]); }
        else if (a == "--damage-deam-sample" && i+1 < argc)     { damage_deam_max_reads = std::stoll(argv[++i]); }
        else if ((a == "-t" || a == "--threads") && i+1 < argc) { threads = static_cast<unsigned>(std::stoi(argv[++i])); }
        else if (a == "--library-type" && i+1 < argc) {
            std::string m = argv[++i];
            if      (m == "ds")   forced_library_type = taph::SampleDamageProfile::LibraryType::DOUBLE_STRANDED;
            else if (m == "ss")   forced_library_type = taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED;
            else if (m == "auto") forced_library_type = taph::SampleDamageProfile::LibraryType::UNKNOWN;
            else { std::cerr << "Error: unknown --library-type '" << m << "'\n"; return 1; }
        }
        else if (a == "--split-model" && i+1 < argc) {
            std::string m = argv[++i];
            if      (m == "auto")      split_model_mode = SplitModelMode::Auto;
            else if (m == "bulk")      split_model_mode = SplitModelMode::Bulk;
            else if (m == "empirical") split_model_mode = SplitModelMode::Empirical;
            else { std::cerr << "Error: unknown --split-model '" << m << "'\n"; return 1; }
        }
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else { std::cerr << "Error: unknown option '" << a << "'\n"; usage(argv[0]); return 1; }
    }

    if (in_path.empty()) {
        std::cerr << "Error: -i INPUT required\n";
        usage(argv[0]);
        return 1;
    }
    if (out_damaged_path.empty() && out_undamaged_path.empty()) {
        std::cerr << "Error: at least one of --out-damaged / --out-undamaged required\n";
        usage(argv[0]);
        return 1;
    }

    init_logger("");
    log_info("=== fqdup split: damage-score classification (no deduplication) ===");
    log_info("fqdup version: " FQDUP_VERSION);
    log_info("Input: " + in_path);
    if (!out_damaged_path.empty())   log_info("Output (damaged):   " + out_damaged_path);
    if (!out_undamaged_path.empty()) log_info("Output (undamaged): " + out_undamaged_path);
    if (damage_deam_max_reads > 0)
        log_info("Damage scan: first " + std::to_string(damage_deam_max_reads) + " reads");

    try {
        // --- Pass 0: damage profile estimation ---
        DamageEstimateOptions opts;
        opts.forced_lib         = forced_library_type;
        opts.qc_enabled         = false;
        opts.max_reads          = damage_deam_max_reads;
        opts.threads            = std::max(1u, threads);

        auto t0 = std::chrono::steady_clock::now();
        DamageEstimate est = estimate_damage_with_qc(in_path, opts);
        double t0s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        {
            char b[64]; std::snprintf(b, sizeof(b), "%.1f s", t0s);
            log_info("Phase timer: damage estimation = " + std::string(b));
        }

        // Adaptive Pass-0 escalation. Identifiability of a low-abundance damaged stratum scales with
        // the number of reads scanned, not with d_max: at the sample cap a real LOW_ABUNDANCE stratum
        // (small pi, small d_max) can sit below the floor and the estimator abstains for lack of power.
        // If the capped scan leaves a present terminal signal unresolved, rescan at full depth so the
        // estimator can separate a real stratum from noise instead of dropping it (deep-coverage
        // libraries are exactly where the cap bites). Resolved or genuinely-flat samples skip this.
        if (opts.max_reads > 0) {
            const auto st = est.profile.pi_state;
            const bool resolved = (st == taph::DamageConfidence::DETECTED ||
                                   st == taph::DamageConfidence::LOW_ABUNDANCE);
            const float sig = std::max(est.profile.d_max_5prime, est.profile.d_max_3prime);
            // lsd_hist counts reads that passed the length filter, so it lands below the raw cap even
            // when the cap was hit; use a fraction of the cap as "the cap bound the scan (more reads
            // exist)" rather than ==, which never fires once any read is length-filtered.
            int64_t n_scanned = 0;
            for (uint64_t h : est.lsd_hist) n_scanned += static_cast<int64_t>(h);
            if (!resolved && sig > 0.02f && n_scanned >= opts.max_reads * 4 / 5) {
                char m[160];
                std::snprintf(m, sizeof(m),
                    "Pass 0 inconclusive at %lld reads with d_max=%.3f present; rescanning at full depth.",
                    (long long)n_scanned, sig);
                log_info(m);
                auto t0b = std::chrono::steady_clock::now();
                DamageEstimateOptions full_opts = opts;
                full_opts.max_reads = 0;
                est = estimate_damage_with_qc(in_path, full_opts);
                char b[64]; std::snprintf(b, sizeof(b), "%.1f s",
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0b).count());
                log_info("Phase timer: full-depth rescan = " + std::string(b));
            }
        }

        DamageProfile profile = est.profile;
        log_info("Split model: " +
            std::string(profile.d_max_5prime > 0.01 ? "empirical candidate" : "bulk exponential"));

        // --- Optional LSD scan for empirical per-length-bin model ---
        LengthStratifiedDamageProfile lsd_data;
        const bool run_empirical =
            split_model_mode != SplitModelMode::Bulk &&
            (split_model_mode == SplitModelMode::Empirical ||
             // max(d5,d3): a 3'-only DS signal (5' ->G-overcall artifact-dead, e.g. FLB57md) was
             // gated out of the empirical model by the 5'-only test and fell to the bulk fallback.
             std::max(est.profile.d_max_5prime, est.profile.d_max_3prime) > 0.01);

        if (run_empirical) {
            auto t1 = std::chrono::steady_clock::now();
            const std::vector<uint64_t>* hist_ptr =
                est.lsd_hist.empty() ? nullptr : &est.lsd_hist;
            lsd_data = estimate_damage_split_model(in_path, est.profile, hist_ptr);
            double t1s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
            char b[64]; std::snprintf(b, sizeof(b), "%.1f s", t1s);
            log_info("Phase timer: length-stratified LLR model = " +
                std::string(b) + " (" + std::to_string(lsd_data.bins.size()) + " bins)");
        } else {
            log_info("Split model: bulk exponential (--split-model bulk or no damage detected)");
        }

        DamageSplitModel split_model = DamageSplitModel::build(lsd_data, profile);

        // Posterior-weighted threshold from the VALIDATED reference-free pi (libtaph SplitPolicy):
        // posterior ≥ 0.5 ↔ llr ≥ −log(pi/(1−pi)), so the base cut shifts by the prior log-odds and
        // low-pi libraries aren't flooded with reads that barely cross llr=0. The pi here is the gated
        // estimate (recovery-aware for diluted samples, ss-correct), NOT the old d_max/d_anc proxy.
        // State gates behaviour: ABSTAIN/BELOW_FLOOR (blanks, no identifiable damaged stratum) ⇒ no
        // damaged output — route every read to undamaged rather than manufacture a stratum from noise.
        taph::DamageEstimate pi_est;
        pi_est.point = profile.pi_point;
        pi_est.lo    = profile.pi_lo;
        pi_est.hi    = profile.pi_hi;
        pi_est.state = profile.pi_state;
        const taph::SplitPolicy pol = taph::split_policy(pi_est);

        // Split strategy is AUTO-SELECTED by the estimator's own regime (pi_state):
        //   DETECTED      -> posterior-threshold: cut = base - log(pi/(1-pi)). Reads are individually
        //                    confident, so the realized damaged fraction tracks pi naturally.
        //   LOW_ABUNDANCE -> yield-locked: the recovery regime where pi is small and the per-read signal
        //                    is weak, so a posterior cut mathematically collapses to ~0 (no read clears
        //                    the high prior-driven threshold). Instead pick the cut at the pi-quantile of
        //                    the LLR distribution -> route the top pi*N best-ranked reads as the damaged
        //                    enrichment set. Costs one extra histogram pass over the input.
        //   ABSTAIN/BELOW_FLOOR -> no identifiable stratum: route everything undamaged.
        float effective_threshold = split_threshold;
        if (!pol.splittable) {
            effective_threshold = std::numeric_limits<float>::infinity();
            log_info("Damage pi not identifiable (no damaged stratum): routing all reads to undamaged.");
        } else if (pol.state == taph::DamageConfidence::LOW_ABUNDANCE) {
            // Yield-locked: histogram pass to find the LLR cut whose upper tail mass ≈ pi.
            constexpr int    NB = 20001;
            constexpr double LO = -100.0, HI = 100.0;
            std::vector<uint64_t> hist(NB, 0);
            uint64_t n_scan = 0;
            { auto hr = make_fastq_reader(in_path); FastqRecord hrec;
              while (hr->read(hrec)) {
                  float s = split_model.score(hrec.seq);
                  int bi = static_cast<int>((s - LO) / (HI - LO) * (NB - 1));
                  bi = std::clamp(bi, 0, NB - 1);
                  ++hist[bi]; ++n_scan;
              } }
            const uint64_t target = static_cast<uint64_t>(pol.pi * static_cast<double>(n_scan));
            uint64_t cum = 0; int tb = NB - 1;
            for (; tb > 0; --tb) { cum += hist[tb]; if (cum >= target) break; }
            effective_threshold = static_cast<float>(LO + (HI - LO) * tb / (NB - 1));
            char b[176];
            std::snprintf(b, sizeof(b),
                "%.4f  (yield-locked: top pi=%.4f [%.4f,%.4f] of %llu reads ~= %llu damaged)",
                effective_threshold, pol.pi, pol.pi_lo, pol.pi_hi,
                (unsigned long long)n_scan, (unsigned long long)target);
            log_info("LOW_ABUNDANCE auto-detected -> yield-locked split threshold: " + std::string(b));
            log_info("NOTE: split is pi-calibrated ENRICHMENT, not a pure partition "
                     "(reference-free per-read AUC ~0.59).");
        } else {
            effective_threshold = static_cast<float>(split_threshold - pol.log_prior_odds);
            char b[160];
            std::snprintf(b, sizeof(b),
                "%.4f  (pi=%.4f [%.4f,%.4f]  log_prior_odds=%.4f)",
                effective_threshold, pol.pi, pol.pi_lo, pol.pi_hi, pol.log_prior_odds);
            log_info("DETECTED -> posterior-threshold split: " + std::string(b));
            log_info("NOTE: split is pi-calibrated ENRICHMENT, not a pure partition "
                     "(reference-free per-read AUC ~0.59).");
        }

        // --- Streaming classification pass ---
        unsigned writer_threads = std::min(threads, 16u);
        auto is_gz = [](const std::string& p) {
            return p.size() > 3 && p.substr(p.size()-3) == ".gz";
        };

        std::unique_ptr<FastqWriter> writer_dam, writer_und;
        if (!out_damaged_path.empty())
            writer_dam = std::make_unique<FastqWriter>(
                out_damaged_path, is_gz(out_damaged_path), static_cast<int>(writer_threads));
        if (!out_undamaged_path.empty())
            writer_und = std::make_unique<FastqWriter>(
                out_undamaged_path, is_gz(out_undamaged_path), static_cast<int>(writer_threads));

        auto t2 = std::chrono::steady_clock::now();
        auto reader = make_fastq_reader(in_path);
        FastqRecord rec;
        uint64_t n_total = 0, n_damaged = 0, n_undamaged = 0;

        while (reader->read(rec)) {
            ++n_total;
            float llr = split_model.score(rec.seq);
            if (llr >= effective_threshold) {
                ++n_damaged;
                if (writer_dam) writer_dam->write(rec);
            } else {
                ++n_undamaged;
                if (writer_und) writer_und->write(rec);
            }
        }

        double t2s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t2).count();
        {
            char b[64]; std::snprintf(b, sizeof(b), "%.1f s", t2s);
            log_info("Phase timer: classification pass = " + std::string(b));
        }

        log_info("=== Split complete ===");
        log_info("Total reads: " + std::to_string(n_total));
        log_info("Damaged:     " + std::to_string(n_damaged) +
                 " (" + [&]{ char b[16]; std::snprintf(b,sizeof(b),"%.1f%%",
                     n_total ? 100.0*n_damaged/n_total : 0.0); return std::string(b); }() + ")");
        log_info("Undamaged:   " + std::to_string(n_undamaged) +
                 " (" + [&]{ char b[16]; std::snprintf(b,sizeof(b),"%.1f%%",
                     n_total ? 100.0*n_undamaged/n_total : 0.0); return std::string(b); }() + ")");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    shutdown_logger();
    return 0;
}
