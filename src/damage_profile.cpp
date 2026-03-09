// damage_profile.cpp
// estimate_damage() — DART full pipeline damage estimation.
// Called by both 'derep' (Pass 0) and 'extend' (Pass 0).

#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"

#include "dart/frame_selector_decl.hpp"
// dart/sample_damage_profile.hpp is included transitively via damage_profile.hpp

#include <string>

static constexpr double MIN_COV_DP = 100.0;

static DamageProfile estimate_damage_impl(
    FastqReaderBase& reader,
    const std::string& path,
    double mask_threshold,
    dart::SampleDamageProfile::LibraryType forced_lib,
    int64_t max_reads)
{
    FastqRecord rec;
    dart::SampleDamageProfile dart_profile;
    dart_profile.forced_library_type = forced_lib;

    int      reads_scanned = 0;
    int      typical_len   = 0;
    int64_t  record_pos    = 0;

    while (reader.read(rec)) {
        record_pos++;
        if (max_reads > 0 && record_pos > max_reads) break;
        int L = static_cast<int>(rec.seq.size());
        if (L < 30) continue;
        if (L > typical_len) typical_len = L;
        dart::FrameSelector::update_sample_profile(dart_profile, rec.seq);
        reads_scanned++;
    }

    dart::FrameSelector::finalize_sample_profile(dart_profile);

    double d_max_5   = dart_profile.d_max_5prime;
    double d_max_3   = dart_profile.d_max_3prime;
    double lambda_5  = dart_profile.lambda_5prime;
    double lambda_3  = dart_profile.lambda_3prime;
    double bg_5      = dart_profile.fit_baseline_5prime;
    double bg_3      = dart_profile.fit_baseline_3prime;
    double background = (bg_5 + bg_3) / 2.0;
    double d_max_combined = dart_profile.d_max_combined;

    const bool is_ss = (dart_profile.library_type ==
                        dart::SampleDamageProfile::LibraryType::SINGLE_STRANDED);

    DamageProfile profile;
    profile.d_max_5prime   = d_max_5;
    profile.d_max_3prime   = d_max_3;
    profile.lambda_5prime  = lambda_5;
    profile.lambda_3prime  = lambda_3;
    profile.background     = background;
    profile.mask_threshold = mask_threshold;
    profile.ss_mode        = is_ss;
    profile.enabled        = (d_max_combined > 0.02 || d_max_5 > 0.02 || d_max_3 > 0.02
                              || dart_profile.damage_validated);

    // Empirical per-position mask from DART's normalized frequencies.
    // DS: max(5' C→T excess, 3' G→A excess) > threshold.
    // SS: max(5' C→T excess, 3' T/(T+C) excess) > threshold — both ends have C→T.
    double bg_tc = dart_profile.baseline_t_freq /
                   (dart_profile.baseline_t_freq + dart_profile.baseline_c_freq + 1e-9);
    for (int p = 0; p < DamageProfile::MASK_POSITIONS; ++p) {
        double excess_5 = 0.0, excess_3 = 0.0;
        if (dart_profile.tc_total_5prime[p] >= MIN_COV_DP)
            excess_5 = dart_profile.t_freq_5prime[p] - bg_5;
        if (is_ss) {
            if (dart_profile.tc_total_3prime[p] >= MIN_COV_DP) {
                double tc3 = dart_profile.tc_total_3prime[p];
                double ct3_freq = (tc3 > 0)
                    ? dart_profile.t_freq_3prime[p] / tc3
                    : bg_tc;
                excess_3 = ct3_freq - bg_tc;
            }
        } else {
            if (dart_profile.ag_total_3prime[p] >= MIN_COV_DP)
                excess_3 = dart_profile.a_freq_3prime[p] - bg_3;
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
    log_info("Pass 0: DART damage estimation — " +
             std::to_string(reads_scanned) + " reads in " + path + sample_note);
    log_info("  Library type: " + std::string(dart_profile.library_type_str()) +
             (dart_profile.library_type_auto_detected ? " (auto-detected)" : " (forced)"));
    log_info("  5'-end: d_max=" + std::to_string(d_max_5) +
             " lambda=" + std::to_string(lambda_5) +
             " bg=" + std::to_string(bg_5));
    log_info("  3'-end: d_max=" + std::to_string(d_max_3) +
             " lambda=" + std::to_string(lambda_3) +
             " bg=" + std::to_string(bg_3));
    log_info("  d_max_combined=" + std::to_string(d_max_combined) +
             " (source=" + dart_profile.d_max_source_str() + ")" +
             (dart_profile.mixture_converged ? " [mixture]" : "") +
             (dart_profile.damage_validated  ? " [validated]" : "") +
             (dart_profile.damage_artifact   ? " [ARTIFACT]" : ""));
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
    dart::SampleDamageProfile::LibraryType forced_lib,
    int64_t max_reads)
{
    auto reader = make_fastq_reader(path);
    return estimate_damage_impl(*reader, path, mask_threshold, forced_lib, max_reads);
}
