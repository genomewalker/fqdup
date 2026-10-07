#pragma once
#include "fqdup/quality_compatibility.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace fqdup {

// Development primitive for separately supported adapter clipping only.
// A single contiguous indel may occur in either observed mate prefix; this
// never creates a gapped consensus or infers the inserted base's origin.
// The caller must account for all anchor-pair and gap-placement hypotheses in
// trials and must reject competing adapter-boundary pairs separately.
inline bool one_run_gap_supported(std::string_view first, std::string_view first_quality,
                                  std::string_view second, std::string_view second_quality,
                                  std::size_t minimum, std::size_t trials,
                                  std::size_t minimum_flank = 1,
                                  std::size_t maximum_mismatches = static_cast<std::size_t>(-1),
                                  std::size_t minimum_right_flank = 1,
                                  std::size_t maximum_gap_start = static_cast<std::size_t>(-1)) {
    if (first.size() != first_quality.size() || second.size() != second_quality.size() ||
        first.size() == second.size() || minimum == 0 || trials == 0)
        return false;
    const auto longer = first.size() > second.size() ? first : second;
    const auto longer_quality = first.size() > second.size() ? first_quality : second_quality;
    const auto shorter = first.size() > second.size() ? second : first;
    const auto shorter_quality = first.size() > second.size() ? second_quality : first_quality;
    if (shorter.size() < minimum) return false;
    auto canonical = [](char b) { return b == 'A' || b == 'C' || b == 'G' || b == 'T'; };
    if (!std::all_of(first.begin(), first.end(), canonical) ||
        !std::all_of(second.begin(), second.end(), canonical)) return false;
    const auto gap = longer.size() - shorter.size();
    // Sequence evidence later uses an order-zero background of at least 1/4.
    // Its divergence against 1/4 bounds the possible evidence from above, so
    // placements failing this inexpensive check cannot pass the full model.
    const double required = std::log(1e6) + std::log(double(trials));
    auto can_pass_sequence_floor = [&](std::size_t mismatches) {
        const double count = double(shorter.size());
        const double fraction = (count - double(mismatches)) / count;
        if (fraction <= .25) return false;
        double divergence = fraction * std::log(fraction / .25);
        if (fraction < 1)
            divergence += (1 - fraction) * std::log((1 - fraction) / .75);
        return count * divergence + 1e-10 >= required;
    };
    // For a gap after `left` bases, comparisons before it use longer[i],
    // while comparisons after it use longer[i + gap]. Count both parts once
    // so placements rejected by the sequence floor allocate no strings.
    const std::size_t n = shorter.size();
    std::vector<std::size_t> prefix(n + 1), suffix(n + 1);
    for (std::size_t i = 0; i < n; ++i)
        prefix[i + 1] = prefix[i] + (longer[i] != shorter[i]);
    for (std::size_t i = n; i-- > 0;)
        suffix[i] = suffix[i + 1] + (longer[i + gap] != shorter[i]);
    // At least one aligned base on each side establishes an internal gap;
    // terminal overhangs do not locate a biological insertion/deletion.
    for (std::size_t left = minimum_flank;
         left + minimum_right_flank <= shorter.size() && left <= maximum_gap_start; ++left) {
        const std::size_t mismatches = prefix[left] + suffix[left];
        if (mismatches > maximum_mismatches ||
            !can_pass_sequence_floor(mismatches)) continue;
        std::string aligned, quality;
        aligned.reserve(n); quality.reserve(n);
        aligned.append(longer.substr(0, left));
        aligned.append(longer.substr(left + gap));
        quality.append(longer_quality.substr(0, left));
        quality.append(longer_quality.substr(left + gap));
        if (adaptive_alignment_supported(aligned.data(), quality.data(),
                                         shorter.data(), shorter_quality.data(),
                                         static_cast<int>(shorter.size()),
                                         static_cast<int>(minimum), trials, true))
            return true;
    }
    return false;
}

// A full-length biological pair can contain one sequencing indel and a
// compensating terminal truncation when both reads have the same cycle count.
// Such an explanation defeats an internal adapter motif or an ungapped merge.
// Bounded gap search is a conservative veto, not a repaired consensus.
inline bool full_read_gap_alternative(std::string_view first, std::string_view first_quality,
                                      std::string_view second, std::string_view second_quality,
                                      std::size_t minimum, std::size_t trials,
                                      std::size_t maximum_gap_start,
                                      std::size_t max_gap = 8) {
    if (first.size() != second.size() || first.size() != first_quality.size() ||
        second.size() != second_quality.size() || first.size() <= minimum ||
        minimum == 0 || trials == 0)
        return false;
    const std::size_t upper = std::min(max_gap, first.size() - minimum);
    const auto canonical = [](char b) { return b == 'A' || b == 'C' || b == 'G' || b == 'T'; };
    if (!std::all_of(first.begin(), first.end(), canonical) ||
        !std::all_of(second.begin(), second.end(), canonical)) return false;
    const double required = std::log(1e6) + std::log(double(trials));
    // The three original flank/mismatch policies differ only in their
    // admissible placement. Search their union once per gap and orientation.
    auto orientation = [&](std::string_view longer, std::string_view longer_quality,
                           std::string_view shorter, std::string_view shorter_quality,
                           std::size_t gap) {
        const std::size_t n = shorter.size();
        if (n < minimum) return false;
        std::vector<std::size_t> prefix(n + 1), suffix(n + 1);
        for (std::size_t i = 0; i < n; ++i)
            prefix[i + 1] = prefix[i] + (longer[i] != shorter[i]);
        for (std::size_t i = n; i-- > 0;)
            suffix[i] = suffix[i + 1] + (longer[i + gap] != shorter[i]);
        for (std::size_t left = 1; left < n && left <= maximum_gap_start; ++left) {
            const std::size_t mismatches = prefix[left] + suffix[left];
            const std::size_t right = n - left;
            if (!((left >= 4 && right >= 4 && mismatches == 0) ||
                  (left >= 4 && right >= minimum && mismatches <= 1) ||
                  (left >= minimum && right >= minimum))) continue;
            const double fraction = (double(n) - double(mismatches)) / double(n);
            if (fraction <= .25) continue;
            double divergence = fraction * std::log(fraction / .25);
            if (fraction < 1)
                divergence += (1 - fraction) * std::log((1 - fraction) / .75);
            if (double(n) * divergence + 1e-10 < required) continue;
            std::string aligned, quality;
            aligned.reserve(n); quality.reserve(n);
            aligned.append(longer.substr(0, left));
            aligned.append(longer.substr(left + gap));
            quality.append(longer_quality.substr(0, left));
            quality.append(longer_quality.substr(left + gap));
            if (adaptive_alignment_supported(aligned.data(), quality.data(),
                                             shorter.data(), shorter_quality.data(),
                                             static_cast<int>(n), static_cast<int>(minimum),
                                             trials, true)) return true;
        }
        return false;
    };
    for (std::size_t gap = 1; gap <= upper; ++gap) {
        const std::size_t n = first.size() - gap;
        if (orientation(first, first_quality, second.substr(0, n),
                        second_quality.substr(0, n), gap) ||
            orientation(second, second_quality, first.substr(0, n),
                        first_quality.substr(0, n), gap)) return true;
    }
    return false;
}

}  // namespace fqdup

