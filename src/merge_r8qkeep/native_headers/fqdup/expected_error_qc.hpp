#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string_view>

namespace fqdup {
// Expectations are additive without independent-base assumptions. The Phred
// marginal error model must nevertheless be calibrated for biological claims.
inline std::size_t longest_error_budget_prefix(std::string_view sequence,
        std::string_view quality, double max_expected_errors = -1,
        double max_error_rate = -1) {
    if (sequence.size() != quality.size() || !std::isfinite(max_expected_errors) ||
        !std::isfinite(max_error_rate) || (max_expected_errors < 0 && max_expected_errors != -1) ||
        (max_error_rate < 0 && max_error_rate != -1) || max_error_rate > 1)
        throw std::invalid_argument("invalid expected-error QC inputs");
    if (max_expected_errors < 0 && max_error_rate < 0) return sequence.size();
    static const std::array<long double, 94> errors = [] {
        std::array<long double, 94> result{};
        for (std::size_t q = 0; q < result.size(); ++q)
            result[q] = std::pow(10.0L, -static_cast<long double>(q) / 10.0L);
        return result;
    }();
    long double expected = 0;
    std::size_t keep = 0;
    auto within = [](long double observed, long double limit) {
        const long double rounding = 32 * std::numeric_limits<double>::epsilon() *
                                     std::max(observed, limit);
        return observed <= limit + rounding;
    };
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const char base = sequence[i];
        const int q = static_cast<unsigned char>(quality[i]) - 33;
        // Unknown bases and invalid quality observations carry maximal error.
        expected += (base != 'A' && base != 'C' && base != 'G' && base != 'T') ||
                    q < 0 || q > 93 ? 1.0L : errors[static_cast<std::size_t>(q)];
        const bool total_passes = max_expected_errors < 0 || within(expected, max_expected_errors);
        const bool rate_passes = max_error_rate < 0 ||
            within(expected, static_cast<long double>(max_error_rate) * (i + 1));
        if (total_passes && rate_passes) keep = i + 1;
        // A later high-quality segment can restore the mean-error criterion;
        // do not stop at an earlier violation.
    }
    return keep;
}
} // namespace fqdup

