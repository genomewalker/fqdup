#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace fqdup {
// Development primitive only. Independence and calibrated, symmetric
// substitution errors are assumptions; this is not pairing evidence.
struct QualityCompatibility {
    double expected_disagreements;
    double log_tail_bound;
    double log_budget;
    bool excess_disagreement_rejected;
    bool unexpected_agreement_requires_calibration;
};

inline QualityCompatibility quality_compatibility(
    const unsigned char* q1, const unsigned char* q2, std::size_t length,
    std::size_t disagreements, double alpha = 1e-6, std::size_t trials = 1) {
    if (!q1 || !q2 || length == 0 || disagreements > length ||
        !std::isfinite(alpha) || alpha <= 0 || alpha >= 1 || trials == 0)
        throw std::invalid_argument("invalid quality compatibility observations");
    static const std::array<double, 94> error_probability = [] {
        std::array<double, 94> probabilities{};
        for (std::size_t q = 0; q < probabilities.size(); ++q)
            probabilities[q] = std::pow(10.0, -double(q) / 10.0);
        return probabilities;
    }();
    double expected = 0;
    for (std::size_t i = 0; i < length; ++i) {
        if (q1[i] > 93 || q2[i] > 93)
            throw std::invalid_argument("invalid Phred quality");
        const double e1 = error_probability[q1[i]];
        const double e2 = error_probability[q2[i]];
        expected += e1 + e2 - (4.0 / 3.0) * e1 * e2;
    }
    const double count = double(disagreements);
    const double rate = disagreements == 0 ? expected :
        count * std::log(count / expected) - count + expected;
    // Separate logarithms avoid underflow when splitting a very small budget.
    const double budget = std::log(alpha) - std::log(double(trials));
    return {expected, -rate, budget,
            count > expected && -rate <= budget,
            count < expected && -rate <= budget};
}
inline bool adaptive_alignment_supported(
    const char* sequence1, const char* quality1, const char* sequence2,
    const char* quality2, int length, int minimum, std::size_t trials,
    bool require_sequence_evidence) {
    std::vector<unsigned char> one, two;
    int mismatches = 0;
    auto canonical = [](char b) { return b == 'A' || b == 'C' || b == 'G' || b == 'T'; };
    for (int i = 0; i < length; ++i) {
        if (!canonical(sequence1[i]) || !canonical(sequence2[i])) continue;
        const int q1 = static_cast<unsigned char>(quality1[i]) - 33;
        const int q2 = quality2 ? static_cast<unsigned char>(quality2[i]) - 33 : 93;
        if (q1 < 0 || q1 > 93 || q2 < 0 || q2 > 93) return false;
        one.push_back(static_cast<unsigned char>(q1));
        two.push_back(static_cast<unsigned char>(q2));
        mismatches += sequence1[i] != sequence2[i];
    }
    if (one.empty() || static_cast<int>(one.size()) < minimum) return false;
    const auto fit = quality_compatibility(one.data(), two.data(), one.size(), mismatches, 1e-6, trials);
    if (fit.excess_disagreement_rejected) return false;
    // Unexpected perfect agreement diagnoses the quality model, not a wrong
    // genotype. Exact canonical observations may establish geometry using the
    // separate sequence floor and unique-candidate checks. This never repairs
    // disagreement or makes the resulting Phred model calibrated.
    if (fit.unexpected_agreement_requires_calibration && mismatches != 0) return false;
    if (!require_sequence_evidence) return true;
    // Use the two observed compositions as a conservative order-zero
    // background instead of assuming uniform DNA. This is a development
    // evidence floor, not a calibrated genomic false-discovery probability.
    auto code = [](char b) { return b == 'A' ? 0 : b == 'C' ? 1 : b == 'G' ? 2 : 3; };
    std::array<std::size_t, 4> bases1{}, bases2{};
    for (int i = 0; i < length; ++i) {
        if (!canonical(sequence1[i]) || !canonical(sequence2[i])) continue;
        ++bases1[code(sequence1[i])];
        ++bases2[code(sequence2[i])];
    }
    double background = 0;
    for (int b = 0; b < 4; ++b)
        background += double(bases1[b]) * bases2[b];
    background = std::max(.25, background / (double(one.size()) * one.size()));
    const double required = std::log(1e6) + std::log(double(trials));
    auto supported = [&](std::size_t count, std::size_t matches, double null_match) {
        if (!count || null_match >= 1) return false;
        const double fraction = double(matches) / count;
        if (fraction <= null_match) return false;
        double divergence = fraction * std::log(fraction / null_match);
        if (fraction < 1)
            divergence += (1 - fraction) * std::log((1 - fraction) / (1 - null_match));
        return count * divergence >= required;
    };
    if (!supported(one.size(), one.size() - mismatches, background)) return false;

    // Nonoverlapping four-base words flag highly repetitive backgrounds whose
    // single-base composition looks diverse (for example ACACAC...). Noisy
    // diverse overlaps keep the base-composition test above; requiring both
    // floors there discards true molecules after several sequencing errors.
    std::array<std::size_t, 256> words1{}, words2{};
    std::size_t blocks = 0, matching_blocks = 0;
    for (int i = 0; i + 4 <= length; i += 4) {
        int word1 = 0, word2 = 0;
        bool valid = true;
        for (int j = 0; j < 4; ++j) {
            if (!canonical(sequence1[i + j]) || !canonical(sequence2[i + j])) {
                valid = false;
                break;
            }
            word1 = 4 * word1 + code(sequence1[i + j]);
            word2 = 4 * word2 + code(sequence2[i + j]);
        }
        if (!valid) continue;
        ++words1[word1]; ++words2[word2];
        ++blocks;
        matching_blocks += word1 == word2;
    }
    if (!blocks) return false;
    double word_background = 0;
    for (std::size_t i = 0; i < words1.size(); ++i)
        word_background += double(words1[i]) * words2[i];
    word_background /= double(blocks) * blocks;
    return word_background < .25 ||
           supported(blocks, matching_blocks, word_background);
}
} // namespace fqdup

