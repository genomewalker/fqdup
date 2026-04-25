#pragma once
// Per-cluster index/arena types and Phase 3 stat counters used by derep.cpp.
// Internal to derep.cpp.

#include "encoding.hpp"
#include "fqdup/logger.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fqdup::derep_detail {

struct ErrCorParams {
    bool     enabled            = false;
    double   snp_threshold      = 0.20;
    uint32_t snp_min_count      = 1;
    uint32_t snp_low_cov_cutoff = 10;
    double   snp_low_cov_factor = 1.75;
    uint32_t bucket_cap         = 0;
    double   pcr_phi            = 5.3e-7;
    double   pcr_rate           = 0.0;
    uint32_t max_h2_count       = 2;
    unsigned threads            = 1;  // 0 = hardware_concurrency
    // T5.8: empirical posterior-odds decision rule. The model fits per-bin
    // P(real | subst×term×occ×damage) and per-occ-bin log(π_pcr/π_real)
    // from B1 candidate edges, then absorbs iff posterior log-odds S > 0.
    // No calibrated threshold — the rule is the definition of posterior odds.
    // Set legacy_veto=true to bypass the empirical model and use the older
    // SNP-veto path instead (kept for ablation studies).
    bool     empirical          = true;
    bool     legacy_veto        = false;  // legacy SNP-only path (no LR/posterior)
    bool     adj_len_probe      = false;  // T5.6: adjacent-length (L±1) probe
};

struct Phase3Stats {
    double decode_hash_parent_ms = 0;
    double insert_ms             = 0;
    double decode_hash_child_ms  = 0;
    double query_ms              = 0;
    double check_ms              = 0;
    uint64_t total_candidates    = 0;
    uint64_t bucket_overflow_drops = 0;
    uint64_t children_scanned    = 0;
    uint64_t parents_indexed     = 0;
    uint64_t children_found      = 0;
    uint64_t snp_protected       = 0;
    uint64_t absorbed            = 0;
    uint64_t short_brute_parents   = 0;
    uint64_t short_brute_evaluated = 0;
    uint64_t short_brute_found     = 0;
    uint64_t short_too_small_skipped = 0;
    // Shadow LR stats — accumulated for every absorption decision (T5.1).
    // Used to calibrate the LR threshold (T5.4) before flipping the decision rule.
    double   lr_sum_absorbed = 0.0;
    uint64_t lr_n_absorbed   = 0;
    double   lr_sum_protected = 0.0;
    uint64_t lr_n_protected  = 0;
    // T5.5: edge-reason classification. Set in B2.
    uint64_t edge_lr_absorbed   = 0;  // absorbed because LR > threshold
    uint64_t edge_lr_protected  = 0;  // protected because LR <= threshold
    uint64_t edge_legacy_absorb = 0;  // absorbed under legacy SNP-veto path
    uint64_t edge_legacy_veto   = 0;  // protected under legacy SNP-veto path
    uint64_t edge_damage_bypass = 0;  // damage-channel bypass fired
    // T5.6 adjacent-length probe stats.
    uint64_t adj_len_evaluated = 0;
    uint64_t adj_len_matched   = 0;  // banded check passed (1-indel, 0-sub)
    uint64_t adj_len_absorbed  = 0;
    uint64_t adj_len_protected = 0;  // matched but LR ≤ threshold (or legacy mode)
    void log() const {
        auto f = [](double ms){ return std::to_string(static_cast<int>(ms)) + " ms"; };
        log_info("Phase 3 timing:");
        log_info("  Parent decode+hash : " + f(decode_hash_parent_ms));
        log_info("  Parent insert      : " + f(insert_ms));
        log_info("  Child decode+hash  : " + f(decode_hash_child_ms));
        log_info("  Child query        : " + f(query_ms));
        log_info("  Candidate check    : " + f(check_ms));
        log_info("  Parents indexed    : " + std::to_string(parents_indexed));
        log_info("  Children scanned   : " + std::to_string(children_scanned));
        log_info("  Children found     : " + std::to_string(children_found));
        log_info("  SNP protected      : " + std::to_string(snp_protected));
        log_info("  Absorbed           : " + std::to_string(absorbed));
        log_info("  Short brute parents: " + std::to_string(short_brute_parents));
        log_info("  Short brute evaluated: " + std::to_string(short_brute_evaluated));
        log_info("  Short brute found  : " + std::to_string(short_brute_found));
        log_info("  Short too small    : " + std::to_string(short_too_small_skipped));
        log_info("  Total candidates   : " + std::to_string(total_candidates));
        if (children_scanned > 0)
            log_info("  Avg cand/child     : " +
                     std::to_string(static_cast<double>(total_candidates)/children_scanned));
        log_info("  Bucket overflow drops: " + std::to_string(bucket_overflow_drops));
        if (lr_n_absorbed + lr_n_protected > 0) {
            auto fmt = [](double v){ char b[32]; std::snprintf(b,sizeof(b),"%.3f",v); return std::string(b); };
            double mean_abs = lr_n_absorbed  ? lr_sum_absorbed  / lr_n_absorbed  : 0.0;
            double mean_pro = lr_n_protected ? lr_sum_protected / lr_n_protected : 0.0;
            log_info("  Shadow LR (absorbed) : mean=" + fmt(mean_abs) + "  n=" + std::to_string(lr_n_absorbed));
            log_info("  Shadow LR (protected): mean=" + fmt(mean_pro) + "  n=" + std::to_string(lr_n_protected));
        }
        if (edge_lr_absorbed + edge_lr_protected + edge_legacy_absorb + edge_legacy_veto + edge_damage_bypass > 0) {
            log_info("Phase 3 edge reasons:");
            log_info("  LR absorbed       : " + std::to_string(edge_lr_absorbed));
            log_info("  LR protected      : " + std::to_string(edge_lr_protected));
            log_info("  Legacy absorbed   : " + std::to_string(edge_legacy_absorb));
            log_info("  Legacy SNP veto   : " + std::to_string(edge_legacy_veto));
            log_info("  Damage bypass     : " + std::to_string(edge_damage_bypass));
        }
        if (adj_len_evaluated > 0) {
            log_info("Phase 3 adjacent-length probe (T5.6):");
            log_info("  Pairs evaluated   : " + std::to_string(adj_len_evaluated));
            log_info("  1-indel matched   : " + std::to_string(adj_len_matched));
            log_info("  Absorbed          : " + std::to_string(adj_len_absorbed));
            log_info("  Protected (LR/legacy): " + std::to_string(adj_len_protected));
        }
        if (bucket_overflow_drops > 0)
            log_warn("Phase 3: " + std::to_string(bucket_overflow_drops) +
                     " pair-key entries dropped — some PCR error candidates were not evaluated. "
                     "Set --errcor-bucket-cap 0 (default) to disable the cap entirely.");
    }
};

struct IndexEntry {
    uint64_t record_index;
    uint64_t count;
    uint32_t seq_id;
    uint8_t  damage_score;

    IndexEntry() : record_index(0), count(1), seq_id(0), damage_score(0) {}
    explicit IndexEntry(uint64_t idx) : record_index(idx), count(1), seq_id(0), damage_score(0) {}
};

// 2-bit packed sequence arena. ~4× memory reduction vs ASCII.
// Non-ACGT bases are encoded as A=0 and the sequence is flagged ineligible
// so Phase 3 skips it.
struct SeqArena {
    std::vector<uint8_t>  packed;
    std::vector<uint64_t> offsets;
    std::vector<uint16_t> lengths;
    std::vector<bool>     eligible;

    uint32_t append(const std::string& seq) {
        if (seq.size() > 65535u)
            throw std::runtime_error("Sequence too long for arena (>65535 bp)");
        if (offsets.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("SeqArena overflow: more than 4 G unique sequences");
        uint32_t id    = static_cast<uint32_t>(offsets.size());
        uint16_t L     = static_cast<uint16_t>(seq.size());
        int      nbytes = (L + 3) / 4;
        offsets.push_back(packed.size());
        lengths.push_back(L);
        packed.resize(packed.size() + nbytes, 0);
        uint8_t* dst = packed.data() + offsets[id];
        bool ok = true;
        int i = 0, byte = 0;
        while (i < L) {
            uint8_t v = 0;
            for (int lane = 0; lane < 4 && i < L; ++lane, ++i) {
                uint8_t b = kEnc2bit[static_cast<uint8_t>(seq[i])];
                if (b == 0xFFu) { ok = false; b = 0; }
                v |= static_cast<uint8_t>(b << (6 - 2 * lane));
            }
            dst[byte++] = v;
        }
        eligible.push_back(ok);
        return id;
    }

    uint32_t append_chars(const char* data_in, int L) {
        if (L > 65535)
            throw std::runtime_error("Sequence too long for arena (>65535 bp)");
        if (offsets.size() >= static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("SeqArena overflow: more than 4 G unique sequences");
        uint32_t id    = static_cast<uint32_t>(offsets.size());
        int      nbytes = (L + 3) / 4;
        offsets.push_back(packed.size());
        lengths.push_back(static_cast<uint16_t>(L));
        packed.resize(packed.size() + nbytes, 0);
        uint8_t* dst = packed.data() + offsets[id];
        bool ok = true;
        int i = 0, byte = 0;
        while (i < L) {
            uint8_t v = 0;
            for (int lane = 0; lane < 4 && i < L; ++lane, ++i) {
                uint8_t b = kEnc2bit[static_cast<uint8_t>(data_in[i])];
                if (b == 0xFFu) { ok = false; b = 0; }
                v |= static_cast<uint8_t>(b << (6 - 2 * lane));
            }
            dst[byte++] = v;
        }
        eligible.push_back(ok);
        return id;
    }

    const uint8_t* data(uint32_t id)        const { return packed.data() + offsets[id]; }
    uint16_t       length(uint32_t id)      const { return lengths[id]; }
    bool           is_eligible(uint32_t id) const { return eligible[id]; }
    size_t         size()                   const { return offsets.size(); }

    void decode_range(uint32_t id, int start, int count, uint8_t* out) const {
        const uint8_t* p = data(id);
        int i = 0;
        int pos = start;

        while (i < count && (pos & 3) != 0) {
            out[i++] = kDec2bit[(p[pos >> 2] >> (6 - 2 * (pos & 3))) & 0x3u];
            ++pos;
        }

        while (i + 4 <= count) {
            const auto& d = kDec4[p[pos >> 2]];
            std::memcpy(out + i, d.data(), 4);
            i += 4;
            pos += 4;
        }

        while (i < count) {
            out[i++] = kDec2bit[(p[pos >> 2] >> (6 - 2 * (pos & 3))) & 0x3u];
            ++pos;
        }
    }
};

// Raw-byte quality arena, parallel to SeqArena. Stores Phred+33 ASCII
// (caller already has it in that form when reading fastq). Phase 3 LR
// scoring reads qualities only at mismatch positions, so the random-access
// cost dominates pack/unpack savings — keep it simple and dense.
//
// Indexed by the same uint32_t id returned by SeqArena.append(). Append
// must occur in lockstep with the seq arena.
struct QualArena {
    std::vector<uint8_t>  bytes;
    std::vector<uint64_t> offsets;
    std::vector<uint16_t> lengths;

    uint32_t append(const std::string& qual) {
        if (qual.size() > 65535u)
            throw std::runtime_error("Quality too long for arena (>65535 bp)");
        uint32_t id = static_cast<uint32_t>(offsets.size());
        offsets.push_back(bytes.size());
        lengths.push_back(static_cast<uint16_t>(qual.size()));
        bytes.insert(bytes.end(),
                     reinterpret_cast<const uint8_t*>(qual.data()),
                     reinterpret_cast<const uint8_t*>(qual.data()) + qual.size());
        return id;
    }

    uint32_t append_chars(const char* data_in, int L) {
        if (L > 65535)
            throw std::runtime_error("Quality too long for arena (>65535 bp)");
        uint32_t id = static_cast<uint32_t>(offsets.size());
        offsets.push_back(bytes.size());
        lengths.push_back(static_cast<uint16_t>(L));
        bytes.insert(bytes.end(),
                     reinterpret_cast<const uint8_t*>(data_in),
                     reinterpret_cast<const uint8_t*>(data_in) + L);
        return id;
    }

    // Push a placeholder when the caller has no qual (e.g. .fa input).
    // Returns id; q_at() will yield kNoQual.
    uint32_t append_empty() {
        uint32_t id = static_cast<uint32_t>(offsets.size());
        offsets.push_back(bytes.size());
        lengths.push_back(0);
        return id;
    }

    static constexpr uint8_t kNoQual = 0;  // Phred 0 → P(err)=1.0; LR will fall back to neutral

    // Phred+33 → integer Q at a single position. Bounds-checked: out-of-range
    // returns kNoQual rather than UB so Phase 3 short-interior fallbacks are safe.
    uint8_t q_at(uint32_t id, int pos) const {
        if (id >= lengths.size()) return kNoQual;
        if (lengths[id] == 0 || pos < 0 || pos >= lengths[id]) return kNoQual;
        uint8_t v = bytes[offsets[id] + pos];
        return v >= 33 ? static_cast<uint8_t>(v - 33) : kNoQual;
    }

    uint16_t length(uint32_t id) const { return id < lengths.size() ? lengths[id] : 0; }
    size_t   size()              const { return offsets.size(); }
};

struct FlatPairIndex {
    std::vector<uint64_t> keys;
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> ids;

    template <typename F>
    void query(uint64_t key, F&& fn) const {
        auto it = std::lower_bound(keys.begin(), keys.end(), key);
        if (it == keys.end() || *it != key) return;
        size_t idx = static_cast<size_t>(it - keys.begin());
        for (uint32_t i = offsets[idx]; i < offsets[idx + 1]; ++i)
            fn(ids[i]);
    }
};

}  // namespace fqdup::derep_detail
