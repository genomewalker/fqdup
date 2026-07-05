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
    // Singleton high-Q SNP veto: block H=1/H=2 absorption when child carries a
    // mismatch base with quality >= this and the mismatch is not damage-bypass.
    // Rationale: empirical S>0 still permits absorption of singleton (sig=1)
    // SNP-bearing reads at low-coverage variant loci. A high-Q non-damage base
    // is direct sequence-divergence evidence; blocking costs at most one dup.
    uint8_t  singleton_qual_min = 25;
    uint32_t snp_low_cov_cutoff = 10;
    double   snp_low_cov_factor = 1.75;
    uint32_t bucket_cap         = 0;
    double   pcr_phi            = 5.3e-7;
    double   pcr_rate           = 0.0;
    uint32_t max_h2_count       = 2;
    unsigned threads            = 0;  // 0 = hardware_concurrency (default)
    // T5.8: empirical posterior-odds decision rule. The model fits per-bin
    // P(real | subst×term×occ×damage) and per-occ-bin log(π_pcr/π_real)
    // from B1 candidate edges, then absorbs iff posterior log-odds S > 0.
    // No calibrated threshold — the rule is the definition of posterior odds.
    // Set legacy_veto=true to bypass the empirical model and use the older
    // SNP-veto path instead (kept for ablation studies).
    bool     empirical          = true;
    bool     legacy_veto        = false;  // legacy SNP-only path (no LR/posterior)
    // B1 damage-adjusted Hamming: terminal C→T (5') and G→A (3', DS only) mismatches
    // are not counted toward the H≤2 admission threshold. h_adj=1 pairs are recorded
    // as H=1 (non-damage mismatch only); h_adj=0 pairs (both terminal damage) are
    // skipped here and require tag-folding to reach Phase 3.
    // Set false for ablation to restore old transversion-only H=2 behaviour.
    bool     b1_damage_adjust       = true;
    // Protect Channel G (C↔G) and Channel H (A↔T) substitutions from H=1/H=2
    // PCR-error absorption. Off by default; enable with --protect-transversions
    // for high-oxidative-damage libraries.
    bool     protect_transversions  = false;
    // Calibration-free LRT for interior C↔T / G↔A candidates admitted by
    // b1_damage_adjust. Bypasses the empirical model for this subclass.
    // f0 = expected child fraction under PCR error (cycle 3-4 ≈ 0.10),
    // f1 = under real variant (0.50). Absorb if log(P(H0)/P(H1)) > log(lrt_T).
    // Hard absorb when parent/child > lrt_hard_ratio regardless of LRT.
    // SNP veto uses Wilson-95 lower CI bound instead of raw frequency.
    // κ gate: ratio of interior C↔T mismatches to interior transversions.
    // κ≈1 → no excess above PCR floor → suppress absorption.
    // κ≥b1_kappa_min → genuine interior deamination signal → allow absorption.
    // Cross-strand veto: child cluster with reads from both orientations (≥cs_min_total total)
    // is treated as a real SNP (strand-symmetric) and protected from absorption.
    double   b1_kappa_min        = 2.0;
    uint32_t b1_cs_min_total     = 3;
    double   lrt_f0              = 0.10;
    double   lrt_f1              = 0.50;
    double   lrt_T               = 10.0;
    double   lrt_hard_ratio      = 30.0;
    double   snp_cp_lb_threshold = 0.10;
    bool     adj_len_probe      = false;  // T5.6: adjacent-length (L±1) probe
    // T8: opt-in indel-rescue path (syncmer-indexed, banded ed≤2).
    bool     rescue_indels      = false;
    int32_t  rescue_min_hits    = 0;      // 0 = auto (3 for cap=16, 2 for cap=8)
    uint32_t rescue_topk        = 32;
    uint32_t rescue_bundle_hot  = 50;
    uint32_t rescue_hash_hot    = 4096;
    double   rescue_alpha_ins   = 2.0;
    double   rescue_alpha_del   = 2.0;
    double   rescue_mask_bonus  = 1.0;
    // Phase B3: damage-aware H>2 merge for heavily damaged ancient DNA.
    // Pairs that differ only in deamination at near-terminal positions
    // (p_damage > b3_deam_threshold) are merged when 3 ≤ H ≤ b3_max_hamming
    // and the count LRT confirms PCR-error origin.
    bool     b3_enabled             = true;
    // Gate: skip B3 if fewer than b3_min_n_elig interior positions have
    // P(damage) > b3_deam_threshold, OR if their sum (expected deamination
    // events in unmasked interior) < b3_min_mass. Both are derived from the
    // fitted profile — no magic d_max cutoff needed.
    int      b3_min_n_elig          = 3;     // min eligible positions (= H_min)
    float    b3_min_mass            = 0.1f;  // min sum(P(damage)) over eligible positions
    float    b3_deam_threshold      = 0.01f; // P(damage) threshold for normalization
    int      b3_max_hamming         = 5;     // max H to consider in B3
    double   b3_count_ratio         = 5.0;   // min parent/child count ratio
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
    // Per-occupancy-bin absorption breakdown (caveat-check: are absorptions
    // concentrated in high-occ bins where the prior favors PCR?).
    std::array<uint64_t, 6> edge_lr_absorbed_by_occ{};
    std::array<uint64_t, 6> edge_lr_protected_by_occ{};
    uint64_t edge_legacy_absorb = 0;  // absorbed under legacy SNP-veto path
    uint64_t edge_legacy_veto   = 0;  // protected under legacy SNP-veto path
    uint64_t edge_damage_bypass = 0;  // damage-channel bypass fired
    // T5.6 adjacent-length probe stats.
    uint64_t adj_len_evaluated = 0;
    uint64_t adj_len_matched   = 0;  // banded check passed (1-indel, 0-sub)
    uint64_t adj_len_absorbed  = 0;
    uint64_t adj_len_protected = 0;  // matched but LR ≤ threshold (or legacy mode)
    // H=3 shadow: near-boundary H=2 edges (|S| < 1.0 nat). Flag-only, never absorbed.
    // Tracks LR score distribution near zero to inform future H=3 absorption threshold.
    uint64_t h3_shadow_absorbed_near_boundary  = 0;
    uint64_t h3_shadow_protected_near_boundary = 0;
    // T8.9 indel-rescue counters (--errcor-rescue-indels).
    double   rescue_index_size_mb     = 0.0;
    uint64_t rescue_index_overflow    = 0;
    uint64_t rescue_indexed_parents   = 0;
    uint64_t rescue_children_examined = 0;
    uint64_t rescue_skip_bundle_hot   = 0;
    uint64_t rescue_index_queries     = 0;
    uint64_t rescue_topk_truncated    = 0;
    uint64_t rescue_pairs_banded      = 0;
    uint64_t rescue_banded_reject     = 0;
    std::array<uint64_t, 3> rescue_banded_ed{};   // ed=0/1/2
    uint64_t rescue_absorbed          = 0;
    uint64_t rescue_protected         = 0;        // S ≤ 0
    std::array<uint64_t, 6> rescue_absorbed_by_occ{};
    // Phase B3 counters.
    uint64_t b3_candidates             = 0;
    uint64_t b3_absorbed               = 0;
    uint64_t b3_protected              = 0;
    uint64_t b3_cross_bundle_protected = 0;
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
            for (int o = 0; o < 6; ++o) {
                log_info("    occ_bin[" + std::to_string(o) + "] absorbed/protected = " +
                         std::to_string(edge_lr_absorbed_by_occ[o]) + " / " +
                         std::to_string(edge_lr_protected_by_occ[o]));
            }
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
        if (h3_shadow_absorbed_near_boundary + h3_shadow_protected_near_boundary > 0) {
            log_info("H=3 shadow (near-boundary H=2, |S|<1.0 nat):");
            log_info("  Absorbed  near-boundary: " + std::to_string(h3_shadow_absorbed_near_boundary));
            log_info("  Protected near-boundary: " + std::to_string(h3_shadow_protected_near_boundary));
        }
        if (rescue_pairs_banded > 0 || rescue_indexed_parents > 0) {
            char b[64];
            std::snprintf(b, sizeof(b), "%.2f", rescue_index_size_mb);
            log_info("Phase 3 indel rescue (T8):");
            log_info("  Index size (MB)   : " + std::string(b));
            log_info("  Index parents     : " + std::to_string(rescue_indexed_parents));
            log_info("  Index hot drops   : " + std::to_string(rescue_index_overflow));
            log_info("  Children examined : " + std::to_string(rescue_children_examined));
            log_info("  Skipped (bundle hot): " + std::to_string(rescue_skip_bundle_hot));
            log_info("  Index queries     : " + std::to_string(rescue_index_queries));
            log_info("  Topk truncated    : " + std::to_string(rescue_topk_truncated));
            log_info("  Pairs banded      : " + std::to_string(rescue_pairs_banded));
            log_info("  Banded reject     : " + std::to_string(rescue_banded_reject));
            log_info("  Banded ed=0/1/2   : " + std::to_string(rescue_banded_ed[0]) +
                     " / " + std::to_string(rescue_banded_ed[1]) +
                     " / " + std::to_string(rescue_banded_ed[2]));
            log_info("  Absorbed          : " + std::to_string(rescue_absorbed));
            log_info("  Protected (S<=0)  : " + std::to_string(rescue_protected));
            for (int o = 0; o < 6; ++o) {
                log_info("    occ_bin[" + std::to_string(o) + "] absorbed = " +
                         std::to_string(rescue_absorbed_by_occ[o]));
            }
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
    uint32_t fwd_count = 0;  // reads arriving in canonical (fwd) orientation

    IndexEntry() : record_index(0), count(1), seq_id(0), damage_score(0), fwd_count(0) {}
    explicit IndexEntry(uint64_t idx) : record_index(idx), count(1), seq_id(0), damage_score(0), fwd_count(0) {}
};

// 2-bit packed sequence arena. ~4× memory reduction vs ASCII.
// Non-ACGT bases are encoded as A=0 and the sequence is flagged ineligible
// so Phase 3 skips it.
struct SeqArena {
    std::vector<uint8_t>  packed;
    std::vector<uint64_t> offsets;
    std::vector<uint16_t> lengths;
    std::vector<uint8_t>  eligible;  // 1 byte/read, no proxy bit-extract

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
        eligible.push_back(ok ? uint8_t{1} : uint8_t{0});
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
        eligible.push_back(ok ? uint8_t{1} : uint8_t{0});
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

    // Open-addressing directory over `keys`, indexed by key (keys are already
    // splitmix64-mixed by pair_key, so the low bits are uniform — no re-hash).
    // Each slot stores the key INLINE alongside idx_plus1 (ki+1 into
    // keys/offsets; 0 = empty). Built once after the CSR is filled, read-only
    // during parallel B1. Inlining the key means a probe reads a single cache
    // line — the previous slots[h]→keys[s-1] pair was TWO dependent cold loads
    // (the 34% child-query cost). The found idx and ids[] emission order are
    // identical, so dedup output is byte-for-byte unchanged.
    // Slot stores a 32-bit fingerprint (key>>32) instead of the full 64-bit key:
    // 8B/slot vs the 16B a {u64,u32} pays after alignment padding, halving `dir`.
    // The low bits of `key` drive the bucket; the high 32 bits are an independent
    // fingerprint. A fp match is only a candidate — query() verifies keys[idx]==key
    // (one load into the resident keys[]) before emitting, so fingerprint collisions
    // cost an extra probe, never an extra candidate. Output is byte-for-byte
    // unchanged. keys[] is already resident, so verification adds no new structure.
    struct Slot {
        uint32_t fp;          // key >> 32; valid only when idx_plus1 != 0
        uint32_t idx_plus1;   // 0 = empty; real idx 0 maps to idx_plus1 == 1
    };
    std::vector<Slot> dir;
    uint64_t slot_mask = 0;

    void build_directory() {
        size_t cap = 8;
        while (cap < keys.size() * 2) cap <<= 1;
        slot_mask = cap - 1;
        dir.assign(cap, Slot{0, 0});
        for (size_t ki = 0; ki < keys.size(); ++ki) {
            size_t h = keys[ki] & slot_mask;
            while (dir[h].idx_plus1) h = (h + 1) & slot_mask;
            dir[h] = Slot{static_cast<uint32_t>(keys[ki] >> 32),
                          static_cast<uint32_t>(ki + 1)};
        }
    }

    template <typename F>
    void query(uint64_t key, F&& fn) const {
        const uint32_t fp = static_cast<uint32_t>(key >> 32);
        for (size_t h = key & slot_mask;; h = (h + 1) & slot_mask) {
            const Slot& sl = dir[h];
            if (!sl.idx_plus1) return;
            if (sl.fp == fp) {
                uint32_t idx = sl.idx_plus1 - 1;
                if (keys[idx] == key) {
                    for (uint32_t i = offsets[idx]; i < offsets[idx + 1]; ++i)
                        fn(ids[i]);
                    return;
                }
            }
        }
    }
};

}  // namespace fqdup::derep_detail
