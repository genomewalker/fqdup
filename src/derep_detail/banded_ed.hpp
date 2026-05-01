#pragma once
// Banded edit distance ≤2 with Ukkonen-bounded DP and a small EditScript.
// Internal to the T8 indel-rescue path.

#include <algorithm>
#include <array>
#include <cstdint>

namespace fqdup::derep_detail {

struct EditScript {
    int8_t n_sub;
    int8_t n_ins;   // bases present in `b` (child) absent in `a` (parent)
    int8_t n_del;   // bases present in `a` absent in `b`
    // Up to 2 events; positions in `a` (parent) coordinates for sub/del,
    // in inserted-base position for ins (insert before parent position p).
    std::array<int16_t, 2> pos;
    std::array<uint8_t, 2> ref_base;  // 2-bit; 0xFF when N/A (e.g. ins)
    std::array<uint8_t, 2> alt_base;  // 2-bit; 0xFF when N/A (e.g. del)
    int8_t kind[2];                   // 0=sub, 1=ins, 2=del

    void clear() {
        n_sub = n_ins = n_del = 0;
        pos = {0, 0};
        ref_base = {0xFF, 0xFF};
        alt_base = {0xFF, 0xFF};
        kind[0] = kind[1] = -1;
    }
    int total() const { return n_sub + n_ins + n_del; }
};

// Banded edit distance with cap = 2.
// Returns -1 if min ed > 2 (early-exit), else returns the actual ed in {0,1,2}.
// Fills `out` with the edit script for the optimal alignment found.
// a/b are 2-bit base arrays (0..3); La and Lb are their lengths.
// Requires |La - Lb| <= 2 (caller filters via length-bin gate).
inline int banded_edit_distance_le2(const uint8_t* a, int La,
                                    const uint8_t* b, int Lb,
                                    EditScript* out) {
    constexpr int kCap = 2;
    if (out) out->clear();
    if (std::abs(La - Lb) > kCap) return -1;

    // Standard DP with band radius kCap. Memory: O(La).
    // dp[i][j] where j ∈ [i-kCap, i+kCap]; we store 2 rows.
    constexpr int kBand = 2 * kCap + 1;  // 5
    auto idx = [&](int j_off) { return j_off + kCap; };  // map j-i to [0,4]

    // prev / curr rows over band offsets.
    // INT_MAX_SAFE > 2*kCap + 1 marks "out of band / impossible".
    constexpr int kInf = 99;
    std::array<int, kBand> prev{}, curr{};
    // Backtrack hint per cell (parent row offset): -1 ins, 0 match/sub, +1 del.
    // We don't keep the full traceback table for ed≤2; instead, on ed==1 or 2
    // we walk a tiny diff routine afterward.

    // Init row i=0: dp[0][j] = j (insertions of b prefix into a).
    for (int j_off = -kCap; j_off <= kCap; ++j_off) {
        int j = j_off;  // i=0 → j = j_off
        if (j < 0 || j > Lb) prev[idx(j_off)] = kInf;
        else prev[idx(j_off)] = j;
    }

    for (int i = 1; i <= La; ++i) {
        int row_min = kInf;
        for (int j_off = -kCap; j_off <= kCap; ++j_off) {
            int j = i + j_off;
            if (j < 0 || j > Lb) { curr[idx(j_off)] = kInf; continue; }

            int best = kInf;
            // diagonal (sub or match): prev[j_off]
            if (j >= 1) {
                int sub_cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
                best = std::min(best, prev[idx(j_off)] + sub_cost);
            } else if (j == 0 && i >= 1) {
                // alignment of empty-prefix b to a[..i] = i deletions
                best = std::min(best, i);
            }
            // delete a[i-1] (consume parent base, no child base): prev[j_off+1]
            if (j_off + 1 <= kCap) {
                best = std::min(best, prev[idx(j_off + 1)] + 1);
            }
            // insert b[j-1] (consume child base, no parent base): curr[j_off-1]
            if (j_off - 1 >= -kCap && j >= 1) {
                best = std::min(best, curr[idx(j_off - 1)] + 1);
            }
            curr[idx(j_off)] = best;
            row_min = std::min(row_min, best);
        }
        if (row_min > kCap) return -1;  // band-wide pruning
        prev = curr;
    }

    // Final value: dp[La][Lb]; j_off = Lb - La
    int j_off_final = Lb - La;
    if (std::abs(j_off_final) > kCap) return -1;
    int ed = prev[idx(j_off_final)];
    if (ed > kCap) return -1;

    // Build script via a second tiny pass: a banded traceback re-runs DP and
    // records moves. For ed ≤ 2 with kBand=5, this is cheap and avoids
    // storing the full DP table during the forward pass.
    if (out && ed > 0) {
        // Re-run the DP, this time storing per-cell direction.
        // dir[i][j_off] ∈ {0:diag, 1:up(=del a[i-1]), 2:left(=ins b[j-1])}
        // Thread-local scratch — grown once, reused across all calls. Avoids
        // hot-loop allocator contention under T8's parallel workers.
        thread_local std::vector<std::array<int8_t, kBand>> dir;
        thread_local std::vector<std::array<int, kBand>> tab;
        if (dir.size() < static_cast<size_t>(La + 1)) dir.resize(La + 1);
        if (tab.size() < static_cast<size_t>(La + 1)) tab.resize(La + 1);
        for (int j_off = -kCap; j_off <= kCap; ++j_off) {
            int j = j_off;
            tab[0][idx(j_off)] = (j >= 0 && j <= Lb) ? j : kInf;
            dir[0][idx(j_off)] = (j > 0) ? 2 : 0;
        }
        for (int i = 1; i <= La; ++i) {
            for (int j_off = -kCap; j_off <= kCap; ++j_off) {
                int j = i + j_off;
                if (j < 0 || j > Lb) { tab[i][idx(j_off)] = kInf; dir[i][idx(j_off)] = 0; continue; }
                int best = kInf, bd = 0;
                if (j >= 1) {
                    int c = tab[i-1][idx(j_off)] + ((a[i-1] == b[j-1]) ? 0 : 1);
                    if (c < best) { best = c; bd = 0; }
                } else if (j == 0) {
                    if (i < best) { best = i; bd = 1; }
                }
                if (j_off + 1 <= kCap) {
                    int c = tab[i-1][idx(j_off + 1)] + 1;
                    if (c < best) { best = c; bd = 1; }
                }
                if (j_off - 1 >= -kCap && j >= 1) {
                    int c = tab[i][idx(j_off - 1)] + 1;
                    if (c < best) { best = c; bd = 2; }
                }
                tab[i][idx(j_off)] = best;
                dir[i][idx(j_off)] = bd;
            }
        }

        // Walk back from (La, Lb).
        int i = La, j = Lb;
        int slot = 0;
        while ((i > 0 || j > 0) && slot < 2) {
            int j_off = j - i;
            int8_t d = (j_off >= -kCap && j_off <= kCap) ? dir[i][idx(j_off)] : 0;
            if (d == 0) {
                if (i > 0 && j > 0 && a[i - 1] != b[j - 1]) {
                    out->kind[slot] = 0;
                    out->pos[slot] = static_cast<int16_t>(i - 1);
                    out->ref_base[slot] = a[i - 1];
                    out->alt_base[slot] = b[j - 1];
                    ++out->n_sub;
                    ++slot;
                }
                --i; --j;
            } else if (d == 1) {
                // deletion: a[i-1] consumed, b unchanged
                out->kind[slot] = 2;
                out->pos[slot] = static_cast<int16_t>(i - 1);
                out->ref_base[slot] = a[i - 1];
                out->alt_base[slot] = 0xFF;
                ++out->n_del; ++slot; --i;
            } else {
                // insertion: b[j-1] consumed, a unchanged
                out->kind[slot] = 1;
                out->pos[slot] = static_cast<int16_t>(i);  // insertion before parent pos i
                out->ref_base[slot] = 0xFF;
                out->alt_base[slot] = b[j - 1];
                ++out->n_ins; ++slot; --j;
            }
        }
    }

    return ed;
}

}  // namespace fqdup::derep_detail
