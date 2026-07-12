        if (errcor_.adj_len_probe && errcor_.legacy_veto)
            log_warn("--errcor-adj-len: requires the empirical posterior-odds rule (legacy SNP veto has no indel notion); probe disabled.");
        if (errcor_.adj_len_probe && !errcor_.legacy_veto && !fqcl_path_.empty())
            log_warn("--errcor-adj-len + --cluster-format: indel-edge absorptions are NOT recorded in the .fqcl genealogy "
                     "(ChildMismatch encodes substitutions only). Dedup output is correct, but the genealogy will under-report "
                     "edges. Disable --errcor-adj-len if exact edge accounting is required.");
        log_info(errcor_.legacy_veto
                 ? "Phase 3 decision rule: legacy SNP veto"
                 : "Phase 3 decision rule: empirical posterior odds (S > 0)");
        if (arena_.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
            throw std::runtime_error("Phase 3: arena size exceeds uint32_t range");
        const uint32_t N = static_cast<uint32_t>(arena_.size());
        if (N == 0) return;

        // The trailing pad byte extract_packed_part needs is now structural:
        // BlockBytes keeps kSlack readable bytes at the end of every block
        // (arena.hpp), so the last record's one-byte over-read is in bounds
        // without appending to the arena here.

        // T5.3: precompute bundle_key + bundle occupancy for LR prior. The
        // bundle_key (start_kmer ⊕ end_kmer) is our reference-free locus
        // proxy — clusters sharing a key landed at the same capture site, so
        // a high-occupancy bundle is evidence that mismatches in this group
        // could be PCR siblings rather than independent capture events.
        std::vector<uint64_t> bundle_key_of(N, 0);
        // For the LR prior we only need bundle occupancy (a count). The full
        // member list is only needed by the T5.6 adj-len probe, which uses a
        // separate length-agnostic map. Keep this one count-only — saves a
        // vector header per bundle (~24 B × ~N/avg_occ on capture data).
        // Occupancy is only ever read through bundle_occ_of() below, which returns 1
        // for any key it fails to find. So ONLY keys with occupancy >= 2 need to be
        // stored: a singleton's absence already encodes its count. On this data
        // 352,925,555 of 354,647,264 keys are singletons (99.5%), so a map over every
        // key spent 12.00 GB recording the value 1. (It was worse before: reserve(N)
        // took 24 GB, because ska's reserve() silently clamps the load factor to 0.5.)
        //
        // Count by sorting a scratch copy of the keys, keep only the repeated ones,
        // then drop the copy. Stored values are unchanged, and the map is never
        // iterated (find/operator[] only), so this is byte-identical by construction.
        ska::flat_hash_map<uint64_t, uint32_t> bundle_occ_map;  // occ >= 2 only
        bundle_occ_map.max_load_factor(0.75f);
        {
            std::vector<uint8_t>  dec;
            std::vector<char>     asc;
            std::vector<uint64_t> keys;
            keys.reserve(N);
            for (uint32_t id = 0; id < N; ++id) {
                if (!arena_.is_eligible(id)) continue;
                int L = arena_.length(id);
                if ((int)dec.size() < L) { dec.resize(L); asc.resize(L); }
                arena_.decode_range(id, 0, L, dec.data());
                for (int i = 0; i < L; ++i) asc[i] = static_cast<char>(dec[i]);
                uint64_t k = fqdup::bundlekey::from_decoded(
                    asc.data(), L, fqdup::bundlekey::kDefaultEndK);
                bundle_key_of[id] = k;
                keys.push_back(k);
            }
            std::sort(keys.begin(), keys.end());
            for (size_t i = 0; i < keys.size();) {
                size_t j = i;
                while (j < keys.size() && keys[j] == keys[i]) ++j;
                if (j - i > 1)
                    bundle_occ_map.emplace(keys[i], static_cast<uint32_t>(j - i));
                i = j;
            }
        }
        auto bundle_occ_of = [&](uint32_t id) -> uint64_t {
            auto it = bundle_occ_map.find(bundle_key_of[id]);
            return it == bundle_occ_map.end() ? 1u : it->second;
        };

        // T8 (rescue): length-agnostic parallel bundle key + occupancy. Built
        // only when --errcor-rescue-indels is set so the standard path pays
        // nothing. Indel pairs (L vs L±1, L vs L±2) collapse to the same key
        // here, so the within-bundle gate in the rescue driver actually fires
        // for them. Occupancy is computed against this key so bundle_hot sees
        // the true PCR cloud size that spans neighbouring lengths.
        std::vector<uint64_t> rescue_bundle_key_of;
        ska::flat_hash_map<uint64_t, uint32_t> rescue_bundle_occ_map;
        if (errcor_.rescue_indels && !errcor_.legacy_veto && errcor_.empirical) {
            rescue_bundle_key_of.assign(N, 0);
            rescue_bundle_occ_map.max_load_factor(0.75f);
            std::vector<uint8_t>  dec;
            std::vector<char>     asc;
            std::vector<uint64_t> keys;   // occ >= 2 only; see bundle_occ_map above
            keys.reserve(N);
            for (uint32_t id = 0; id < N; ++id) {
                if (!arena_.is_eligible(id)) continue;
                int L = arena_.length(id);
                if ((int)dec.size() < L) { dec.resize(L); asc.resize(L); }
                arena_.decode_range(id, 0, L, dec.data());
                for (int i = 0; i < L; ++i) asc[i] = static_cast<char>(dec[i]);
                uint64_t k = fqdup::bundlekey::from_decoded_no_len(
                    asc.data(), L, fqdup::bundlekey::kDefaultEndK);
                rescue_bundle_key_of[id] = k;
                keys.push_back(k);
            }
            std::sort(keys.begin(), keys.end());
            for (size_t i = 0; i < keys.size();) {
                size_t j = i;
                while (j < keys.size() && keys[j] == keys[i]) ++j;
                if (j - i > 1)
                    rescue_bundle_occ_map.emplace(keys[i], static_cast<uint32_t>(j - i));
                i = j;
            }
        }
        auto rescue_bundle_occ_of = [&](uint32_t id) -> uint64_t {
            if (rescue_bundle_key_of.empty()) return 1u;
            auto it = rescue_bundle_occ_map.find(rescue_bundle_key_of[id]);
            return it == rescue_bundle_occ_map.end() ? 1u : it->second;
        };

        is_error_.assign(N, false);

        std::vector<uint32_t> id_count(N, 0);
        std::vector<uint32_t> id_fwd_count(N, 0);
        // Reverse map: seq_id → IndexEntry* for representative propagation during absorption.
        std::vector<IndexEntry*> seq_entry(N, nullptr);
        index_.for_each([&](const SequenceFingerprint&, IndexEntry& entry) {
            id_count[entry.seq_id] = entry.count;
            id_fwd_count[entry.seq_id] = entry.fwd_count;
            seq_entry[entry.seq_id] = &entry;
        });

        // Accumulated mass per cluster: starts as Pass-1 count, grows as children
        // are absorbed.  Used for B3 boundary check only — the SNP veto still uses
        // id_count (original Pass-1 count) so the ratio threshold remains calibrated.
        std::vector<uint32_t> acc_count(id_count);

        // Scratch for packed interior hashing
        std::vector<uint8_t> scratch;
        scratch.reserve(65535);

        // Number of interior positions from each end considered "near-damage zone".
        // Positions just beyond the mask boundary retain residual deamination damage
        // (exponential tail, typically 5-10% at the first unmasked position).
        // C↔T / G↔A mismatches (xr=2) at these positions are treated as damage
        // variants and bypass the SNP veto, rather than being protected as SNPs.
        constexpr int kDamageEdgeMargin = 5;
        // Minimum interior length for a meaningful 3-way pigeonhole split.
        // Below this, splits degenerate (parts of 0-4 bp) and the hash keys
        // become unreliable, risking false absorptions.
        constexpr int kMinInteriorLen = 20;
        uint64_t short_interior_skipped = 0;
        // Short-read policy (unconditional): reads whose masked interior is shorter
        // than kMinInteriorLen skip Phase-3 1-mismatch error-correction and pass
        // through as their own representative. Phase 1 (XXH128 full-seq fingerprint)
        // already exact-dedups them, so no duplicate survives. Below ilen 20 the
        // 4-way pigeonhole split degenerates (parts of 0-4 bp -> unreliable keys) and
        // the only alternative -- an O(bucket^2) same-length brute-force scan -- both
        // avalanches on the huge short-read buckets AND over-merges distinct short
        // molecules (the H=1 neighbour graph saturates). The threshold 20 = max of
        // pigeonhole viability (4 parts x ~5 bp) and chance-neighbour discriminability
        // (d(l,B) = (B-1).l.(3/4).(1/4)^(l-1) < 1 needs l >= 13-15 at KapK bucket
        // sizes). Keeping short reads as-is is O(n) and lossless in the safe direction.

        struct InteriorLayout {
            int  k5 = 0, ilen = 0;
            int  s0 = 0, s1 = 0, s2 = 0, s3 = 0;
            int  nb0 = 0, nb1 = 0, nb2 = 0, nb3 = 0;
            int  nbytes = 0;
            bool ready = false;
        };
        std::vector<InteriorLayout> layout_cache(65536);
        auto get_layout = [&](int L) -> const InteriorLayout& {
            auto& lay = layout_cache[static_cast<uint16_t>(L)];
            if (!lay.ready) {
                auto [k5, k3] = damage_zone_bounds(L, profile_);
                lay.k5   = k5;
                lay.ilen = std::max(0, L - k5 - k3);
                if (lay.ilen >= kMinInteriorLen) {
                    split4_lens(lay.ilen, lay.s0, lay.s1, lay.s2, lay.s3);
                    lay.nb0    = (lay.s0 + 3) / 4;
                    lay.nb1    = (lay.s1 + 3) / 4;
                    lay.nb2    = (lay.s2 + 3) / 4;
                    lay.nb3    = (lay.s3 + 3) / 4;
                    lay.nbytes = lay.nb0 + lay.nb1 + lay.nb2 + lay.nb3;
                }
                lay.ready = true;
            }
            return lay;
        };

        // Single EC-eligibility predicate: a read may take part in ANY Phase-3
        // error-correction path only if its damage-masked interior is >=
        // kMinInteriorLen. Below that the pigeonhole split degenerates and the
        // key is unreliable, so the read stays its own representative. The main
        // H=1/H=2 loop guards on lay.ilen inline; this same gate is applied to
        // the adj-len and rescue-indels bundle builds (which key on is_eligible
        // alone) so no correction path can silently absorb a short-interior read.
        auto ec_eligible = [&](uint32_t id) -> bool {
            return get_layout(static_cast<int>(arena_.length(id))).ilen >= kMinInteriorLen;
        };

        Phase3Stats stats;
        using clk = std::chrono::steady_clock;

        // ── Phase A: index parents (count > min_parent_count) ──────────────
        struct BuildEntry { uint64_t key; uint32_t id; };
        // Tags 0-5: pair-hash for all C(4,2)=6 combinations of 4 parts.
        // 4-way pigeonhole: if H(child,parent)≤2, at least 2 of 4 parts match →
        // one of the 6 pair-keys fires.  Covers both H=1 and H=2.
        static constexpr int kPairA[6] = {0, 0, 0, 1, 1, 2};
        static constexpr int kPairB[6] = {1, 2, 3, 2, 3, 3};
        struct LenShard { std::array<FlatPairIndex, 6> pi; };

        // Group eligible ids by interior length so each ilen owns a contiguous
        // window in `order`. Sorted by (ilen, id) ascending, so an ilen's window
        // lists its ids in exactly the order the old [0,N) scan visited them --
        // which is what lets build_shard() fill that ilen's postings just-in-time
        // and still produce bit-identical FlatPairIndex buckets.
        std::vector<uint32_t> order;
        struct IlenSeg { int ilen; uint32_t lo, hi; };
        std::vector<IlenSeg> ilen_segs;
        {
            std::vector<std::pair<int, uint32_t>> tmp;
            tmp.reserve(N);
            for (uint32_t id = 0; id < N; ++id) {
                if (!arena_.is_eligible(id)) continue;
                int il = get_layout(arena_.length(id)).ilen;
                tmp.emplace_back(il, id);  // include tiny ids too: the worker visits
            }                              // and counts them exactly as the [0,N) scan did
            std::sort(tmp.begin(), tmp.end(),
                      [](const std::pair<int, uint32_t>& a,
                         const std::pair<int, uint32_t>& b) {
                          return a.first != b.first ? a.first < b.first
                                                    : a.second < b.second;
                      });
            order.reserve(tmp.size());
            for (size_t i = 0; i < tmp.size();) {
                int il = tmp[i].first;
                uint32_t lo = static_cast<uint32_t>(order.size());
                while (i < tmp.size() && tmp[i].first == il) {
                    order.push_back(tmp[i].second);
                    ++i;
                }
                ilen_segs.push_back({il, lo, static_cast<uint32_t>(order.size())});
            }
        }

        uint64_t n_parents = 0;
        for (const auto& seg : ilen_segs)
            if (seg.ilen >= kMinInteriorLen) n_parents += seg.hi - seg.lo;
        stats.parents_indexed = n_parents;
        log_info("Phase 3: indexed " + std::to_string(n_parents) +
                 " sequences (directed monotone ascent, all eligible)");

        // Streaming Phase-3: instead of building every length-shard up front and
        // holding them all resident (the ~33GB FlatPairIndex peak), the B1 driver
        // below builds ONE ilen's shard just-in-time, runs that ilen's children
        // against it, then frees it before the next ilen. Candidate edges are
        // strictly intra-ilen (B1 requires parent.len==child.len), and the output
        // is invariant to the order in which edges are collected (proven: identical
        // md5 across thread counts / work-stealing orders), so per-ilen streaming
        // is byte-for-byte identical to the all-resident build. `build_map` is kept
        // live and consumed one ilen at a time. Peak collapses to floor +
        // largest-single-shard.
        ska::flat_hash_map<int, LenShard> shards;
        std::array<uint64_t, 8> bhist{};

        // Largest single shard's FlatPairIndex footprint, sampled while that shard is
        // still live. The mem breakdown below runs at B2, after every shards.erase(),
        // so it would otherwise report FlatPairIndex=0 and hide ~20% of the peak.
        uint64_t fpi_peak_b = 0;

        // Distinct-4mer complexity of an interior, read directly off the 2-bit packed
        // arena (this streaming engine packs sequences 4 bases/byte; the old ASCII
        // arena is gone). Repetitive runs (homopolymer/STR/telomere/minisatellite)
        // collapse to few distinct 4-mers regardless of length; real DNA keeps most
        // positions distinct. Used to decide whether an oversized candidate bucket is
        // a technical explosion (cap it) or genuinely abundant biology (keep every
        // posting). `packed` is the read's packed arena data, `start` the interior's
        // first base index (= layout.k5); eligible ids are pure ACGT, so every window
        // is valid and the distinct-4mer mapping is a bijection of the base encoding.
        auto interior_low_complexity = [](const uint8_t* packed, int start, int ilen,
                                          double frac) -> bool {
            if (ilen < 8) return false;              // too short to judge → never cap
            uint64_t seen[4] = {0,0,0,0};
            uint32_t kmer = 0; int valid = 0;
            for (int p = 0; p < ilen; ++p) {
                int pos = start + p;
                uint32_t c = (packed[pos >> 2] >> (6 - 2 * (pos & 3))) & 0x3u;
                kmer = ((kmer << 2) | c) & 0xFFu;    // trailing 4 bases (8 bits)
                if (p >= 3) { seen[kmer >> 6] |= (1ull << (kmer & 63)); ++valid; }
            }
            if (valid <= 0) return false;
            int distinct = __builtin_popcountll(seen[0]) + __builtin_popcountll(seen[1])
                         + __builtin_popcountll(seen[2]) + __builtin_popcountll(seen[3]);
            return (double)distinct < frac * (double)valid;
        };
        // Entropy-cap visibility, accumulated across all per-ilen build_shard calls:
        // oversized buckets seen, capped (low-complexity), spared (high-complexity →
        // real abundant biology, left uncapped).
        int64_t ent_oversized = 0, ent_capped = 0, ent_spared = 0;

        // Largest single ilen's posting arrays, sampled while live (same reason as
        // fpi_peak_b: they are gone by the time the mem breakdown runs).
        uint64_t post_peak_b = 0;

        auto build_shard = [&](const IlenSeg& seg) -> LenShard& {
            const int ilen = seg.ilen;
            auto& sh = shards[ilen];

            // Postings for THIS ilen only, built just-in-time from the ilen's window
            // in `order`. Materialising every ilen's postings up front cost
            // 6 x sizeof(BuildEntry)(16 B) x n_parents(350,534,710) = 33.65 GB before
            // std::vector growth slack -- the largest single term in the run, and
            // invisible to the mem breakdown because it was fully drained by the time
            // that logger fires. Peak collapses to the largest single ilen.
            std::array<std::vector<BuildEntry>, 6> tag_entries;
            {
                std::vector<uint8_t> scratch;
                for (uint32_t oi = seg.lo; oi < seg.hi; ++oi) {
                    const uint32_t id  = order[oi];
                    const auto&    lay = get_layout(arena_.length(id));

                    auto t0 = clk::now();
                    if (scratch.size() < static_cast<size_t>(lay.nbytes))
                        scratch.resize(lay.nbytes);
                    const uint8_t* psrc = arena_.data(id);
                    int starts[4] = {lay.k5,
                                     lay.k5 + lay.s0,
                                     lay.k5 + lay.s0 + lay.s1,
                                     lay.k5 + lay.s0 + lay.s1 + lay.s2};
                    int sizes[4]  = {lay.s0, lay.s1, lay.s2, lay.s3};
                    int nbs[4]    = {lay.nb0, lay.nb1, lay.nb2, lay.nb3};
                    uint8_t* parts[4] = {scratch.data(),
                                         scratch.data() + lay.nb0,
                                         scratch.data() + lay.nb0 + lay.nb1,
                                         scratch.data() + lay.nb0 + lay.nb1 + lay.nb2};
                    for (int p = 0; p < 4; ++p)
                        extract_packed_part(psrc, starts[p], sizes[p], parts[p]);
                    uint64_t h[4];
                    for (int p = 0; p < 4; ++p)
                        h[p] = XXH3_64bits(parts[p], nbs[p]);
                    auto t1 = clk::now();

                    for (int t = 0; t < 6; ++t)
                        tag_entries[t].push_back(
                            {pair_key(h[kPairA[t]], h[kPairB[t]], t, lay.ilen), id});
                    auto t2 = clk::now();

                    stats.decode_hash_parent_ms +=
                        std::chrono::duration<double, std::milli>(t1 - t0).count();
                    stats.insert_ms +=
                        std::chrono::duration<double, std::milli>(t2 - t1).count();
                }
                uint64_t pb = 0;
                for (const auto& ev : tag_entries)
                    pb += static_cast<uint64_t>(ev.capacity()) * sizeof(BuildEntry);
                post_peak_b = std::max(post_peak_b, pb);
            }

            for (int t = 0; t < 6; ++t) {
                auto& ev = tag_entries[t];
                std::sort(ev.begin(), ev.end(),
                          [](const BuildEntry& a, const BuildEntry& b){ return a.key < b.key; });
                auto& pi = sh.pi[t];
                size_t i = 0;
                pi.offsets.push_back(0);
                while (i < ev.size()) {
                    uint64_t k = ev[i].key;
                    pi.keys.push_back(k);
                    // Bucket extent [i, j): all postings sharing this part-key.
                    size_t j = i;
                    while (j < ev.size() && ev[j].key == k) ++j;
                    // Decide the effective cap for THIS bucket. In entropy-aware mode an
                    // oversized bucket is capped only when its representative interior is
                    // low-complexity; complex buckets (real abundant reads) stay uncapped.
                    uint32_t cap = errcor_.bucket_cap;  // 0 = no cap
                    if (cap && errcor_.bucket_cap_lowcomplexity_only && (j - i) > cap) {
                        ++ent_oversized;
                        const auto& rlay = get_layout(arena_.length(ev[i].id));
                        if (!interior_low_complexity(arena_.data(ev[i].id), rlay.k5,
                                                     rlay.ilen, errcor_.bucket_cap_complexity_frac)) {
                            cap = 0;  // high-complexity → keep every posting
                            ++ent_spared;
                        } else {
                            ++ent_capped;
                        }
                    }
                    uint32_t cnt = 0;
                    for (; i < j; ++i) {
                        if (cap == 0 || cnt < cap) {
                            pi.ids.push_back(ev[i].id); cnt++;
                        } else {
                            stats.bucket_overflow_drops++;
                        }
                    }
                    pi.offsets.push_back(static_cast<uint32_t>(pi.ids.size()));
                    uint32_t blen = pi.offsets.back() - pi.offsets[pi.offsets.size() - 2];
                    unsigned b = blen == 0 ? 0 : 31 - __builtin_clz(blen);
                    bhist[std::min(b, 7u)]++;
                }
                pi.build_directory();  // O(1) key lookup, replaces per-probe lower_bound
                ev.clear();
                ev.shrink_to_fit();    // drop this tag's postings as soon as it is indexed
            }
            return sh;
        };

        // ── Phase B1: collect ChildMismatch records (parallel) ──────────────
        // Pre-warm layout_cache before workers start so it becomes read-only, and
        // build a per-id interior slab: extract each eligible long-interior id's
        // full interior ONCE here, so the per-candidate B1 loop reads it (a pointer
        // deref) instead of re-extracting it per child. That re-extraction was the
        // 38% __memmove_evex in the perf profile. Output-preserving: at the B1 use
        // site the parent's length == the child's length (checked there), so the
        // id's own layout used here is the exact layout used at the comparison.
        std::vector<uint64_t> interior_off(static_cast<size_t>(N) + 1, 0);
        for (uint32_t id = 0; id < N; ++id) {
            if (!arena_.is_eligible(id)) continue;
            const auto& lay = get_layout(arena_.length(id));
            if (lay.ilen >= kMinInteriorLen)
                interior_off[id + 1] = static_cast<uint64_t>((lay.ilen + 3) / 4);
        }
        for (uint32_t id = 0; id < N; ++id)
            interior_off[id + 1] += interior_off[id];
        std::vector<uint8_t> interior_slab(interior_off[N]);
        for (uint32_t id = 0; id < N; ++id) {
            if (!arena_.is_eligible(id)) continue;
            const auto& lay = get_layout(arena_.length(id));
            if (lay.ilen < kMinInteriorLen) continue;
            extract_packed_part(arena_.data(id), lay.k5, lay.ilen,
                                interior_slab.data() + interior_off[id]);
        }
        auto interior_ptr = [&](uint32_t id) -> const uint8_t* {
            return interior_slab.data() + interior_off[id];
        };

        // `order` / `ilen_segs` are built above, before the parent index, because
        // build_shard() now drives its just-in-time posting fill off each ilen's
        // window.

        unsigned n_threads = errcor_.threads;
        if (n_threads == 0) n_threads = std::max(1u, std::thread::hardware_concurrency());

        struct B1LocalStats {
            uint64_t children_scanned       = 0;
            uint64_t total_candidates       = 0;
            uint64_t children_found         = 0;
            uint64_t short_interior_skipped = 0;
            uint64_t short_brute_evaluated  = 0;
            uint64_t short_brute_found      = 0;
            uint64_t short_too_small_skipped = 0;
            double   decode_hash_child_ms   = 0;
            double   query_ms               = 0;
            double   check_ms               = 0;
        };

        std::vector<std::vector<ChildMismatch>> per_thread_mm(n_threads);
        std::vector<B1LocalStats>              per_thread_stats(n_threads);

        constexpr uint32_t kChunkSize = 256;
        // Per-ilen dispatch window into `order`; the driver sets [seg_lo,seg_hi)
        // and seg_chunks and resets next_chunk before launching each ilen's workers.
        uint32_t seg_lo = 0, seg_hi = 0, seg_chunks = 0;
        std::atomic<uint32_t> next_chunk{0};

        auto b1_worker = [&](unsigned tid) {
            auto& local_mm = per_thread_mm[tid];
            auto& ls       = per_thread_stats[tid];
            local_mm.reserve(std::min(static_cast<size_t>(N) * 2 / n_threads + 1024,
                                      static_cast<size_t>(1 << 22)));

            std::vector<uint8_t>  scratch;
            scratch.reserve(65535);
            std::vector<uint32_t> cand_storage;
            if (errcor_.bucket_cap > 0)
                cand_storage.resize(12u * errcor_.bucket_cap);
            else
                cand_storage.reserve(256);

            while (true) {
                uint32_t chunk = next_chunk.fetch_add(1, std::memory_order_relaxed);
                if (chunk >= seg_chunks) break;
                uint32_t pos_lo = seg_lo + chunk * kChunkSize;
                uint32_t pos_hi = std::min(pos_lo + kChunkSize, seg_hi);

                for (uint32_t pos = pos_lo; pos < pos_hi; ++pos) {
                    uint32_t cid = order[pos];
                    if (is_error_[cid]) continue;
                    if (!arena_.is_eligible(cid)) continue;

                    // Timing here is diagnostic only (decode_hash_child_ms/query_ms/check_ms
                    // feed a log summary, never an algorithmic decision -- confirmed by reading
                    // every consumer of B1LocalStats). clk::now() was being called 4x per read
                    // on the full-scale hot path (profiled: ~60%+ of this worker's CPU time was
                    // in __vdso_clock_gettime, not the actual comparison logic). Sample 1-in-64
                    // and extrapolate (x kTimingStride) to keep the aggregate estimate roughly
                    // unbiased while cutting timestamp-call overhead ~64x.
                    constexpr uint32_t kTimingStride = 64;
                    const bool do_time = (cid % kTimingStride) == 0;

                    int L = arena_.length(cid);
                    const auto& lay = get_layout(L);
                    if (lay.ilen < kMinInteriorLen) {
                        // Unconditional short-read skip: Phase 1 already exact-deduped
                        // this read; below ilen 20 the pigeonhole index is unreliable,
                        // so it passes through as its own representative. See the
                        // policy note at the Phase-A index build.
                        ls.short_interior_skipped++;
                        continue;
                    }

                    auto shard_it = shards.find(lay.ilen);
                    if (shard_it == shards.end()) continue;

                    auto t0 = do_time ? clk::now() : clk::time_point{};
            // Scratch layout (per child):
            //   [0 .. nbytes)               : ci_parts  (4-part canonical, nb0+nb1+nb2+nb3 bytes)
            //   [nbytes .. nbytes+nf)       : ci_full   (full canonical interior)
            //   [nbytes+nf .. nbytes+2nf)   : crc_full  (RC of ci_full)
            //   [nbytes+2nf .. 2*nbytes+2nf): crc_parts (4-part RC, nb0+nb1+nb2+nb3 bytes)
            // (parent interior is read from the per-id slab, not scratch.)
            int nf = (lay.ilen + 3) / 4;
            size_t total_scratch = static_cast<size_t>(2 * lay.nbytes + 2 * nf);
            if (scratch.size() < total_scratch)
                scratch.resize(total_scratch);
            uint8_t* ci_parts  = scratch.data();
            uint8_t* ci_full   = ci_parts  + lay.nbytes;
            uint8_t* crc_full  = ci_full   + nf;
            uint8_t* crc_parts = crc_full  + nf;

            const uint8_t* psrc_c = arena_.data(cid);
            int starts[4] = {lay.k5,
                             lay.k5 + lay.s0,
                             lay.k5 + lay.s0 + lay.s1,
                             lay.k5 + lay.s0 + lay.s1 + lay.s2};
            int sizes[4]  = {lay.s0, lay.s1, lay.s2, lay.s3};
            int nbs[4]    = {lay.nb0, lay.nb1, lay.nb2, lay.nb3};
            uint8_t* ci_p[4] = {ci_parts,
                                ci_parts + lay.nb0,
                                ci_parts + lay.nb0 + lay.nb1,
                                ci_parts + lay.nb0 + lay.nb1 + lay.nb2};
            for (int p = 0; p < 4; ++p)
                extract_packed_part(psrc_c, starts[p], sizes[p], ci_p[p]);
            uint64_t h[4], rh[4];
            for (int p = 0; p < 4; ++p) h[p] = XXH3_64bits(ci_p[p], nbs[p]);

            extract_packed_part(psrc_c, lay.k5, lay.ilen, ci_full);

            // RC of full interior for orientation-aware comparison and RC hash keys.
            // A 1-base change can flip which orientation produces the smaller canonical hash,
            // so parent and child may be stored in opposite orientations in the arena.
            compute_interior_rc(ci_full, lay.ilen, crc_full);

            int rc_starts[4] = {0, lay.s0, lay.s0+lay.s1, lay.s0+lay.s1+lay.s2};
            uint8_t* crc_p[4] = {crc_parts,
                                 crc_parts + lay.nb0,
                                 crc_parts + lay.nb0 + lay.nb1,
                                 crc_parts + lay.nb0 + lay.nb1 + lay.nb2};
            for (int p = 0; p < 4; ++p)
                extract_packed_part(crc_full, rc_starts[p], sizes[p], crc_p[p]);
            for (int p = 0; p < 4; ++p) rh[p] = XXH3_64bits(crc_p[p], nbs[p]);
            if (do_time) {
                auto t1 = clk::now();
                ls.decode_hash_child_ms +=
                    std::chrono::duration<double, std::milli>(t1 - t0).count() * kTimingStride;
            }

            ls.children_scanned++;

            // Collect unique parent candidates — query both canonical and RC hash keys.
            // With bucket_cap=0 (default), the storage grows; with bucket_cap>0 it is
            // pre-sized and any overflow is counted (paired with bucket_overflow_drops).
            //
            // Tried and reverted: dedup-on-insert (linear scan before append), reasoning
            // that avg candidates/child measured at 1.23 on a real 5M-read run meant
            // std::sort+std::unique's overhead should dominate actual work. Measured
            // result was the opposite: 98.3s -> 112.7s wall clock, cache-miss rate
            // 10.17% -> 11.41%, LLC misses +29% (perf stat, same dataset, byte-identical
            // correctness both ways). Likely cause: introsort's small-N fast path
            // (insertion sort below its threshold) is already a fast, cache-friendly
            // single sequential pass; the interleaved linear-scan-per-insert during the
            // 12 hash lookups appears to hurt cache behavior more than the sort/unique
            // pass it replaced. Kept as sort+unique -- a small average candidate count
            // does not by itself justify removing a well-optimized generic algorithm.
            uint32_t n_cands = 0;
            // Entropy-aware mode leaves high-complexity buckets uncapped, so a query can
            // legitimately return more than 12*bucket_cap postings — grow on demand there.
            const bool unbounded = (errcor_.bucket_cap == 0) || errcor_.bucket_cap_lowcomplexity_only;
            auto collect = [&](uint32_t pid) {
                if (unbounded) {
                    if (n_cands >= cand_storage.size()) cand_storage.resize(n_cands + 1);
                    cand_storage[n_cands++] = pid;
                } else if (n_cands < static_cast<uint32_t>(cand_storage.size())) {
                    cand_storage[n_cands++] = pid;
                }
            };
            uint32_t* cand_buf = cand_storage.data();

            auto tq0 = do_time ? clk::now() : clk::time_point{};
            auto& sh = shard_it->second;
            for (int t = 0; t < 6; ++t) {
                sh.pi[t].query(pair_key(h[kPairA[t]],  h[kPairB[t]],  t, lay.ilen), collect);
                sh.pi[t].query(pair_key(rh[kPairA[t]], rh[kPairB[t]], t, lay.ilen), collect);
            }
            if (do_time) {
                auto tq1 = clk::now();
                ls.query_ms +=
                    std::chrono::duration<double, std::milli>(tq1 - tq0).count() * kTimingStride;
            }

            // Refresh: collect() may have reallocated cand_storage in unbounded mode.
            cand_buf = cand_storage.data();
            std::sort(cand_buf, cand_buf + n_cands);
            n_cands = static_cast<uint32_t>(
                std::unique(cand_buf, cand_buf + n_cands) - cand_buf);
            ls.total_candidates += n_cands;

            auto tc0 = do_time ? clk::now() : clk::time_point{};
            // Two-pass software-pipelined gather to deepen memory-level
            // parallelism. Each pi_buf gather below scatters into the 3.6GB
            // slab and is the dominant LLC miss (57% of B1). Pass A applies the
            // cheap guards and issues a prefetch for EVERY surviving parent's
            // interior, so all ~21 cold misses go outstanding together (vs one
            // K-ahead before, which also wasted slots on guard-failing pids).
            // Pass B runs the compares once the lines are in flight. Output-
            // neutral: survivors are compacted into cand_buf in candidate order
            // (n_surv <= ci, no clobber of unread entries) and the compare body
            // is unchanged, so the recorded edges are byte-identical.
            uint32_t n_surv = 0;
            for (uint32_t ci = 0; ci < n_cands; ++ci) {
                uint32_t pid = cand_buf[ci];
                if (is_error_[pid]) continue;
                if (!arena_.is_eligible(pid)) continue;
                if (arena_.length(pid) != static_cast<uint16_t>(L)) continue;
                if (pid == cid) continue;
                // Directed edge condition: absorb into higher-count sequences, OR into
                // equal-count sequences using seq_id as tiebreak (lower id = parent).
                // The tiebreak makes the equal-count DAG acyclic (edges always flow toward
                // lower seq_id) while still allowing singleton→singleton chains, which
                // enables H=2 absorption via intermediates: A(1)→B(1)→T(30).
                // Equal-count true SNP variants are still protected by the SNP veto
                // (both halves of the pair share the same mismatch → sig/parent = 100%).
                if (id_count[pid] < id_count[cid]) continue;
                if (id_count[pid] == id_count[cid] && pid >= cid) continue;
                __builtin_prefetch(interior_ptr(pid), 0, 0);
                cand_buf[n_surv++] = pid;
            }
            for (uint32_t si = 0; si < n_surv; ++si) {
                uint32_t pid = cand_buf[si];
                // Parent's full interior: read from the per-id slab (extracted
                // once). pid has length == L here (checked in Pass A), so its
                // slab layout equals lay -> byte-identical to a live extract.
                const uint8_t* pi_buf = interior_ptr(pid);

                // Try direct comparison (H=1 — same canonical orientation)
                MismatchInfo mm = packed_find_mismatch(pi_buf, ci_full, 0, lay.ilen, errcor_.protect_transversions);
                if (mm.found) {
                    local_mm.push_back({pid, cid, mm.position, mm.base_b, mm.base_a,
                                          0, 0, 0, 1, {}});
                    ls.children_found++;
                    continue;
                }
                // Try RC comparison (H=1)
                mm = packed_find_mismatch(pi_buf, crc_full, 0, lay.ilen, errcor_.protect_transversions);
                if (mm.found) {
                    uint16_t cpos = static_cast<uint16_t>(lay.ilen - 1 - mm.position);
                    local_mm.push_back({pid, cid, cpos,
                                          static_cast<uint8_t>(mm.base_b ^ 0x3u),
                                          static_cast<uint8_t>(mm.base_a ^ 0x3u),
                                          0, 0, 0, 1, {}});
                    ls.children_found++;
                    continue;
                }
                // H=2 path: only for reads with count ≤ max_h2_count.
                // With b1_damage_adjust=true (default): terminal damage mismatches
                // (C→T at 5', G→A at 3' for DS; both channels at both ends for SS)
                // are not counted toward H_adj. h_adj=2 → record both; h_adj=1 →
                // record only the non-damage mismatch as H=1; h_adj=0 → skip (both
                // terminal damage, handled by Phase 2 tag-folding).
                // With b1_damage_adjust=false: reverts to old transversion-only rule.
                if (id_count[cid] <= errcor_.max_h2_count) {
                    auto admit_h2 = [&](MismatchInfo2 mm2, bool rc) -> bool {
                        if (mm2.count != 2) return false;
                        // Transversion guard: A↔T / C↔G (xr==3) protected regardless
                        // of b1_damage_adjust when --protect-transversions is active.
                        if (errcor_.protect_transversions) {
                            if ((mm2.base_a[0] ^ mm2.base_b[0]) == 3u) return false;
                            if ((mm2.base_a[1] ^ mm2.base_b[1]) == 3u) return false;
                        }
                        bool dmg0, dmg1;
                        if (errcor_.b1_damage_adjust) {
                            dmg0 = is_terminal_damage_mismatch(
                                       mm2.base_a[0], mm2.base_b[0],
                                       mm2.pos[0], lay.ilen, profile_.ss_mode);
                            dmg1 = is_terminal_damage_mismatch(
                                       mm2.base_a[1], mm2.base_b[1],
                                       mm2.pos[1], lay.ilen, profile_.ss_mode);
                        } else {
                            dmg0 = is_damage_sub_packed(mm2.base_a[0], mm2.base_b[0], errcor_.protect_transversions);
                            dmg1 = is_damage_sub_packed(mm2.base_a[1], mm2.base_b[1], errcor_.protect_transversions);
                            if (dmg0 || dmg1) return false; // old: reject any damage type
                        }
                        int h_adj = 2 - (int)dmg0 - (int)dmg1;
                        if (h_adj == 0) return false; // both terminal damage → skip (Phase 2)
                        auto cvt = [&](int i, uint16_t p, uint8_t ab, uint8_t pb) {
                            if (rc) {
                                p  = static_cast<uint16_t>(lay.ilen - 1 - p);
                                ab = static_cast<uint8_t>(ab ^ 0x3u);
                                pb = static_cast<uint8_t>(pb ^ 0x3u);
                            }
                            return std::make_tuple(p, ab, pb);
                        };
                        if (h_adj == 1) {
                            int ri = dmg0 ? 1 : 0;
                            auto [p, ab, pb] = cvt(ri, mm2.pos[ri], mm2.base_b[ri], mm2.base_a[ri]);
                            local_mm.push_back({pid, cid, p, ab, pb, 0, 0, 0, 1, {}});
                        } else {
                            auto [p0, ab0, pb0] = cvt(0, mm2.pos[0], mm2.base_b[0], mm2.base_a[0]);
                            auto [p1, ab1, pb1] = cvt(1, mm2.pos[1], mm2.base_b[1], mm2.base_a[1]);
                            local_mm.push_back({pid, cid, p0, ab0, pb0, p1, ab1, pb1, 2, {}});
                        }
                        ls.children_found++;
                        return true;
                    };
                    MismatchInfo2 mm2 = packed_find_mismatches2(pi_buf, ci_full, 0, lay.ilen);
                    if (admit_h2(mm2, false)) continue;
                    mm2 = packed_find_mismatches2(pi_buf, crc_full, 0, lay.ilen);
                    if (admit_h2(mm2, true)) { /* admitted */ }
                }
            }
            if (do_time) {
                ls.check_ms +=
                    std::chrono::duration<double, std::milli>(clk::now() - tc0).count() * kTimingStride;
            }
                }  // end for cid in chunk
            }  // end while chunk dispatcher
        };  // end b1_worker

        // Streaming driver: process one ilen at a time. Build that ilen's shard
        // (long-interior ilens only; short ilens have no shard and their children
        // skip EC and pass through), dispatch B1 over its contiguous child window in `order`, then
        // free the shard before the next ilen. FlatPairIndex residency collapses
        // from sum-of-all-shards to the single largest shard. Edges are intra-ilen
        // and output is order-invariant, so this is byte-identical to the old
        // all-shards-resident pass.
        for (const auto& seg : ilen_segs) {
            seg_lo = seg.lo;
            seg_hi = seg.hi;
            seg_chunks = (seg_hi - seg_lo + kChunkSize - 1) / kChunkSize;
            const bool has_shard = seg.ilen >= kMinInteriorLen;
            if (has_shard) {
                const LenShard& ls_mem = build_shard(seg);
                uint64_t b = 0;
                for (const auto& pi : ls_mem.pi)
                    b += static_cast<uint64_t>(pi.keys.capacity())    * sizeof(*pi.keys.data()) +
                         static_cast<uint64_t>(pi.offsets.capacity()) * sizeof(*pi.offsets.data()) +
                         static_cast<uint64_t>(pi.ids.capacity())     * sizeof(*pi.ids.data()) +
                         static_cast<uint64_t>(pi.dir.capacity())     * sizeof(*pi.dir.data());
                fpi_peak_b = std::max(fpi_peak_b, b);
            }
            next_chunk.store(0, std::memory_order_relaxed);
            if (n_threads == 1) {
                b1_worker(0);
            } else {
                std::vector<std::thread> ts;
                ts.reserve(n_threads - 1);
                for (unsigned t = 1; t < n_threads; ++t) ts.emplace_back(b1_worker, t);
                b1_worker(0);
                for (auto& th : ts) th.join();
            }
            if (has_shard) shards.erase(seg.ilen);
        }

        if (errcor_.bucket_cap && errcor_.bucket_cap_lowcomplexity_only)
            log_info("Phase 3 entropy-cap: oversized buckets=" + std::to_string(ent_oversized)
                     + " capped(low-complexity)=" + std::to_string(ent_capped)
                     + " spared(high-complexity)=" + std::to_string(ent_spared)
                     + " postings_dropped=" + std::to_string(stats.bucket_overflow_drops));

        // Merge per-thread state. B2 sorts mismatches deterministically, but
        // we concat in tid order so the same -j N gives byte-identical output.
        size_t total_mm = 0;
        for (const auto& v : per_thread_mm) total_mm += v.size();
        std::vector<ChildMismatch> mismatches;
        mismatches.reserve(total_mm);
        for (auto& v : per_thread_mm) {
            mismatches.insert(mismatches.end(), v.begin(), v.end());
            std::vector<ChildMismatch>().swap(v);
        }
        for (const auto& ls : per_thread_stats) {
            stats.children_scanned       += ls.children_scanned;
            stats.total_candidates       += ls.total_candidates;
            stats.children_found         += ls.children_found;
            stats.decode_hash_child_ms   += ls.decode_hash_child_ms;
            stats.query_ms               += ls.query_ms;
            stats.check_ms               += ls.check_ms;
            short_interior_skipped       += ls.short_interior_skipped;
            stats.short_brute_evaluated  += ls.short_brute_evaluated;
            stats.short_brute_found      += ls.short_brute_found;
            stats.short_too_small_skipped += ls.short_too_small_skipped;
        }

        log_info("Phase 3 B1: " + std::to_string(n_threads) + " thread" +
                 (n_threads == 1 ? "" : "s") + ", found " +
                 std::to_string(stats.children_found) +
                 " directed edges from " + std::to_string(stats.children_scanned) + " sequences scanned");

        // ── Phase B2: directed ascent — sort high-count parents first, apply SNP veto ──
        // Processing highest-count parents first enables chain rerouting: when a parent
        // is itself absorbed (early PCR error), its children are redirected to the root's
        // count for the SNP veto.  This handles H=2 transitivity: A→B→C absorbs A even
        // when A is H=2 from C, because A→B uses B's count and B→C uses C's count.
        std::vector<uint32_t> parent_chain(N, UINT32_MAX);
        auto find_root_chain = [&](uint32_t id) -> uint32_t {
            while (parent_chain[id] != UINT32_MAX) id = parent_chain[id];
            return id;
        };

        std::sort(mismatches.begin(), mismatches.end(),
                  [&id_count](const ChildMismatch& a, const ChildMismatch& b) {
                      uint64_t ca = id_count[a.parent_id], cb = id_count[b.parent_id];
                      if (ca != cb) return ca > cb;  // descending count: absorb into highest count first
                      if (a.parent_id != b.parent_id) return a.parent_id < b.parent_id;
                      if (a.mismatch_pos != b.mismatch_pos) return a.mismatch_pos < b.mismatch_pos;
                      return a.alt_base < b.alt_base;
                  });

        // ── T5.8: empirical posterior-odds model ────────────────────────────
        // Build EdgeCandidate features for every B1 edge and fit the model
        // (per-bin P_real from cross-bundle recurrence; per-occ-bin
        // log(π_pcr/π_real) from singleton-vs-recurring split). The fitted
        // model replaces lr_threshold + kP_real_uniform + kBundlePrior_alpha.
        constexpr int kDamageEdgeMargin_T58 = 5;
        auto build_edge = [&](const ChildMismatch& cm) -> fqdup::errcor_emp::EdgeCandidate {
            const int L  = static_cast<int>(arena_.length(cm.child_id));
            const int k5 = get_layout(L).k5;
            const int ilen = get_layout(L).ilen;
            fqdup::errcor_emp::EdgeCandidate e{};
            e.child_id   = cm.child_id;
            e.parent_id  = cm.parent_id;
            e.bundle_key = bundle_key_of[cm.parent_id];
            e.bundle_occ = bundle_occ_of(cm.parent_id);
            e.child_len  = static_cast<uint16_t>(L);
            e.n_mm       = cm.hamming;
            auto fill = [&](int i, uint16_t mm_pos, uint8_t alt, uint8_t pb) {
                uint16_t pos_full = static_cast<uint16_t>(mm_pos + k5);
                e.mm[i].pos        = pos_full;
                e.mm[i].parent_2b  = pb;
                e.mm[i].alt_2b     = alt;
                e.mm[i].qual       = qual_arena_.q_at(cm.child_id, pos_full);
                const uint32_t dmg_xr = alt ^ pb;
                const bool is_ct_e  = (dmg_xr == 2u) && ((alt & pb) == 1u);
                const bool is_ga_e  = (dmg_xr == 2u) && ((alt & pb) == 0u);
                const bool at_5_e   = mm_pos < kDamageEdgeMargin_T58;
                const bool at_3_e   = mm_pos >= static_cast<uint16_t>(ilen - kDamageEdgeMargin_T58);
                e.mm[i].damage_chan = (profile_.ss_mode
                    ? (is_ct_e && (at_5_e || at_3_e))
                    : ((is_ct_e && at_5_e) || (is_ga_e && at_3_e))) ? 1 : 0;
                e.mm[i].p_damage    = (dmg_xr && profile_.enabled)
                                      ? profile_.p_damage_at(pos_full, L) : 0.0;
            };
            fill(0, cm.mismatch_pos, cm.alt_base, cm.parent_base);
            if (cm.hamming == 2)
                fill(1, cm.mismatch_pos2, cm.alt_base2, cm.parent_base2);
            return e;
        };

        // Cross-bundle recurrence map for interior-transition candidates.
        // Key: mismatch_pos * 4 + alt_base. Value: number of distinct bundle_keys
        // that have an interior-transition H=1 mismatch at (pos, alt).
        // Built once before B2; used in LRT path to veto absorption when ≥2
        // independent fragments (bundles) share the same interior mismatch.
        std::unordered_map<uint32_t, uint32_t> interior_trans_bundle_count;
        uint64_t n_interior_ct = 0;
        uint64_t n_interior_tv = 0;
        if (errcor_.b1_damage_adjust) {
            std::unordered_map<uint32_t, std::unordered_set<uint64_t>> key_to_bundles;
            for (const auto& cm : mismatches) {
                if (cm.hamming != 1) continue;
                const int pid_len = get_layout(static_cast<int>(arena_.length(cm.parent_id))).ilen;
                if (is_terminal_damage_mismatch(cm.parent_base, cm.alt_base,
                                                static_cast<int>(cm.mismatch_pos),
                                                pid_len, profile_.ss_mode))
                    continue;
                const uint32_t xr = cm.alt_base ^ cm.parent_base;
                if (xr == 2u) {
                    ++n_interior_ct;
                    uint32_t key = static_cast<uint32_t>(cm.mismatch_pos) * 4u + cm.alt_base;
                    key_to_bundles[key].insert(bundle_key_of[cm.parent_id]);
                } else {
                    ++n_interior_tv;
                }
            }
            interior_trans_bundle_count.reserve(key_to_bundles.size());
            for (const auto& [k, bset] : key_to_bundles)
                interior_trans_bundle_count[k] = static_cast<uint32_t>(bset.size());
        }
        // κ = interior C↔T rate / interior transversion rate.
        // κ≈1: no excess above PCR error floor → absorption would destroy real variants.
        // κ≥b1_kappa_min: genuine interior deamination signal → absorption is safe.
        double kappa = (n_interior_tv > 0)
            ? static_cast<double>(n_interior_ct) / static_cast<double>(n_interior_tv)
            : (n_interior_ct > 0 ? std::numeric_limits<double>::infinity() : 0.0);
        bool b1_active = errcor_.b1_damage_adjust && (kappa >= errcor_.b1_kappa_min);
        {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%.3f", kappa);
            std::string verdict = !errcor_.b1_damage_adjust ? "" :
                (b1_active ? " → active" : " → suppressed (κ < " +
                    std::to_string(errcor_.b1_kappa_min).substr(0, 4) + ")");
            log_info("Phase 3: interior κ=" + std::string(buf) +
                     " (CT=" + std::to_string(n_interior_ct) +
                     " TV=" + std::to_string(n_interior_tv) + ")" + verdict);
        }

        fqdup::errcor_emp::ErrCorEmpiricalModel emp_model;
        if (errcor_.empirical && !errcor_.legacy_veto && !mismatches.empty()) {
            // WIN 1: stream edges into the empirical model instead of
            // materializing the full all_edges vector (~43GB at full scale).
            // Each thread folds its slice of mismatches[] into a private per-bin
            // accumulator; the set-union merge is order-independent, so the
            // fitted model is bit-identical to the old build-then-fit path.
            // build_edge() is a pure function of one ChildMismatch, so B2
            // recomputes each edge inline where it needs it.
            unsigned n_edge_threads = errcor_.threads;
            if (n_edge_threads == 0)
                n_edge_threads = std::max(1u, std::thread::hardware_concurrency());
            const size_t NE = mismatches.size();
            std::vector<std::unordered_set<uint64_t>> local_bundles(n_edge_threads);
            std::vector<std::vector<std::unordered_set<uint64_t>>> local_bins(n_edge_threads);
            auto edge_worker = [&](unsigned t) {
                const size_t lo = NE * t / n_edge_threads;
                const size_t hi = NE * (t + 1) / n_edge_threads;
                auto& lb = local_bundles[t];
                auto& bins = local_bins[t];
                bins.assign(fqdup::errcor_emp::kNumBins, {});
                for (size_t j = lo; j < hi; ++j) {
                    emp_model.accumulate_edge(bins, build_edge(mismatches[j]));
                    lb.insert(bundle_key_of[mismatches[j].parent_id]);
                }
            };
            if (n_edge_threads == 1) {
                edge_worker(0);
            } else {
                std::vector<std::thread> ets;
                ets.reserve(n_edge_threads - 1);
                for (unsigned t = 1; t < n_edge_threads; ++t)
                    ets.emplace_back(edge_worker, t);
                edge_worker(0);
                for (auto& th : ets) th.join();
            }
            std::unordered_set<uint64_t> distinct_bundles;
            for (auto& lb : local_bundles)
                distinct_bundles.insert(lb.begin(), lb.end());
            // ── Phase-3 peak memory accounting (LOG-ONLY; byte-identical) ─────
            // Emitted at the Phase-3 resident peak: every large B1 structure is
            // still live (read arena, FlatPairIndex shards, interior slab,
            // mismatches[]) alongside the just-built distinct-bundle sets and the
            // per-thread fit accumulators, so accounted bytes reconcile against
            // the process VmHWM. Reads only container sizes — no algorithm or
            // output change.
            {
                auto vcap = [](const auto& v) -> uint64_t {
                    return static_cast<uint64_t>(v.capacity()) * sizeof(*v.data());
                };
                auto set_bytes = [](const std::unordered_set<uint64_t>& s) -> uint64_t {
                    // node ≈ 8B key + 8B next-ptr + malloc overhead (~16B); plus
                    // the bucket pointer array. Approximate, labelled "~".
                    return static_cast<uint64_t>(s.bucket_count()) * 8ull +
                           static_cast<uint64_t>(s.size()) * 32ull;
                };
                // Block-chunked (arena.hpp): allocated == resident, and there is no
                // doubling transient hiding on top of these numbers.
                uint64_t arena_b = arena_.packed.bytes_alloc() +
                                   arena_.offsets.bytes_alloc() +
                                   arena_.lengths.bytes_alloc() +
                                   arena_.eligible.bytes_alloc();
                uint64_t qual_b  = qual_arena_.bytes.bytes_alloc() +
                                   qual_arena_.offsets.bytes_alloc() +
                                   qual_arena_.lengths.bytes_alloc();
                uint64_t counts_b = vcap(id_count) + vcap(acc_count) +
                                    vcap(bundle_key_of) + vcap(rescue_bundle_key_of);
                // The dedup index is the largest single structure in the run and was
                // absent from this breakdown until now, leaving ~38% of VmHWM
                // unattributed. Slot is sherwood_v3_entry, not the bare pair.
                using IndexSlot = ska::detailv3::sherwood_v3_entry<
                    std::pair<SequenceFingerprint, IndexEntry>>;
                uint64_t index_b = static_cast<uint64_t>(index_.bucket_count()) *
                                   sizeof(IndexSlot);
                // Streamed: only one shard is ever live, so report the largest one
                // (sampled at build time) rather than the live sum, which is 0 here.
                uint64_t fpi_b  = fpi_peak_b;
                uint64_t mm_b   = vcap(mismatches);
                uint64_t slab_b = vcap(interior_slab) + vcap(interior_off);
                uint64_t db_b   = set_bytes(distinct_bundles);
                uint64_t bins_b = 0;
                for (const auto& lbset : local_bundles) bins_b += set_bytes(lbset);
                for (const auto& bins : local_bins)
                    for (const auto& s : bins) bins_b += set_bytes(s);
                // The bundle occupancy maps were the last blind spot: reserve(N)
                // took 24 GB each and never appeared here, so 48 GB of the peak
                // was attributed to a "doubling transient" that did not exist.
                using OccSlot = ska::detailv3::sherwood_v3_entry<
                    std::pair<uint64_t, uint32_t>>;
                uint64_t occ_b = static_cast<uint64_t>(
                                     bundle_occ_map.bucket_count() +
                                     rescue_bundle_occ_map.bucket_count()) *
                                 sizeof(OccSlot);
                // Streamed per-ilen like FlatPairIndex, so report the largest single
                // shard's postings. This term used to be built for ALL ilens at once
                // (33.65 GB) and was drained before this logger ran -- it was the
                // entire 39 GB residual between TOTAL_accounted and VmHWM.
                uint64_t post_b = post_peak_b;
                uint64_t total_b = index_b + arena_b + qual_b + counts_b + fpi_b +
                                   mm_b + slab_b + db_b + bins_b + occ_b + post_b;
                unsigned long vmhwm_kb = 0;
                if (FILE* f = std::fopen("/proc/self/status", "r")) {
                    char line[256];
                    while (std::fgets(line, sizeof(line), f))
                        if (std::sscanf(line, "VmHWM: %lu kB", &vmhwm_kb) == 1) break;
                    std::fclose(f);
                }
                auto gb = [](uint64_t b) {
                    char x[32];
                    std::snprintf(x, sizeof(x), "%.2f", static_cast<double>(b) / 1073741824.0);
                    return std::string(x);
                };
                auto fld = [&](const char* n, uint64_t b) {
                    return std::string(n) + "=" + std::to_string(b) + "|" + gb(b) + "GB ";
                };
                char lf[64];
                std::snprintf(lf, sizeof(lf), "index_[%zu slots @load %.2f] ",
                              index_.bucket_count(),
                              index_.size() / static_cast<double>(
                                  index_.bucket_count() ? index_.bucket_count() : 1));
                log_info("Phase 3 mem breakdown: " + std::string(lf) +
                         fld("index", index_b) +
                         fld("read_arena", arena_b) + fld("qual_arena", qual_b) +
                         fld("counts(id+acc+bundle_key)", counts_b) +
                         fld("FlatPairIndex_peakshard", fpi_b) +
                         fld("postings_peakshard", post_b) + fld("mismatches", mm_b) +
                         fld("interior_slab", slab_b) +
                         fld("bundle_occ_maps", occ_b) +
                         "[keys=" + std::to_string(bundle_occ_map.size()) + "+" +
                         std::to_string(rescue_bundle_occ_map.size()) + " of N=" +
                         std::to_string(N) + "] " +
                         fld("distinct_bundles~", db_b) + fld("fit_accum~", bins_b) +
                         "| TOTAL_accounted=" + std::to_string(total_b) + "|" + gb(total_b) + "GB" +
                         "  VmHWM=" + std::to_string(static_cast<uint64_t>(vmhwm_kb) * 1024ull) +
                         "|" + gb(static_cast<uint64_t>(vmhwm_kb) * 1024ull) + "GB");
            }
            for (auto& bins : local_bins)
                emp_model.fit_merge(bins);
            emp_model.fit_finalize(distinct_bundles.size());
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.4g", emp_model.p_real_global);
            log_info("Phase 3 empirical model: " + std::to_string(NE) +
                     " edges, " + std::to_string(distinct_bundles.size()) +
                     " distinct bundles, p_real_global=" + buf);
            std::snprintf(buf, sizeof(buf), "%.3f", emp_model.log_pi_ratio_global);
            log_info("  log_pi_ratio_global=" + std::string(buf) + " nats; per-occ:");
            for (int o = 0; o < fqdup::errcor_emp::kNumOccBins; ++o) {
                std::snprintf(buf, sizeof(buf), "%.3f", emp_model.log_pi_ratio_occ[o]);
                log_info("    occ_bin[" + std::to_string(o) + "] = " + buf);
            }
        }
        auto joint_dispersion_adj = [](
            const std::vector<std::pair<uint32_t,uint64_t>>& pac,
            uint64_t veto_count,
            double snp_threshold) -> double
        {
            if (pac.empty() || veto_count < 4) return 0.0;
            uint64_t total_w = 0, max_w = 0;
            for (const auto& [k, w] : pac) { total_w += w; max_w = std::max(max_w, w); }
            if (total_w == 0) return 0.0;
            double max_frac = static_cast<double>(max_w) / static_cast<double>(total_w);
            if (max_frac >= snp_threshold) return 0.0;
            if (max_frac > 0.10)
                return -1.5 * (max_frac - 0.10) / (snp_threshold - 0.10);
            if (pac.size() >= 3) return 0.3;
            return 0.0;
        };

        size_t i = 0;
        uint64_t n_joint_adj_pos = 0;
        uint64_t n_joint_adj_neg = 0;
        uint64_t n_parents_processed = 0;
        while (i < mismatches.size()) {
            uint32_t pid = mismatches[i].parent_id;
            // Reroute to effective root if this parent was itself absorbed.
            uint32_t eff_pid = pid;
            if (is_error_[pid]) eff_pid = find_root_chain(pid);
            uint64_t parent_count = id_count[eff_pid];
            // SNP veto denominator uses pid's original count: sig was accumulated from
            // pid's children, so the correct ratio is sig/pid_count, not sig/eff_pid_count.
            // Using eff_pid's inflated count deflates the ratio and causes over-absorption
            // of real SNPs that cluster under an absorbed intermediate parent.
            uint64_t veto_count = id_count[pid];

            // Collect all children for this parent
            size_t parent_start = i;
            while (i < mismatches.size() && mismatches[i].parent_id == pid) ++i;
            size_t parent_end = i;

            // Accumulate sig_count per (pos, alt_base) for SNP detection.
            // Key: mismatch_pos * 4 + alt_base. Sort+collapse a small vector
            // (typical fanout ≤ 32) — beats unordered_map allocations + cache misses.
            std::vector<std::pair<uint32_t, uint64_t>> pos_alt_counts;
            pos_alt_counts.reserve(2 * (parent_end - parent_start));
            for (size_t j = parent_start; j < parent_end; ++j) {
                uint64_t cc = id_count[mismatches[j].child_id];
                pos_alt_counts.emplace_back(
                    static_cast<uint32_t>(mismatches[j].mismatch_pos) * 4u
                        + mismatches[j].alt_base,
                    cc);
                if (mismatches[j].hamming == 2) {
                    pos_alt_counts.emplace_back(
                        static_cast<uint32_t>(mismatches[j].mismatch_pos2) * 4u
                            + mismatches[j].alt_base2,
                        cc);
                }
            }
            std::sort(pos_alt_counts.begin(), pos_alt_counts.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
            // Collapse runs of equal key in-place.
            size_t pac_w = 0;
            for (size_t r = 0; r < pos_alt_counts.size(); ++r) {
                if (pac_w > 0 && pos_alt_counts[pac_w - 1].first == pos_alt_counts[r].first)
                    pos_alt_counts[pac_w - 1].second += pos_alt_counts[r].second;
                else
                    pos_alt_counts[pac_w++] = pos_alt_counts[r];
            }
            pos_alt_counts.resize(pac_w);
            double joint_adj = joint_dispersion_adj(pos_alt_counts, veto_count, errcor_.snp_threshold);
            ++n_parents_processed;
            if (joint_adj > 0.0)  ++n_joint_adj_pos;
            else if (joint_adj < 0.0) ++n_joint_adj_neg;
            auto find_pac = [&](uint32_t k) -> uint64_t {
                auto it = std::lower_bound(
                    pos_alt_counts.begin(), pos_alt_counts.end(), k,
                    [](const std::pair<uint32_t, uint64_t>& p, uint32_t v) {
                        return p.first < v;
                    });
                return (it != pos_alt_counts.end() && it->first == k) ? it->second : 0ull;
            };

            // For each child, decide: SNP-protected or absorb
            const int pid_ilen = get_layout(static_cast<int>(arena_.length(pid))).ilen;
            for (size_t j = parent_start; j < parent_end; ++j) {
                const ChildMismatch& cm = mismatches[j];
                uint32_t key = static_cast<uint32_t>(cm.mismatch_pos) * 4u + cm.alt_base;
                uint64_t sig = find_pac(key);

                // H=2 path: both mismatches already verified non-damage in B1.
                // SNP veto: protect if either mismatch position has population support.
                if (cm.hamming == 2) {
                    uint32_t key2 = static_cast<uint32_t>(cm.mismatch_pos2) * 4u + cm.alt_base2;
                    uint64_t sig2 = find_pac(key2);

                    double eff_snp = errcor_.snp_threshold *
                        (veto_count < errcor_.snp_low_cov_cutoff ? errcor_.snp_low_cov_factor : 1.0);
                    // On strict-dominance edges relax snp_min_count to 1; keep 2 on equal-count
                    // edges where singleton-vs-singleton false protection is a real risk.
                    uint32_t eff_min_count = (parent_count > id_count[cm.child_id]) ? 1u : errcor_.snp_min_count;
                    bool snp_veto_h2 = (sig  >= eff_min_count &&
                                        static_cast<double>(sig)  >= eff_snp * static_cast<double>(veto_count)) ||
                                       (sig2 >= eff_min_count &&
                                        static_cast<double>(sig2) >= eff_snp * static_cast<double>(veto_count));
                    // T5.8: empirical posterior-odds for H=2.
                bool   absorb_h2;
                float  fqcl_lr_h2 = std::numeric_limits<float>::quiet_NaN();
                if (errcor_.empirical && !errcor_.legacy_veto) {
                    auto e = build_edge(cm);
                    e.bundle_key = bundle_key_of[eff_pid];
                    e.bundle_occ = bundle_occ_of(eff_pid);
                    double S = emp_model.score(e);
                    // Count-asymmetry prior: PCR error families have
                    // parent_count >> child_count (the error is rare per
                    // amplification, so the descendant lineage is small).
                    // Real biology has parent_count ~ child_count (two
                    // independent molecules drawn from the same coverage).
                    // log((p+1)/(c+1)) shifts the log-odds toward PCR when
                    // the parent dominates. Scale by n_mm: each independent
                    // PCR error compounds the parent-vs-child likelihood
                    // ratio (probability of N errors going the same direction
                    // is the product), so the H=2 hypothesis needs 2× the
                    // count-asymmetry support that H=1 needs to overcome the
                    // 2× per-mismatch penalty.
                    S += static_cast<double>(e.n_mm) *
                         std::log(static_cast<double>(parent_count + 1) /
                                  static_cast<double>(id_count[cm.child_id] + 1));
                    S += joint_adj;
                    absorb_h2 = (S > 0.0);
                    fqcl_lr_h2 = static_cast<float>(S);
                    // H=3 shadow: track near-boundary LR scores before veto logic.
                    constexpr double kH3ShadowNat = 1.0;
                    if      (S >  0.0 && S <  kH3ShadowNat) ++stats.h3_shadow_absorbed_near_boundary;
                    else if (S <= 0.0 && S > -kH3ShadowNat) ++stats.h3_shadow_protected_near_boundary;
                    // Count-based co-occurrence clamp. The ratio test in
                    // snp_veto_h2 fails at hyper-cluster scale (sig must clear
                    // ~snp_threshold * inflated parent_count). When ≥K
                    // independent reads share the same (pos, alt), that is
                    // direct population evidence of a real SNP — block
                    // absorption regardless of S. K = max(snp_min_count, 2):
                    // sig==1 is the child alone; sig≥2 is genuine co-occurrence.
                    // H=2 mismatches were pre-filtered against damage in B1, so
                    // no damage_bypass exemption applies here.
                    uint32_t hard_veto_min_h2 = std::max(errcor_.snp_min_count, 2u);
                    // Singleton high-Q clamp. H=2 mismatches were pre-filtered
                    // against damage in B1, so any high-Q mismatch is variant
                    // evidence regardless of co-occurrence count.
                    bool hi_q_singleton_h2 =
                        (e.mm[0].qual >= errcor_.singleton_qual_min) ||
                        (e.mm[1].qual >= errcor_.singleton_qual_min);
                    if (absorb_h2 && (sig >= hard_veto_min_h2 ||
                                      sig2 >= hard_veto_min_h2 ||
                                      hi_q_singleton_h2))
                        absorb_h2 = false;
                    if (snp_veto_h2) { stats.lr_sum_protected += S; ++stats.lr_n_protected; }
                    else             { stats.lr_sum_absorbed  += S; ++stats.lr_n_absorbed;  }
                    if (absorb_h2) ++stats.edge_lr_absorbed;
                    else           ++stats.edge_lr_protected;
                    int ob = fqdup::errcor_emp::occ_bin(e.bundle_occ);
                    if (absorb_h2) ++stats.edge_lr_absorbed_by_occ[ob];
                    else           ++stats.edge_lr_protected_by_occ[ob];
                } else {
                    absorb_h2 = !snp_veto_h2;
                    if (snp_veto_h2) ++stats.edge_legacy_veto;
                    else             ++stats.edge_legacy_absorb;
                }

                if (!absorb_h2) {
                        stats.snp_protected++;
                    } else {
                        if (!is_error_[cm.child_id]) {
                            is_error_[cm.child_id] = true;
                            parent_chain[cm.child_id] = eff_pid;
                            acc_count[eff_pid] += acc_count[cm.child_id];
                            acc_count[cm.child_id] = 0;
                            if (seq_entry[cm.child_id] && seq_entry[eff_pid] &&
                                seq_entry[cm.child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                                seq_entry[eff_pid]->record_index = seq_entry[cm.child_id]->record_index;
                                seq_entry[eff_pid]->damage_score = seq_entry[cm.child_id]->damage_score;
                            }
                            stats.absorbed++;
                            errcor_absorbed_++;
                            if (!fqcl_path_.empty()) {
                                auto cm_s = cm; cm_s.lr_score = fqcl_lr_h2;
                                fqcl_mismatches_.push_back(cm_s);
                            }
                        }
                    }
                    continue;  // skip the H=1 block below
                }

                // Damage-aware bypass: residual deamination at terminal positions.
                // DS libraries: C→T at 5' end, G→A at 3' end.
                // SS libraries: C→T at both ends (G→A is not a damage pattern in SS).
                // Distinguish C↔T (xr=2, AND=1) from G↔A (xr=2, AND=0).
                const uint32_t dmg_xr_h1 = cm.alt_base ^ cm.parent_base;
                const bool is_ct_h1 = (dmg_xr_h1 == 2u) && ((cm.alt_base & cm.parent_base) == 1u);
                const bool is_ga_h1 = (dmg_xr_h1 == 2u) && ((cm.alt_base & cm.parent_base) == 0u);
                const bool at_5_h1  = cm.mismatch_pos < kDamageEdgeMargin;
                const bool at_3_h1  = cm.mismatch_pos >=
                                          static_cast<uint16_t>(pid_ilen - kDamageEdgeMargin);
                bool damage_bypass = profile_.enabled && (
                    profile_.ss_mode
                    ? (is_ct_h1 && (at_5_h1 || at_3_h1))
                    : ((is_ct_h1 && at_5_h1) || (is_ga_h1 && at_3_h1)));

                double eff_snp_threshold = errcor_.snp_threshold *
                    (veto_count < errcor_.snp_low_cov_cutoff ? errcor_.snp_low_cov_factor : 1.0);
                // On strict-dominance edges relax snp_min_count to 1.
                // damage_bypass mismatches are still absorbed unless population support is strong
                // enough (sig >= 1 AND ratio >= threshold) indicating a real SNP at the margin.
                uint32_t eff_min_count = (parent_count > id_count[cm.child_id]) ? 1u : errcor_.snp_min_count;
                bool snp_veto = (sig >= (damage_bypass ? 1u : eff_min_count)) &&
                                (static_cast<double>(sig) >= eff_snp_threshold *
                                                              static_cast<double>(veto_count));
                if (damage_bypass) ++stats.edge_damage_bypass;

                // Interior transition: C↔T or G↔A at a non-terminal position.
                // Admitted by b1_damage_adjust but requires calibration-free LRT
                // (not the empirical model, which was fit on transversion candidates).
                const bool is_interior_trans_h1 = !damage_bypass && (dmg_xr_h1 == 2u) && b1_active;

                // T5.8: empirical posterior-odds for H=1.
                bool absorb_h1;
                float fqcl_lr_h1 = std::numeric_limits<float>::quiet_NaN();
                if (errcor_.empirical && !errcor_.legacy_veto) {
                    auto e = build_edge(cm);
                    e.bundle_key = bundle_key_of[eff_pid];
                    e.bundle_occ = bundle_occ_of(eff_pid);

                    if (is_interior_trans_h1) {
                        // === Calibration-free LRT path for interior transitions ===
                        uint64_t n_child_count = id_count[cm.child_id];
                        double log_lambda = interior_transition_lrt(
                            parent_count, n_child_count, errcor_.lrt_f0, errcor_.lrt_f1);
                        bool lrt_absorb = (log_lambda > std::log(errcor_.lrt_T));
                        bool cross_strand_veto = false;
                        bool hard_absorb = (parent_count >
                                            static_cast<uint64_t>(errcor_.lrt_hard_ratio * n_child_count))
                                        && (e.mm[0].qual < errcor_.singleton_qual_min);
                        // Wilson CI SNP veto: require sig ≥ 2 (same as empirical hard_veto_min).
                        double wlb = wilson_lower95(sig, veto_count);
                        bool wilson_veto = (sig >= 2 && wlb >= errcor_.snp_cp_lb_threshold);
                        // Cross-bundle recurrence veto: if ≥2 independent bundles show
                        // this same (pos, alt) interior mismatch it is biological, not PCR.
                        auto cb_it = interior_trans_bundle_count.find(key);
                        bool cross_bundle_veto = (cb_it != interior_trans_bundle_count.end()
                                                  && cb_it->second >= 2);
                        // Cross-strand veto: ancient deamination is strand-specific; a
                        // deamination artifact child has reads from one orientation only.
                        // A true SNP child has reads from both fwd and RC orientations.
                        // Require ≥b1_cs_min_total total child reads to avoid false vetoes
                        // where orientation balance is random at very low coverage.
                        {
                            uint64_t cf = id_fwd_count[cm.child_id];
                            uint64_t ct = id_count[cm.child_id];
                            uint64_t cr = (ct >= cf) ? ct - cf : 0;
                            cross_strand_veto = (cf > 0 && cr > 0 &&
                                                 ct >= errcor_.b1_cs_min_total);
                        }
                        bool hi_q = (e.mm[0].qual >= errcor_.singleton_qual_min);
                        absorb_h1 = (hard_absorb || lrt_absorb)
                                    && !wilson_veto && !cross_bundle_veto
                                    && !cross_strand_veto && !hi_q;
                        ++stats.edge_lr_absorbed;  // reuse counter for logging
                    } else {
                    double S = emp_model.score(e);
                    // See H=2 site for rationale on the count-asymmetry term.
                    S += std::log(static_cast<double>(parent_count + 1) /
                                  static_cast<double>(id_count[cm.child_id] + 1));
                    S += joint_adj;
                    absorb_h1 = (S > 0.0);
                    fqcl_lr_h1 = static_cast<float>(S);
                    // Count-based co-occurrence clamp. snp_veto's ratio test
                    // fails at hyper-cluster scale (sig must clear ~snp_threshold
                    // * inflated parent_count). When ≥K independent reads share
                    // the same (pos, alt), that is direct population evidence of
                    // a real SNP — block absorption regardless of S.
                    // K = max(snp_min_count, 2): sig==1 is the child alone;
                    // sig≥2 is genuine co-occurrence. damage_bypass mismatches
                    // at terminal positions remain exempt (residual deamination).
                    uint32_t hard_veto_min = std::max(errcor_.snp_min_count, 2u);
                    // Singleton high-Q clamp. A single high-Q non-damage
                    // mismatch is direct sequence-divergence evidence; the
                    // empirical S>0 path otherwise allows absorption of
                    // singleton variant alleles at low-coverage loci.
                    // damage_bypass remains exempt (residual deamination).
                    bool hi_q_singleton =
                        (e.mm[0].qual >= errcor_.singleton_qual_min);
                    if (absorb_h1 && !damage_bypass &&
                        (sig >= hard_veto_min || hi_q_singleton))
                        absorb_h1 = false;
                    if (snp_veto) { stats.lr_sum_protected += S; ++stats.lr_n_protected; }
                    else          { stats.lr_sum_absorbed  += S; ++stats.lr_n_absorbed;  }
                    if (absorb_h1) ++stats.edge_lr_absorbed;
                    else           ++stats.edge_lr_protected;
                    int ob = fqdup::errcor_emp::occ_bin(e.bundle_occ);
                    if (absorb_h1) ++stats.edge_lr_absorbed_by_occ[ob];
                    else           ++stats.edge_lr_protected_by_occ[ob];
                    }  // end empirical path (non-interior-transition)
                } else {
                    absorb_h1 = !snp_veto;
                    if (snp_veto) ++stats.edge_legacy_veto;
                    else          ++stats.edge_legacy_absorb;
                }

                if (!absorb_h1) {
                    stats.snp_protected++;
                } else {
                    if (!is_error_[cm.child_id]) {
                        is_error_[cm.child_id] = true;
                        parent_chain[cm.child_id] = eff_pid;  // track chain for rerouting descendants
                        acc_count[eff_pid] += acc_count[cm.child_id];  // propagate accumulated mass (incl. child's absorbed descendants)
                        acc_count[cm.child_id] = 0;  // zero to prevent double-counting if re-traversed
                        // Propagate representative: adopt child's record if it has stronger damage signal.
                        if (seq_entry[cm.child_id] && seq_entry[eff_pid] &&
                            seq_entry[cm.child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                            seq_entry[eff_pid]->record_index = seq_entry[cm.child_id]->record_index;
                            seq_entry[eff_pid]->damage_score = seq_entry[cm.child_id]->damage_score;
                        }
                        stats.absorbed++;
                        errcor_absorbed_++;
                        if (!fqcl_path_.empty()) {
                            auto cm_s = cm; cm_s.lr_score = fqcl_lr_h1;
                            fqcl_mismatches_.push_back(cm_s);
                        }
                    }
                }
            }
        }

        // ── Phase B3: damage-aware H>2 merge ───────────────────────────────
        // For heavily damaged ancient DNA (d_max ≥ b3_min_dmax), PCR copies of
        // the same molecule can accumulate 3–5 deamination differences in the
        // near-terminal interior zone (just beyond the mask). Phase 3 pigeonhole
        // only covers H≤2 — this pass fills that gap.
        // Algorithm: K_deam hash buckets (normalize damage-probable positions)
        // → same-length pairs within bucket → all-deam-consistent diff check
        // → count LRT (same threshold as B2).
        if (errcor_.b3_enabled && profile_.enabled && !errcor_.legacy_veto) {
            // Profile-derived gate: count interior positions where P(damage) >
            // b3_deam_threshold using the median read length from the arena.
            // Skip if n_elig < b3_min_n_elig (B3 cannot produce H>=3 absorptions)
            // or sum(P) < b3_min_mass (events too rare to matter).
            int rep_L = 0;
            {
                std::vector<uint16_t> ls(arena_.lengths.size());
                for (size_t i = 0; i < ls.size(); ++i) ls[i] = arena_.lengths[i];
                if (!ls.empty()) {
                    std::nth_element(ls.begin(), ls.begin() + ls.size() / 2, ls.end());
                    rep_L = static_cast<int>(ls[ls.size() / 2]);
                }
            }
            int    b3_n_elig = 0;
            double b3_mass   = 0.0;
            if (rep_L > 0) {
                auto [k5, k3] = damage_zone_bounds(rep_L, profile_);
                int interior_end = rep_L - k3;
                for (int p = k5; p < interior_end; ++p) {
                    double q = profile_.p_damage_at(p, rep_L);
                    if (q > static_cast<double>(errcor_.b3_deam_threshold)) {
                        ++b3_n_elig;
                        b3_mass += q;
                    }
                }
            }
            double dmax_avg = 0.5 * (profile_.d_max_5prime + profile_.d_max_3prime);
            if (b3_n_elig < errcor_.b3_min_n_elig ||
                b3_mass   < static_cast<double>(errcor_.b3_min_mass)) {
                log_info("Phase B3: skipped (n_elig=" + std::to_string(b3_n_elig) +
                         " mass=" + std::to_string(b3_mass) +
                         " d_max_avg=" + std::to_string(dmax_avg) + ")");
            } else {
                log_info("Phase B3: running (n_elig=" + std::to_string(b3_n_elig) +
                         " mass=" + std::to_string(b3_mass) +
                         " d_max_avg=" + std::to_string(dmax_avg) + ")");
                // Cross-bundle veto map: (pos_abs*4 + base) → count of distinct
                // bundle_keys that carry that (position, base) combination at a
                // damage-eligible interior position. If ≥2 independent loci share
                // the same mismatch, it is a real SNP, not PCR deamination noise.
                ska::flat_hash_map<uint32_t, uint32_t> b3_bundle_veto;
                {
                    ska::flat_hash_map<uint32_t, ska::flat_hash_set<uint64_t>> key_to_bundles;
                    for (uint32_t id = 0; id < N; ++id) {
                        if (is_error_[id]) continue;
                        if (!arena_.is_eligible(id)) continue;
                        int L = arena_.length(id);
                        const auto& lay = get_layout(L);
                        if (lay.ilen < kMinInteriorLen) continue;
                        const uint8_t* data = arena_.data(id);
                        uint64_t bkey = bundle_key_of[id];
                        for (int ip = 0; ip < lay.ilen; ++ip) {
                            int pos_abs = lay.k5 + ip;
                            if (profile_.p_damage_at(pos_abs, L) < static_cast<double>(errcor_.b3_deam_threshold))
                                continue;
                            uint8_t base = (data[pos_abs >> 2] >> (6 - 2 * (pos_abs & 3))) & 0x3u;
                            uint32_t key = static_cast<uint32_t>(pos_abs) * 4u + base;
                            key_to_bundles[key].insert(bkey);
                        }
                    }
                    b3_bundle_veto.reserve(key_to_bundles.size());
                    for (const auto& [k, bset] : key_to_bundles)
                        b3_bundle_veto[k] = static_cast<uint32_t>(bset.size());
                }
                // Build K_deam buckets.
                ska::flat_hash_map<uint64_t, std::vector<uint32_t>> b3_buckets;
                b3_buckets.reserve(N / 2);
                {
                    std::vector<uint8_t> scratch;
                    for (uint32_t id = 0; id < N; ++id) {
                        if (is_error_[id]) continue;
                        if (!arena_.is_eligible(id)) continue;
                        int L = arena_.length(id);
                        const auto& lay = get_layout(L);
                        if (lay.ilen < kMinInteriorLen) continue;
                        scratch.resize(static_cast<size_t>((lay.ilen + 3) / 4));
                        uint64_t kd = fqdup::b3::kdamage_hash(
                            arena_.data(id), lay.k5, lay.ilen, L,
                            profile_, errcor_.b3_deam_threshold, scratch.data());
                        // Gate by locus: combine bundle_key (same fragment endpoints)
                        // with kdamage_hash so only reads from the same capture site
                        // are ever compared.
                        uint64_t b3k = XXH3_64bits_withSeed(&kd, sizeof(kd), bundle_key_of[id]);
                        b3_buckets[b3k].push_back(id);
                    }
                }
                // Pairwise check within each bucket.
                std::vector<uint8_t> pi_buf, ci_buf;
                for (auto& [kd, bucket] : b3_buckets) {
                    if (bucket.size() < 2) continue;
                    std::sort(bucket.begin(), bucket.end(),
                              [&](uint32_t a, uint32_t b_) {
                                  return acc_count[a] > acc_count[b_];
                              });
                    for (size_t pidx = 0; pidx < bucket.size(); ++pidx) {
                        uint32_t pid = bucket[pidx];
                        if (is_error_[pid]) continue;
                        uint32_t eff_pid = find_root_chain(pid);
                        int Lp = arena_.length(eff_pid);
                        const auto& layp = get_layout(Lp);
                        if (layp.ilen < kMinInteriorLen) continue;
                        int nbp = (layp.ilen + 3) / 4;
                        pi_buf.resize(static_cast<size_t>(nbp));
                        extract_packed_part(arena_.data(eff_pid), layp.k5, layp.ilen,
                                            pi_buf.data());
                        uint64_t cp = acc_count[eff_pid];
                        for (size_t cidx = pidx + 1; cidx < bucket.size(); ++cidx) {
                            uint32_t cid = bucket[cidx];
                            if (is_error_[cid]) continue;
                            if (arena_.length(cid) != Lp) continue;
                            uint64_t cc = acc_count[cid];
                            if (cc == 0) continue;
                            // Require parent >= b3_count_ratio x child abundance
                            // (a real sequence dominates its deaminated sibling).
                            // The bucket is sorted count-descending so cp >= cc;
                            // the old cc/cp test could never exceed 1, so the gate
                            // never fired and no ratio filter was applied.
                            if (static_cast<double>(cp) <
                                static_cast<double>(cc) * errcor_.b3_count_ratio)
                                continue;
                            ++stats.b3_candidates;
                            ci_buf.resize(static_cast<size_t>(nbp));
                            extract_packed_part(arena_.data(cid), layp.k5, layp.ilen,
                                                ci_buf.data());
                            fqdup::b3::DiffList diffs;
                            fqdup::b3::packed_diff_list(pi_buf.data(), ci_buf.data(),
                                                         layp.ilen, diffs);
                            if (diffs.count < 3 || diffs.count > errcor_.b3_max_hamming)
                                continue;
                            if (!fqdup::b3::all_diffs_deam_consistent(
                                    diffs, layp.k5, Lp, profile_,
                                    errcor_.b3_deam_threshold))
                                continue;
                            // Cross-bundle veto: if any diff position+base is seen
                            // in ≥2 independent bundle_keys it is a real SNP, not
                            // a PCR deamination artifact — block absorption.
                            {
                                bool cross_veto = false;
                                for (int di = 0; di < diffs.count && !cross_veto; ++di) {
                                    uint32_t key = static_cast<uint32_t>(layp.k5 + diffs.pos[di]) * 4u
                                                   + diffs.base_c[di];
                                    auto it = b3_bundle_veto.find(key);
                                    if (it != b3_bundle_veto.end() && it->second >= 2)
                                        cross_veto = true;
                                }
                                if (cross_veto) {
                                    ++stats.b3_cross_bundle_protected;
                                    continue;
                                }
                            }
                            double lrt = fqdup::b3::b3_count_lrt(
                                cp, cc, diffs.count, errcor_.lrt_f0, errcor_.lrt_f1);
                            if (lrt <= std::log(errcor_.lrt_T)) {
                                ++stats.b3_protected;
                                continue;
                            }
                            if (!is_error_[cid]) {
                                is_error_[cid] = true;
                                parent_chain[cid] = eff_pid;
                                acc_count[eff_pid] += cc;
                                acc_count[cid] = 0;
                                if (seq_entry[cid] && seq_entry[eff_pid] &&
                                    seq_entry[cid]->damage_score >
                                    seq_entry[eff_pid]->damage_score) {
                                    seq_entry[eff_pid]->record_index =
                                        seq_entry[cid]->record_index;
                                    seq_entry[eff_pid]->damage_score =
                                        seq_entry[cid]->damage_score;
                                }
                                ++stats.b3_absorbed;
                                ++errcor_absorbed_;
                                if (!fqcl_path_.empty()) {
                                    ChildMismatch cm_b3{};
                                    cm_b3.parent_id  = eff_pid;
                                    cm_b3.child_id   = cid;
                                    cm_b3.hamming    = static_cast<uint8_t>(diffs.count);
                                    cm_b3.lr_score   = std::numeric_limits<float>::quiet_NaN();
                                    if (diffs.count >= 1) {
                                        cm_b3.mismatch_pos = static_cast<uint16_t>(layp.k5 + diffs.pos[0]);
                                        cm_b3.alt_base     = diffs.base_c[0];
                                        cm_b3.parent_base  = diffs.base_p[0];
                                    }
                                    if (diffs.count >= 2) {
                                        cm_b3.mismatch_pos2 = static_cast<uint16_t>(layp.k5 + diffs.pos[1]);
                                        cm_b3.alt_base2     = diffs.base_c[1];
                                        cm_b3.parent_base2  = diffs.base_p[1];
                                    }
                                    fqcl_mismatches_.push_back(cm_b3);
                                }
                            }
                        }
                    }
                }
                log_info("Phase B3: candidates=" + std::to_string(stats.b3_candidates) +
                         " absorbed=" + std::to_string(stats.b3_absorbed) +
                         " protected=" + std::to_string(stats.b3_protected) +
                         " cross_bundle_protected=" + std::to_string(stats.b3_cross_bundle_protected));
            }
        }

        // ── T5.6: Adjacent-length (L±1) probe ──────────────────────────────
        // Catches PCR siblings differing by a single insertion/deletion at any
        // position (damage-truncated reads, polymerase slippage). Restricted to
        // bundle members so candidate-set is locus-bounded — never N² over all
        // sequences. Per-pair cost: O(L) banded check; absorbs only on exact
        // 1-indel + 0-substitution match. LR rule applies as for H=1.
        if (errcor_.adj_len_probe && !errcor_.legacy_veto) {
            // Build a length-agnostic bundle map so L vs L±1 reads share keys.
            // Only constructed when the probe is enabled — zero overhead off-path.
            ska::flat_hash_map<uint64_t, std::vector<uint32_t>> adj_bundles;
            adj_bundles.reserve(N);
            {
                std::vector<uint8_t> dec;
                std::vector<char>    asc;
                for (uint32_t id = 0; id < N; ++id) {
                    if (!arena_.is_eligible(id)) continue;
                    if (!ec_eligible(id)) continue;   // no short-interior read enters adj-len absorption
                    int L = arena_.length(id);
                    if ((int)dec.size() < L) { dec.resize(L); asc.resize(L); }
                    arena_.decode_range(id, 0, L, dec.data());
                    for (int i = 0; i < L; ++i) asc[i] = static_cast<char>(dec[i]);
                    uint64_t k = fqdup::bundlekey::from_decoded_no_len(
                        asc.data(), L, fqdup::bundlekey::kDefaultEndK);
                    adj_bundles[k].push_back(id);
                }
            }

            // Per-bundle parent shortlist cap. On capture hotspots a bundle can
            // hold hundreds of clusters; pairing all-vs-all is O(M²). We restrict
            // candidate parents to the top-K members by id_count (the only ones
            // that can plausibly absorb others), then test every potential child
            // of adjacent length against those parents. Work is O(K·M_adj).
            constexpr size_t kAdjLenParentCap = 16;

            std::vector<uint8_t> abuf, bbuf;
            for (auto& kv : adj_bundles) {
                auto& members = kv.second;
                if (members.size() < 2) continue;
                std::vector<std::pair<int, uint32_t>> by_len;
                by_len.reserve(members.size());
                for (uint32_t mid : members)
                    by_len.emplace_back(static_cast<int>(arena_.length(mid)), mid);
                std::sort(by_len.begin(), by_len.end());
                const size_t M = by_len.size();

                // Build parent shortlist: top-K members by id_count.
                std::vector<uint32_t> parents;
                if (M <= kAdjLenParentCap) {
                    parents.reserve(M);
                    for (auto& p : by_len) parents.push_back(p.second);
                } else {
                    parents.reserve(M);
                    for (auto& p : by_len) parents.push_back(p.second);
                    std::partial_sort(
                        parents.begin(),
                        parents.begin() + kAdjLenParentCap,
                        parents.end(),
                        [&](uint32_t a, uint32_t b) {
                            if (id_count[a] != id_count[b]) return id_count[a] > id_count[b];
                            return a < b;  // stable tie-break for determinism
                        });
                    parents.resize(kAdjLenParentCap);
                }
                ska::flat_hash_set<uint32_t> parent_set(parents.begin(), parents.end());

                for (size_t i = 0; i < M; ++i) {
                    int Li = by_len[i].first;
                    for (size_t j = i + 1; j < M; ++j) {
                        int Lj = by_len[j].first;
                        if (Lj - Li > 1) break;
                        if (Lj - Li != 1) continue;  // only adjacent (skip same-length, handled by B2)
                        // At least one side must be in the parent shortlist —
                        // otherwise neither could plausibly be the absorber.
                        if (!parent_set.count(by_len[i].second) &&
                            !parent_set.count(by_len[j].second)) continue;
                        uint32_t a_id = by_len[i].second;
                        uint32_t b_id = by_len[j].second;  // longer
                        ++stats.adj_len_evaluated;

                        if ((int)abuf.size() < Li) abuf.resize(Li);
                        if ((int)bbuf.size() < Lj) bbuf.resize(Lj);
                        arena_.decode_range(a_id, 0, Li, abuf.data());
                        arena_.decode_range(b_id, 0, Lj, bbuf.data());

                        // Test: is `a` (shorter) exactly a 1-deletion of `b` (longer)?
                        int d = 0;
                        while (d < Li && bbuf[d] == abuf[d]) ++d;
                        bool ok = true;
                        for (int p = d; p < Li; ++p)
                            if (bbuf[p + 1] != abuf[p]) { ok = false; break; }
                        if (!ok) continue;
                        ++stats.adj_len_matched;

                        // Higher count wins parent role (matches B2 convention).
                        uint32_t parent_id, child_id;
                        if (id_count[a_id] >= id_count[b_id]) {
                            parent_id = a_id; child_id = b_id;
                        } else {
                            parent_id = b_id; child_id = a_id;
                        }
                        if (is_error_[child_id]) continue;
                        uint32_t eff_pid = parent_id;
                        if (is_error_[eff_pid]) eff_pid = find_root_chain(eff_pid);

                        // T5.8: posterior odds for adj-len indel edge. Treated
                        // as a single-mismatch edge at the indel position.
                        if (errcor_.empirical && !errcor_.legacy_veto) {
                            int L_child = static_cast<int>(arena_.length(child_id));
                            uint16_t pos16 = static_cast<uint16_t>(std::min(d, L_child - 1));
                            fqdup::errcor_emp::EdgeCandidate e{};
                            e.child_id   = child_id;
                            e.parent_id  = parent_id;
                            e.bundle_key = bundle_key_of[eff_pid];
                            e.bundle_occ = members.size();
                            e.child_len  = static_cast<uint16_t>(L_child);
                            e.n_mm       = 1;
                            e.mm[0].pos        = pos16;
                            e.mm[0].parent_2b  = 0;
                            e.mm[0].alt_2b     = 1;  // arbitrary non-equal pair: flagged as transversion
                            e.mm[0].qual       = qual_arena_.q_at(child_id, pos16);
                            e.mm[0].damage_chan = 0;
                            e.mm[0].p_damage    = 0.0;
                            double S = emp_model.score(e);
                            // Count-asymmetry term (see B2 H=2 site).
                            S += std::log(static_cast<double>(id_count[parent_id] + 1) /
                                          static_cast<double>(id_count[child_id]  + 1));
                            if (S <= 0.0) {
                                ++stats.adj_len_protected;
                                continue;
                            }
                        }

                        is_error_[child_id] = true;
                        parent_chain[child_id] = eff_pid;
                        acc_count[eff_pid] += acc_count[child_id];
                        acc_count[child_id] = 0;
                        if (seq_entry[child_id] && seq_entry[eff_pid] &&
                            seq_entry[child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                            seq_entry[eff_pid]->record_index = seq_entry[child_id]->record_index;
                            seq_entry[eff_pid]->damage_score = seq_entry[child_id]->damage_score;
                        }
                        ++stats.absorbed;
                        ++errcor_absorbed_;
                        ++stats.adj_len_absorbed;
                        // No fqcl record: ChildMismatch encodes substitutions only.
                        // The .fqcl genealogy will lack the indel edge but the
                        // dedup output is correct (is_error_ drops the read).
                    }
                }
            }
        }

        // ── T8.7: Phase3RescueIndels ──────────────────────────────────────────
        // Opt-in indel-rescue path. T8 rescue is strictly within-bundle
        // (the bundle filter is a hard correctness constraint, not a
        // scoring preference), so candidate generation iterates bundles
        // directly. For each bundle with ≥2 eligible members:
        //   - direct path (default): pairwise banded-ed over same-bundle
        //     parents with |Δlen|≤1
        //   - syncmer path (only above ~32k directed pairs): bundle-local
        //     ephemeral SyncmerIndex
        // Decisions go to per-thread vectors and are merged in (child_id,
        // score desc, parent_id asc) order, applying at most one
        // absorption per child.
        if (errcor_.rescue_indels && !errcor_.legacy_veto && errcor_.empirical) {
            using fqdup::derep_detail::SyncmerIndex;
            using fqdup::derep_detail::build_syncmer_index_from_pids;
            using fqdup::derep_detail::EditScript;
            using fqdup::derep_detail::banded_edit_distance_le2;
            using fqdup::errcor_emp::EdgeCandidate;
            using fqdup::errcor_emp::IndelScoreParams;

            IndelScoreParams ip;
            ip.alpha_ins  = errcor_.rescue_alpha_ins;
            ip.alpha_del  = errcor_.rescue_alpha_del;
            ip.mask_bonus = errcor_.rescue_mask_bonus;

            auto child_bundle_occ_of = rescue_bundle_occ_of;

            // T8 Step 1 — Pre-group eligible ids by rescue_bundle_key. Bundle
            // membership is the hard correctness gate; iterating bundles
            // directly makes the SyncmerIndex unnecessary for the common case
            // (mean bundle size ~3 → direct banded-ed beats build+query+sort
            // overhead). The pre-grouping pass also (a) accounts the bundle-
            // hot skip into stats.rescue_skip_bundle_hot (previously silent;
            // GPT-5.5 audit), and (b) computes pairs_est used to drive LPT
            // scheduling and the syncmer-path threshold.
            ska::flat_hash_map<uint64_t, std::vector<uint32_t>> bundle_members;
            bundle_members.reserve(N / 4 + 16);
            for (uint32_t id = 0; id < N; ++id) {
                if (!arena_.is_eligible(id)) continue;
                if (!ec_eligible(id)) continue;   // no short-interior read enters rescue-indel absorption
                if (is_error_[id]) continue;
                uint32_t occ = static_cast<uint32_t>(child_bundle_occ_of(id));
                if (occ < 2u) continue;
                if (occ >= errcor_.rescue_bundle_hot) {
                    // Fix silent recall loss: bundle-hot skip was uncounted.
                    ++stats.rescue_skip_bundle_hot;
                    continue;
                }
                bundle_members[rescue_bundle_key_of[id]].push_back(id);
            }

            struct BundleTask {
                uint64_t                key;
                std::vector<uint32_t>*  members;     // borrowed
                uint32_t                M;
                uint64_t                pairs_est;   // est. directed pairs |Δlen|≤2
            };
            std::vector<BundleTask> tasks;
            tasks.reserve(bundle_members.size());
            for (auto& kv : bundle_members) {
                if (kv.second.size() < 2) continue;
                std::sort(kv.second.begin(), kv.second.end());
                const auto& bm = kv.second;
                const uint32_t M = static_cast<uint32_t>(bm.size());
                // Length histogram. M is bounded by rescue_bundle_hot (≤50
                // by default) so a flat hashmap is overkill; sort the lens
                // and run-length-encode.
                std::vector<int> lens; lens.reserve(M);
                for (uint32_t i = 0; i < M; ++i) lens.push_back(arena_.length(bm[i]));
                std::sort(lens.begin(), lens.end());
                uint64_t pairs_est = 0;
                size_t i = 0;
                while (i < lens.size()) {
                    int Lv = lens[i];
                    size_t j = i + 1;
                    while (j < lens.size() && lens[j] == Lv) ++j;
                    uint64_t cnt = j - i;
                    // Window neighbours within ±2 (inclusive) of Lv.
                    uint64_t neigh = 0;
                    for (size_t k = 0; k < lens.size(); ++k)
                        if (std::abs(lens[k] - Lv) <= 2) ++neigh;
                    // Exclude self from neighbour count.
                    pairs_est += cnt * (neigh - 1);
                    i = j;
                }
                tasks.push_back(BundleTask{kv.first, &kv.second, M, pairs_est});
            }
            // LPT scheduling: biggest-pair-count bundles start first so tail
            // latency is minimised. Tie-break by key for determinism.
            std::sort(tasks.begin(), tasks.end(),
                      [](const BundleTask& a, const BundleTask& b) {
                          if (a.pairs_est != b.pairs_est) return a.pairs_est > b.pairs_est;
                          return a.key < b.key;
                      });

            struct RescueDecision {
                uint32_t child_id;
                uint32_t parent_id;
                double   score;
                uint8_t  ed;
            };

            unsigned n_threads = errcor_.threads;
            if (n_threads == 0) n_threads = std::max(1u, std::thread::hardware_concurrency());

            struct alignas(64) T8Local {
                std::vector<RescueDecision> dec;
                uint64_t children_examined = 0;
                uint64_t index_queries     = 0;
                uint64_t topk_truncated    = 0;
                uint64_t pairs_banded      = 0;
                uint64_t banded_reject     = 0;
                uint64_t protected_        = 0;
                uint64_t banded_ed[3]      = {0, 0, 0};
                char _pad[8] = {};
            };
            std::vector<T8Local> per_thread(n_threads);
            std::vector<SyncmerIndex::QueryScratch> per_thread_qs(n_threads);

            // Above this directed-pair count per bundle, switch to a bundle-
            // local SyncmerIndex. Below, direct banded-ed is cheaper than
            // build+query overhead. Long-tail bundles trigger the syncmer
            // path; mean-bundle case never does.
            constexpr uint64_t kSyncmerPairThreshold = 32000;

            // T8 Step 2 — atomic dynamic scheduler (replaces static round-robin).
            // Combined with LPT ordering this keeps tail-thread idle time small
            // even when bundle work is heavily skewed.
            std::atomic<size_t> next_bundle{0};

            auto worker = [&](unsigned tid) {
                T8Local& T = per_thread[tid];
                SyncmerIndex::QueryScratch& qs = per_thread_qs[tid];
                // T8 Step 3 / Step 7 — per-thread reusable scratch.
                thread_local std::vector<uint8_t>  decoded;
                thread_local std::vector<uint32_t> off;          // size = M+1
                thread_local std::vector<int>      L;            // per-member length
                thread_local std::vector<uint32_t> order_by_len; // local indices sorted by L
                thread_local std::vector<uint64_t> ch_h;
                thread_local std::vector<uint16_t> ch_p;
                thread_local ska::flat_hash_map<uint32_t, uint32_t> pid_to_local;
                // Step 7 — reusable bundle-local SyncmerIndex (zero allocs after warmup).
                thread_local SyncmerIndex bidx;
                thread_local std::vector<std::pair<uint64_t, uint32_t>> pairs_scratch;

                while (true) {
                    size_t bi = next_bundle.fetch_add(1, std::memory_order_relaxed);
                    if (bi >= tasks.size()) break;
                    const BundleTask& task = tasks[bi];
                    const auto& bm = *task.members;
                    const size_t M = task.M;
                    if (M < 2) continue;

                    // Step 3 — cache lengths and compute decode offsets in one pass.
                    L.assign(M, 0);
                    off.assign(M + 1, 0);
                    for (size_t i = 0; i < M; ++i) {
                        L[i] = arena_.length(bm[i]);
                        off[i + 1] = off[i] + L[i];
                    }
                    if (decoded.size() < off.back()) decoded.resize(off.back());
                    for (size_t i = 0; i < M; ++i)
                        arena_.decode_range(bm[i], 0, L[i], decoded.data() + off[i]);

                    // Local indices sorted by (length asc, parent_id asc) — used
                    // for the 2-pointer / binary-search Δlen window in the
                    // direct path. Determinism via the parent_id tiebreak.
                    order_by_len.resize(M);
                    for (uint32_t i = 0; i < M; ++i) order_by_len[i] = i;
                    std::sort(order_by_len.begin(), order_by_len.end(),
                              [&](uint32_t a, uint32_t b) {
                                  if (L[a] != L[b]) return L[a] < L[b];
                                  return bm[a] < bm[b];
                              });

                    // T8 Step 1 / Step 3 — pairs_est drives the path switch.
                    // Old code recomputed an O(M²) probe here; deleted entirely.
                    const bool over_threshold = (task.pairs_est > kSyncmerPairThreshold);

                    if (over_threshold) {
                        // Step 7 — build into reusable thread_local index.
                        build_syncmer_index_from_decoded(
                            bidx,
                            decoded.data(), off.data(), L.data(),
                            bm.data(), M,
                            pairs_scratch,
                            errcor_.rescue_hash_hot);
                        pid_to_local.clear();
                        if (pid_to_local.bucket_count() < M * 2)
                            pid_to_local.reserve(M * 2);
                        for (size_t i = 0; i < M; ++i)
                            pid_to_local[bm[i]] = static_cast<uint32_t>(i);
                    }

                    // T8 Step 5 — refactored try_pair. Returns through the
                    // best_* refs; the per-child best is pushed to T.dec ONCE
                    // after the parent loop (Step 4 / GPT-5.5 P1.6).
                    auto try_pair = [&](uint32_t child_id, uint32_t parent_id,
                                        int Lc, int Lp,
                                        const uint8_t* child_dec_p,
                                        const uint8_t* parent_dec_p,
                                        int k5_c, int k3_c,
                                        double& best_score, uint32_t& best_parent,
                                        uint8_t& best_ed) {
                        // |Δlen|≤2 already filtered by callers; defensive.
                        if (std::abs(Lp - Lc) > 2) return;

                        EditScript es{};
                        int ed = banded_edit_distance_le2(
                            parent_dec_p, Lp, child_dec_p, Lc, &es);
                        ++T.pairs_banded;
                        if (ed < 0) { ++T.banded_reject; return; }
                        if (ed >= 0 && ed < 3) ++T.banded_ed[ed];

                        EdgeCandidate e{};
                        e.child_id   = child_id;
                        e.parent_id  = parent_id;
                        e.bundle_key = bundle_key_of[parent_id];
                        e.bundle_occ = bundle_occ_of(parent_id);
                        e.child_len  = static_cast<uint16_t>(Lc);
                        e.n_mm       = 0;
                        e.n_ins      = static_cast<uint8_t>(es.n_ins);
                        e.n_del      = static_cast<uint8_t>(es.n_del);
                        e.indel_in_mask = 0;
                        for (int s = 0; s < 2; ++s) {
                            if (es.kind[s] != 1 && es.kind[s] != 2) continue;
                            int p = es.pos[s];
                            if (p < k5_c || p >= Lc - k3_c) {
                                e.indel_in_mask = 1; break;
                            }
                        }
                        for (int s = 0; s < 2 && e.n_mm < 2; ++s) {
                            if (es.kind[s] != 0) continue;
                            int p = es.pos[s];
                            auto& mm = e.mm[e.n_mm];
                            mm.pos        = static_cast<uint16_t>(p);
                            mm.parent_2b  = es.ref_base[s];
                            mm.alt_2b     = es.alt_base[s];
                            mm.qual       = qual_arena_.q_at(child_id, p);
                            mm.damage_chan = 0;
                            mm.p_damage    = 0.0;
                            ++e.n_mm;
                        }

                        double S = emp_model.score_with_indel(e, ip);
                        if (!(S > 0.0)) { ++T.protected_; return; }
                        // Track per-child best (score desc, parent_id asc).
                        if (S > best_score ||
                            (S == best_score && parent_id < best_parent)) {
                            best_score  = S;
                            best_parent = parent_id;
                            best_ed     = static_cast<uint8_t>(ed);
                        }
                    };

                    // id_count gate hoisted out of try_pair (cheap; pre-prunes
                    // banded-ed work).
                    auto id_count_admits = [&](uint32_t child_id, uint32_t parent_id) {
                        if (parent_id == child_id) return false;
                        if (id_count[parent_id] < id_count[child_id]) return false;
                        if (id_count[parent_id] == id_count[child_id] &&
                            parent_id > child_id) return false;
                        return true;
                    };

                    for (size_t ci_local = 0; ci_local < M; ++ci_local) {
                        const uint32_t child_id = bm[ci_local];
                        ++T.children_examined;
                        const int Lc = L[ci_local];
                        const uint8_t* child_dec_p = decoded.data() + off[ci_local];
                        auto [k5_c, k3_c] = damage_zone_bounds(Lc, profile_);

                        double   best_score  = 0.0;
                        uint32_t best_parent = UINT32_MAX;
                        uint8_t  best_ed     = 0;

                        if (!over_threshold) {
                            // T8 Step 4 — direct path. |Δlen|≤2 (was ≤1; the
                            // banded_ed cap supports 2 — silent recall miss
                            // flagged by GPT-5.5). Iterate via order_by_len
                            // window; M ≤ rescue_bundle_hot (≤50 default), so
                            // std::lower_bound is fine.
                            auto lo = std::lower_bound(
                                order_by_len.begin(), order_by_len.end(), Lc - 2,
                                [&](uint32_t a, int v) { return L[a] < v; });
                            auto hi = std::upper_bound(
                                order_by_len.begin(), order_by_len.end(), Lc + 2,
                                [&](int v, uint32_t a) { return v < L[a]; });
                            for (auto it = lo; it != hi; ++it) {
                                const uint32_t pi_local = *it;
                                if (pi_local == ci_local) continue;
                                const uint32_t parent_id = bm[pi_local];
                                if (!id_count_admits(child_id, parent_id)) continue;
                                try_pair(child_id, parent_id,
                                         Lc, L[pi_local],
                                         child_dec_p,
                                         decoded.data() + off[pi_local],
                                         k5_c, k3_c,
                                         best_score, best_parent, best_ed);
                            }
                        } else {
                            fqdup::SyncmerParams sp = fqdup::syncmer_params_for_length(Lc);
                            uint16_t min_hits = errcor_.rescue_min_hits > 0
                                ? static_cast<uint16_t>(errcor_.rescue_min_hits)
                                : (sp.cap >= 16 ? uint16_t{3} : uint16_t{2});
                            fqdup::compute_sketch(child_dec_p, Lc, ch_h, ch_p);
                            if (ch_h.empty()) continue;
                            ++T.index_queries;
                            // T8 Step 6 — filtered query. Pre-tally rejection
                            // by length AND id_count avoids same-bundle
                            // wrong-length parents stealing topk slots from
                            // valid same-length ones (silent recall loss).
                            auto accept = [&](uint32_t pid) -> bool {
                                auto pit = pid_to_local.find(pid);
                                if (pit == pid_to_local.end()) return false;
                                int Lp = L[pit->second];
                                if (std::abs(Lp - Lc) > 2) return false;
                                return id_count_admits(child_id, pid);
                            };
                            const auto& hits = bidx.query_filtered(
                                qs, ch_h.data(), (int)ch_h.size(),
                                min_hits, errcor_.rescue_topk, accept);
                            if (hits.size() == errcor_.rescue_topk) ++T.topk_truncated;
                            for (auto& h : hits) {
                                auto pit = pid_to_local.find(h.parent_id);
                                if (pit == pid_to_local.end()) continue;
                                const uint32_t pi_local = pit->second;
                                try_pair(child_id, h.parent_id,
                                         Lc, L[pi_local],
                                         child_dec_p,
                                         decoded.data() + off[pi_local],
                                         k5_c, k3_c,
                                         best_score, best_parent, best_ed);
                            }
                        }

                        // Push at most ONE decision per child — the global
                        // sort + first-wins below reduces to identity.
                        if (best_parent != UINT32_MAX) {
                            T.dec.push_back(RescueDecision{
                                child_id, best_parent, best_score, best_ed});
                        }
                    }
                }
            };

            if (n_threads == 1) {
                worker(0);
            } else {
                std::vector<std::thread> ts;
                ts.reserve(n_threads - 1);
                for (unsigned t = 1; t < n_threads; ++t) ts.emplace_back(worker, t);
                worker(0);
                for (auto& t : ts) t.join();
            }

            // Fold per-thread sub-counters and decisions into globals.
            std::vector<RescueDecision> decisions;
            for (auto& T : per_thread) {
                stats.rescue_children_examined += T.children_examined;
                stats.rescue_index_queries     += T.index_queries;
                stats.rescue_topk_truncated    += T.topk_truncated;
                stats.rescue_pairs_banded      += T.pairs_banded;
                stats.rescue_banded_reject     += T.banded_reject;
                stats.rescue_protected         += T.protected_;
                for (int e = 0; e < 3; ++e)
                    stats.rescue_banded_ed[e]  += T.banded_ed[e];
                decisions.insert(decisions.end(),
                                 T.dec.begin(), T.dec.end());
            }

            // Deterministic merge: at most one absorption per child, picking
            // the highest-scoring parent (ties broken by parent_id asc).
            std::sort(decisions.begin(), decisions.end(),
                      [](const RescueDecision& a, const RescueDecision& b) {
                          if (a.child_id != b.child_id) return a.child_id < b.child_id;
                          if (a.score    != b.score)    return a.score    > b.score;
                          return a.parent_id < b.parent_id;
                      });
            uint32_t prev_child = UINT32_MAX;
            for (auto& d : decisions) {
                if (d.child_id == prev_child) continue;
                prev_child = d.child_id;
                if (is_error_[d.child_id]) continue;
                uint32_t eff_pid = d.parent_id;
                if (is_error_[eff_pid]) eff_pid = find_root_chain(eff_pid);
                if (eff_pid == d.child_id) continue;
                is_error_[d.child_id] = true;
                parent_chain[d.child_id] = eff_pid;
                acc_count[eff_pid]      += acc_count[d.child_id];
                acc_count[d.child_id]    = 0;
                if (seq_entry[d.child_id] && seq_entry[eff_pid] &&
                    seq_entry[d.child_id]->damage_score > seq_entry[eff_pid]->damage_score) {
                    seq_entry[eff_pid]->record_index = seq_entry[d.child_id]->record_index;
                    seq_entry[eff_pid]->damage_score = seq_entry[d.child_id]->damage_score;
                }
                ++stats.absorbed;
                ++stats.rescue_absorbed;
                int ob = static_cast<int>(fqdup::errcor_emp::occ_bin(rescue_bundle_occ_of(eff_pid)));
                if (ob >= 0 && ob < 6) ++stats.rescue_absorbed_by_occ[ob];
                ++errcor_absorbed_;
                // No fqcl record (indel edge — same constraint as adj-len).
            }
        }

        // Write accumulated cluster counts back into the index so pass2() and the
        // cluster TSV reflect absorbed PCR error reads, not just the original
        // representative count.  acc_count is a local — must be written back before
        // phase3_error_correct() returns.
        index_.for_each([&](const SequenceFingerprint&, IndexEntry& entry) {
            if (entry.seq_id < static_cast<uint32_t>(acc_count.size()) &&
                !is_error_[entry.seq_id])
                entry.count = acc_count[entry.seq_id];
        });

        stats.log();
        // Persist for write_fqcl_ → metadata.
        loss_bucket_overflow_drops_   = stats.bucket_overflow_drops;
        loss_short_brute_evaluated_   = stats.short_brute_evaluated;   // 0: brute-force fallback removed
        loss_short_brute_found_       = stats.short_brute_found;       // 0
        loss_short_too_small_skipped_ = stats.short_too_small_skipped; // 0
        // Short interiors are no longer a loss: they pass through as their own
        // Phase-1-exact-deduped representatives. This is a kept-count, not a drop.
        loss_short_interior_skipped_ = short_interior_skipped;
        if (short_interior_skipped > 0) {
            double pct_short = 100.0 * static_cast<double>(short_interior_skipped) /
                               static_cast<double>(N);
            log_info("Phase 3: " + std::to_string(short_interior_skipped) +
                     " reads have interior < " + std::to_string(kMinInteriorLen) +
                     " bp (" + std::to_string(pct_short).substr(0, 4) + "%); "
                     "kept as own representatives (Phase-1 exact-deduped, EC skipped)");
            if (pct_short > 10.0)
                log_warn("Over 10% of reads have interiors too short for the pigeonhole index "
                         "after damage masking. Consider reducing --mask-threshold or "
                         "adjusting damage zone parameters.");
        }
        log_info("Phase 3 complete: absorbed " + std::to_string(errcor_absorbed_) +
                 " sequences (H=1 directed-ascent + H=2 direct-detection, snp_threshold=" +
                 std::to_string(errcor_.snp_threshold) + ")");
        log_info("  Joint adj fired: neg=" + std::to_string(n_joint_adj_neg) +
                 " pos=" + std::to_string(n_joint_adj_pos) +
                 " (of " + std::to_string(n_parents_processed) + " parents)");

        if (!fqcl_path_.empty()) {
            fqcl_parent_chain_ = std::move(parent_chain);
            write_fqcl_(seq_entry);
        }
