#pragma once
// Shard-of-hash-maps wrapper for the derep dedup index.
//
// The index is the largest single structure in a derep run: 32 GB of slots at
// 355 M uniques. As one ska::flat_hash_map, crossing the load ceiling rehashes
// the whole table at once — the new slot array is allocated while the old one is
// still live, so a 2^29 -> 2^30 growth needs 32 + 64 = 96 GB resident at that
// instant, on top of the arenas. Pre-sizing does not remove the risk, only moves
// it: the estimate comes from the input file size, but the number of *unique*
// reads depends on the duplication rate, which is unknowable before Pass 1 runs.
// A library whose uniques overshoot the estimate pays the full spike.
//
// Sharding bounds it. Each shard grows independently, so the transient is one
// shard's doubling (~0.5 GB at 64 shards) rather than the whole table's. Memory
// stops being a cliff and becomes a slope.
//
// Shard selection uses the HIGH bits of the hash. ska's fibonacci_hash_policy
// derives its slot index from the high bits of (hash * 2^64/phi), a product that
// depends on every bit of the input, so the two are not correlated and the shards
// do not clump.
//
// Byte-identity: iteration order across shards differs from a single map's, and
// that is safe here because it never reaches the output. The one consumer that
// feeds the emitted FASTQ (derep.cpp pass2) collects entries and immediately
// sorts them by record_index; the other two only scatter into arrays indexed by
// seq_id. Gated on md5 of the emitted FASTQ + clusters regardless.

#include "flat_hash_map.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace fqdup::derep_detail {

template <typename K, typename V, typename H, unsigned ShardBits = 6>
class ShardedIndex {
  public:
    using Map      = ska::flat_hash_map<K, V, H>;
    using iterator = typename Map::iterator;

    static constexpr size_t kShardBits = ShardBits;
    static constexpr size_t kShards    = size_t{1} << ShardBits;

    void max_load_factor(float f) {
        for (auto& m : shards_) m.max_load_factor(f);
    }

    // total_slots is the whole-index slot budget; split it evenly. Undershooting
    // is now cheap — an overflowing shard regrows alone.
    void rehash(size_t total_slots) {
        const size_t per = (total_slots + kShards - 1) / kShards;
        for (auto& m : shards_) m.rehash(per);
    }

    std::pair<iterator, bool> emplace(const K& key, const V& value) {
        return shards_[shard_of(key)].emplace(key, value);
    }

    size_t size() const {
        size_t n = 0;
        for (const auto& m : shards_) n += m.size();
        return n;
    }

    size_t bucket_count() const {
        size_t n = 0;
        for (const auto& m : shards_) n += m.bucket_count();
        return n;
    }

    template <typename F>
    void for_each(F&& f) {
        for (auto& m : shards_)
            for (auto& kv : m) f(kv.first, kv.second);
    }

    template <typename F>
    void for_each(F&& f) const {
        for (const auto& m : shards_)
            for (const auto& kv : m) f(kv.first, kv.second);
    }

    // Release the slot arrays. ska has no shrinking clear(), so swap each shard
    // against a fresh map.
    void release() {
        for (auto& m : shards_) {
            Map empty;
            std::swap(m, empty);
        }
    }

  private:
    size_t shard_of(const K& key) const {
        return static_cast<size_t>(static_cast<uint64_t>(H{}(key)) >>
                                   (64 - ShardBits));
    }

    std::array<Map, kShards> shards_;
};

}  // namespace fqdup::derep_detail
