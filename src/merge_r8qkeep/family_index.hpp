#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace fqdup::r8qkeep {
// Compact postings over distinct 13/14-base parts, not over all reads.
class FamilyPartIndex {
    struct Part { std::uint64_t value; int begin, end; };
    using Entry = std::pair<std::uint32_t, int>;
    std::vector<int> members_;
    std::vector<Part> parts_;
    std::array<std::vector<Entry>, 3> blocks_;
    std::array<int, 4> bounds_{};
    int length_ = 0;
    static int code(char c) {
        switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return 4; }
    }
    std::uint32_t block(std::uint64_t value, int b) const {
        return (value >> (3 * bounds_[b])) & ((std::uint64_t(1) << (3 * (bounds_[b+1] - bounds_[b]))) - 1);
    }
    bool has_n(std::uint64_t value, int b) const {
        for (int p = bounds_[b]; p < bounds_[b+1]; ++p) if (((value >> (3*p)) & 7) == 4) return true;
        return false;
    }
public:
    static std::uint64_t pack(const std::string& key, int begin, int end) {
        std::uint64_t value = 0;
        for (int p = begin; p < end; ++p) value |= std::uint64_t(code(key[p])) << (3 * (p-begin));
        return value;
    }
    void build(const std::vector<std::string>& keys, const std::vector<int>& representatives, int begin, int end) {
        length_ = end - begin;
        bounds_ = {0, 5, length_ == 13 ? 9 : 10, length_};
        members_.resize(representatives.size());
        std::iota(members_.begin(), members_.end(), 0);
        const auto value = [&](int i) { return pack(keys[representatives[i]], begin, end); };
        std::sort(members_.begin(), members_.end(), [&](int a, int b) {
            const auto x = value(a), y = value(b); return x < y || (x == y && a < b);
        });
        std::size_t count = 0;
        for (std::size_t i = 0; i < members_.size(); ++i) if (!i || value(members_[i]) != value(members_[i-1])) ++count;
        parts_.reserve(count);
        for (std::size_t i = 0; i < members_.size();) {
            const auto v = value(members_[i]); const int first = i++;
            while (i < members_.size() && value(members_[i]) == v) ++i;
            parts_.push_back({v, first, int(i)});
        }
        for (int b = 0; b < 3; ++b) {
            blocks_[b].reserve(parts_.size());
            for (std::size_t i = 0; i < parts_.size(); ++i) blocks_[b].emplace_back(block(parts_[i].value, b), i);
            std::sort(blocks_[b].begin(), blocks_[b].end());
        }
    }
    // A query N consumes one mismatch even against stored N. A <=2-distance
    // part therefore has a completely equal, N-free block among its three.
    // The visitor returns false to stop without allocating a candidate list.
    template<class Visitor>
    void visit(std::uint64_t query, Visitor visitor, std::uint64_t& lookups, std::uint64_t& candidates) const {
        for (int b = 0; b < 3; ++b) {
            if (has_n(query, b)) continue;
            const auto target = block(query, b);
            ++lookups;
            const auto& entries = blocks_[b];
            auto it = std::lower_bound(entries.begin(), entries.end(), Entry{target, 0});
            for (; it != entries.end() && it->first == target; ++it) {
                const auto& part = parts_[it->second];
                bool earlier = false;
                for (int a = 0; a < b; ++a) if (!has_n(query, a) && block(query, a) == block(part.value, a)) { earlier = true; break; }
                if (earlier) continue;
                ++candidates;
                int differences = 0;
                for (int p = 0; p < length_ && differences <= 2; ++p) {
                    const int q = (query >> (3*p)) & 7, s = (part.value >> (3*p)) & 7;
                    differences += q == 4 || q != s;
                }
                if (differences > 2) continue;
                for (int i = part.begin; i < part.end; ++i) if (!visitor(members_[i])) return;
            }
        }
    }
    std::size_t cells() const { return parts_.size(); }
    std::size_t max_bucket() const {
        std::size_t n = 0; for (const auto& p : parts_) n = std::max(n, std::size_t(p.end-p.begin)); return n;
    }
    std::size_t bytes() const {
        std::size_t n = members_.capacity()*sizeof(int) + parts_.capacity()*sizeof(Part);
        for (const auto& b : blocks_) n += b.capacity()*sizeof(Entry);
        return n;
    }
};
}
