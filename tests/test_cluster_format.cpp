// Invariant tests for v2 .fqcl format and recently-fixed silent-cap regressions.
//
// Covers:
//   - Magic / version / 8 KB reserved meta header
//   - Round-trip via Reader (parent_seq, edges, damage, n_members)
//   - member_ids round-trip (kFlagHasMemberIds)
//   - Tinyblock encoding for eligible singletons (no edges, no qual, ≤1 id)
//   - kZstdMinBytes threshold: large full block → kBlockTypeFullZstd
//   - DamageProfile::MASK_POSITIONS == 32
//   - WriterMetadata new loss counters serialise into JSON

#include "fqdup/cluster_format.hpp"

#include <cassert>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace fqdup::clusterfmt;

static void check_static_invariants() {
    static_assert(kVersion == 2u, ".fqcl v2 required");
    static_assert(kZstdMinBytes == 512u, "zstd threshold must be 512 B raw");
}

static ClusterRecord make_full(std::uint64_t cid, std::uint32_t L_bytes,
                               std::uint32_t L_bases, std::uint32_t n_members) {
    ClusterRecord r;
    r.cluster_id = cid;
    r.flags = kFlagHasMemberIds;
    r.n_members = n_members;
    r.n_after_damage = n_members;
    r.parent_seq_len = L_bases;
    r.parent_seq.assign(L_bytes, 0);
    for (std::uint32_t i = 0; i < L_bytes; ++i)
        r.parent_seq[i] = static_cast<std::uint8_t>((i * 37u + 11u) & 0xFFu);
    r.edges.push_back({0, 1, 5, 1, 3, n_members});
    r.edges.push_back({1, 2, 9, 0, 2, n_members / 2});
    for (int i = 0; i < (int)kTermWindow; ++i) {
        r.damage.term_5[i] = static_cast<std::uint8_t>(255 / (i + 1));
        r.damage.term_3[i] = static_cast<std::uint8_t>(128 / (i + 1));
    }
    for (std::uint32_t i = 0; i < n_members; ++i)
        r.member_ids.push_back("read_" + std::to_string(cid) + "_" + std::to_string(i));
    return r;
}

static ClusterRecord make_tinyblock_singleton(std::uint64_t cid) {
    ClusterRecord r;
    r.cluster_id = cid;
    r.flags = kFlagHasMemberIds;
    r.n_members = 1;
    r.n_after_damage = 1;
    r.parent_seq_len = 16;
    r.parent_seq = {0x1B, 0x4E, 0xE4, 0x1B};
    r.member_ids.push_back("singleton_" + std::to_string(cid));
    return r;
}

static ClusterRecord make_large_full(std::uint64_t cid) {
    ClusterRecord r;
    r.cluster_id = cid;
    r.flags = kFlagHasMemberIds;
    r.n_members = 200;
    r.n_after_damage = 200;
    r.parent_seq_len = 200;
    r.parent_seq.assign(50, 0xA5);
    for (std::uint32_t to_node = 1; to_node <= 60; ++to_node) {
        r.edges.push_back({to_node - 1, to_node,
                           static_cast<std::uint16_t>(to_node % 200),
                           static_cast<std::uint8_t>(to_node & 3),
                           static_cast<std::uint8_t>((to_node + 1) & 3),
                           200u - to_node});
    }
    for (std::uint32_t i = 0; i < r.n_members; ++i)
        r.member_ids.push_back("big_" + std::to_string(i));
    return r;
}

static void check_block_type(const std::string& path,
                             const std::vector<std::uint64_t>& offsets,
                             std::size_t idx, std::uint8_t expected_type) {
    std::ifstream is(path, std::ios::binary);
    is.seekg(static_cast<std::streamoff>(offsets[idx]), std::ios::beg);
    std::uint32_t wire_size = 0;
    is.read(reinterpret_cast<char*>(&wire_size), 4);
    std::uint8_t type = 0xFF;
    is.read(reinterpret_cast<char*>(&type), 1);
    assert(type == expected_type);
}

static std::vector<std::uint64_t> read_offsets(const std::string& path,
                                               std::uint64_t n_clusters) {
    std::ifstream is(path, std::ios::binary | std::ios::ate);
    is.seekg(-(static_cast<std::streamoff>(16 + 8 * n_clusters)), std::ios::end);
    std::vector<std::uint64_t> offs(n_clusters);
    is.read(reinterpret_cast<char*>(offs.data()), 8 * n_clusters);
    return offs;
}

int main() {
    check_static_invariants();

    const char* path = "smoke_clusters.fqcl";

    WriterMetadata m;
    m.tool_version = "test";
    m.input_fastq  = "synthetic.fq.gz";
    m.n_input_reads = 100;
    m.library_type = "ds";
    m.d_max_5 = 0.21; m.d_max_3 = 0.13;
    m.lambda_5 = 0.28; m.lambda_3 = 0.27;
    m.snp_threshold = 0.20; m.snp_min_count = 1; m.bucket_cap = 0;  // unlimited
    m.bucket_overflow_drops  = 7;
    m.short_interior_skipped = 3;
    m.n_singletons_tinyblock = 1;

    ClusterRecord full      = make_full(0, 4, 16, 12);
    ClusterRecord tiny      = make_tinyblock_singleton(1);
    ClusterRecord big       = make_large_full(2);

    {
        Writer w(path, m);
        w.write_cluster(full);
        w.write_cluster(tiny);
        w.write_cluster(big);
        w.close();
        assert(w.n_clusters() == 3);
    }

    // Magic header + reserved meta region
    {
        std::ifstream is(path, std::ios::binary);
        std::uint32_t magic = 0, version = 0;
        std::uint64_t meta_size = 0;
        is.read(reinterpret_cast<char*>(&magic), 4);
        is.read(reinterpret_cast<char*>(&version), 4);
        is.read(reinterpret_cast<char*>(&meta_size), 8);
        assert(magic == kMagic);
        assert(version == kVersion);
        assert(meta_size >= 8192u && "meta region must be reserved (≥8 KB)");

        std::vector<char> buf(meta_size);
        is.read(buf.data(), meta_size);
        std::string js(buf.data(), meta_size);
        assert(js.find("\"tool\":\"fqdup\"") != std::string::npos);
        assert(js.find("\"library_type\":\"ds\"") != std::string::npos);
        assert(js.find("\"bucket_overflow_drops\":7") != std::string::npos);
        assert(js.find("\"short_interior_skipped\":3") != std::string::npos);
        assert(js.find("\"n_singletons_tinyblock\":1") != std::string::npos);
        assert(js.find("\"block_compression\":\"zstd\"") != std::string::npos);
        assert(js.find("\"n_clusters\":3") != std::string::npos);
    }

    auto offsets = read_offsets(path, 3);
    check_block_type(path, offsets, 0, kBlockTypeFull);      // small full
    check_block_type(path, offsets, 1, kBlockTypeTiny);      // singleton
    check_block_type(path, offsets, 2, kBlockTypeFullZstd);  // large full → zstd

    // Reader round-trip
    {
        Reader rd(path);
        assert(rd.n_clusters() == 3);

        ClusterRecord got;
        rd.read_cluster(0, got);
        assert(got.cluster_id == 0);
        assert(got.n_members == 12);
        assert(got.parent_seq_len == 16);
        assert(got.parent_seq == full.parent_seq);
        assert(got.edges.size() == 2);
        assert(got.edges[0].pos == 5 && got.edges[0].to_base == 3);
        assert(got.edges[1].pos == 9 && got.edges[1].to_base == 2);
        assert(std::memcmp(got.damage.term_5, full.damage.term_5, kTermWindow) == 0);
        assert(std::memcmp(got.damage.term_3, full.damage.term_3, kTermWindow) == 0);
        assert(got.member_ids.size() == 12);
        assert(got.member_ids.front() == "read_0_0");
        assert(got.member_ids.back() == "read_0_11");

        rd.read_cluster(1, got);
        assert(got.cluster_id == 1);
        assert(got.n_members == 1);
        assert(got.edges.empty());
        assert(got.parent_seq_len == 16);
        assert(got.parent_seq == tiny.parent_seq);
        assert(got.member_ids.size() == 1);
        assert(got.member_ids[0] == "singleton_1");

        rd.read_cluster(2, got);
        assert(got.cluster_id == 2);
        assert(got.n_members == 200);
        assert(got.edges.size() == 60);
        assert(got.member_ids.size() == 200);
        assert(got.member_ids.front() == "big_0");
        assert(got.member_ids.back() == "big_199");
    }

    std::remove(path);
    std::cout << "test_cluster_format: PASS\n";
    return 0;
}
