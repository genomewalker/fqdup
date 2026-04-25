#pragma once

// fqdup cluster genealogy format (.fqcl) — writer/reader API.
// Spec: wiki/cluster-format.md
//
// Persists per-cluster: parent sequence, member count, edge tree of single-base
// edits, terminal damage profile. Block-oriented, footer-indexed for random
// access. v1 is uncompressed; zstd in v2.

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace fqdup::clusterfmt {

inline constexpr std::uint32_t kMagic       = 0x4C434651u; // 'FQCL' little-endian
inline constexpr std::uint32_t kVersion     = 2u;
inline constexpr std::size_t   kTermWindow  = 8u;          // damage_term_{5,3}

// Per-block wire format (v2):
//   uint32_t wire_size  — bytes that follow (header + payload)
//   uint8_t  block_type — see kBlockType*
//   payload (varies by type)
enum BlockType : std::uint8_t {
    kBlockTypeFull       = 0,  // full uncompressed serialise_block payload
    kBlockTypeTiny       = 1,  // singleton compact payload (no edges, no qual)
    kBlockTypeFullZstd   = 2,  // full payload, zstd-compressed; preceded by uint32 uncompressed_size
};

// Threshold above which full blocks get zstd-compressed (<=512 B raw → leave uncompressed).
inline constexpr std::size_t kZstdMinBytes = 512u;

// One edit linking a child node to its parent in the cluster's tree.
struct Edge {
    std::uint32_t from_node;   // 0 = parent_seq root
    std::uint32_t to_node;
    std::uint16_t pos;         // position in parent_seq
    std::uint8_t  from_base;   // 2-bit, low bits
    std::uint8_t  to_base;     // 2-bit
    std::uint32_t n_reads;     // size of subtree at to_node
};
static_assert(sizeof(Edge) == 16, "Edge layout drift");

// Per-cluster damage histogram, fixed-window terminal C->T (5') / G->A (3').
struct DamageTrack {
    std::uint8_t term_5[kTermWindow] = {};   // C->T fraction × 255
    std::uint8_t term_3[kTermWindow] = {};   // G->A fraction × 255
};

// In-memory representation of one cluster, fed to Writer::write_cluster.
struct ClusterRecord {
    std::uint64_t                cluster_id     = 0;
    std::uint32_t                flags          = 0;
    std::uint32_t                n_members      = 0;
    std::uint32_t                n_after_damage = 0;
    std::vector<std::uint8_t>    parent_seq;       // 2-bit packed
    std::uint32_t                parent_seq_len = 0;
    std::vector<std::uint8_t>    parent_qual;      // raw Phred, may be empty
    std::vector<Edge>            edges;
    DamageTrack                  damage;
    std::vector<std::string>     member_ids;       // optional; gated by flag
};

enum ClusterFlag : std::uint32_t {
    kFlagHasMemberIds  = 1u << 0,
    kFlagHasQuality    = 1u << 1,
    kFlagRevcompUsed   = 1u << 2,
};

// Run-once metadata serialised as JSON in the file header.
struct WriterMetadata {
    std::string   tool          = "fqdup";
    std::string   tool_version  = "0.0.0";
    std::string   input_fastq;
    std::uint64_t n_input_reads = 0;
    std::uint64_t n_clusters    = 0;     // filled by Writer::close()
    std::string   library_type  = "auto";
    double        d_max_5       = 0.0;
    double        d_max_3       = 0.0;
    double        lambda_5      = 0.0;
    double        lambda_3      = 0.0;
    int           mask_pos_5    = 0;
    int           mask_pos_3    = 0;
    double        snp_threshold = 0.0;
    int           snp_min_count = 0;
    int           bucket_cap    = 0;
    // v2: singletons are always included (tinyblock encoded). Loss counters surfaced here.
    std::uint64_t bucket_overflow_drops    = 0;
    std::uint64_t short_interior_skipped   = 0;
    std::uint64_t short_brute_evaluated    = 0;
    std::uint64_t short_brute_found        = 0;
    std::uint64_t short_too_small_skipped  = 0;
    std::uint64_t n_singletons_tinyblock   = 0;
    std::string   block_compression        = "zstd";  // "none" | "zstd"
};

// Streaming writer. Open(), write_cluster() once per cluster, close().
class Writer {
public:
    explicit Writer(const std::string& path, WriterMetadata meta);
    ~Writer();

    Writer(const Writer&) = delete;
    Writer& operator=(const Writer&) = delete;

    // Append a cluster. cluster_id need not be sequential but must be unique.
    void write_cluster(const ClusterRecord& rec);

    // Finalise: writes footer index + closes file. Idempotent.
    void close();

    std::uint64_t n_clusters() const noexcept { return n_clusters_; }

private:
    void write_magic_header_();
    void write_meta_header_();

    std::ofstream                ofs_;
    std::string                  path_;
    WriterMetadata               meta_;
    std::vector<std::uint64_t>   offsets_;     // file_offset per cluster, in write order
    std::uint64_t                n_clusters_ = 0;
    bool                         closed_ = false;
};

// Helpers exposed for tests / the reader.
std::string serialise_meta_json(const WriterMetadata& m);

// Random-access reader. Opens, validates magic + footer, then reads clusters
// by index or sequentially via an iterator-like interface.
class Reader {
public:
    explicit Reader(const std::string& path);

    std::uint64_t       n_clusters()  const noexcept { return n_clusters_; }
    const std::string&  meta_json()   const noexcept { return meta_json_; }

    // Read cluster at the given index. Throws on out-of-range or corruption.
    void read_cluster(std::uint64_t idx, ClusterRecord& out);

    // Sequential helper: reads clusters [0..n) in storage order.
    template <typename Fn>
    void for_each(Fn&& fn) {
        ClusterRecord r;
        for (std::uint64_t i = 0; i < n_clusters_; ++i) {
            read_cluster(i, r);
            fn(i, r);
        }
    }

private:
    void parse_block_(const std::vector<std::uint8_t>& block, ClusterRecord& out);

    std::ifstream                   ifs_;
    std::string                     path_;
    std::string                     meta_json_;
    std::vector<std::uint64_t>      offsets_;
    std::uint64_t                   n_clusters_ = 0;
};

// Decode a 2-bit packed sequence (MSB-first, A=0 C=1 G=2 T=3) of length L
// into ACGT ASCII. dst must hold at least L bytes.
void decode_2bit(const std::uint8_t* packed, std::uint32_t L, char* dst);

} // namespace fqdup::clusterfmt
