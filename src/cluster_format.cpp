// fqdup cluster genealogy format writer (.fqcl) — see wiki/cluster-format.md.
// v2: per-block type tag (Full/Tiny/FullZstd); singletons emitted as tinyblock,
// large full blocks zstd-compressed individually (preserves O(1) random access).

#include "fqdup/cluster_format.hpp"

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <zstd.h>

namespace fqdup::clusterfmt {

namespace {

template <typename T>
void write_le(std::ofstream& os, T v) {
    static_assert(std::is_trivially_copyable_v<T>);
    os.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

void write_bytes(std::ofstream& os, const void* p, std::size_t n) {
    if (n) os.write(static_cast<const char*>(p), static_cast<std::streamsize>(n));
}

template <typename T>
void push_pod(std::vector<std::uint8_t>& out, T v) {
    const auto* p = reinterpret_cast<const std::uint8_t*>(&v);
    out.insert(out.end(), p, p + sizeof(v));
}

void push_bytes(std::vector<std::uint8_t>& out, const void* p, std::size_t n) {
    const auto* b = static_cast<const std::uint8_t*>(p);
    out.insert(out.end(), b, b + n);
}

// Singleton-eligible: exactly one member, no genealogy/quality/damage data.
// Damage track is left zero by callers for singletons (no per-cluster damage track yet).
bool is_singleton_tiny_eligible(const ClusterRecord& rec) {
    if (rec.n_members != 1) return false;
    if (!rec.edges.empty()) return false;
    if (!rec.parent_qual.empty()) return false;
    if (rec.member_ids.size() > 1) return false;
    if (rec.parent_seq_len > 0xFFFFu) return false;
    if (rec.member_ids.size() == 1 && rec.member_ids[0].size() > 0xFFFFu) return false;
    for (auto b : rec.damage.term_5) if (b) return false;
    for (auto b : rec.damage.term_3) if (b) return false;
    return true;
}

// Tiny block payload (no leading flags/n_members — implicit):
//   uint64 cluster_id
//   uint16 parent_seq_len
//   uint8[(L+3)/4] parent_seq packed
//   uint16 member_id_len  (0 if absent)
//   uint8[member_id_len] member_id
void serialise_tiny(const ClusterRecord& rec, std::vector<std::uint8_t>& out) {
    out.clear();
    push_pod<std::uint64_t>(out, rec.cluster_id);
    push_pod<std::uint16_t>(out, static_cast<std::uint16_t>(rec.parent_seq_len));
    push_bytes(out, rec.parent_seq.data(), rec.parent_seq.size());
    std::uint16_t id_len = rec.member_ids.empty()
        ? 0u
        : static_cast<std::uint16_t>(rec.member_ids[0].size());
    push_pod<std::uint16_t>(out, id_len);
    if (id_len) push_bytes(out, rec.member_ids[0].data(), id_len);
}

void parse_tiny(const std::uint8_t* p, const std::uint8_t* end, ClusterRecord& out) {
    auto need = [&](std::size_t n) {
        if (p + n > end) throw std::runtime_error("cluster_format: tinyblock truncated");
    };
    need(8); std::memcpy(&out.cluster_id, p, 8); p += 8;
    need(2); std::uint16_t L; std::memcpy(&L, p, 2); p += 2;
    out.parent_seq_len = L;
    const std::size_t pseq_bytes = (L + 3) / 4;
    need(pseq_bytes);
    out.parent_seq.assign(p, p + pseq_bytes); p += pseq_bytes;
    need(2); std::uint16_t id_len; std::memcpy(&id_len, p, 2); p += 2;
    out.member_ids.clear();
    if (id_len) {
        need(id_len);
        out.member_ids.emplace_back(reinterpret_cast<const char*>(p), id_len);
        p += id_len;
    }
    out.flags          = id_len ? kFlagHasMemberIds : 0u;
    out.n_members      = 1;
    out.n_after_damage = 1;
    out.parent_qual.clear();
    out.edges.clear();
    std::memset(out.damage.term_5, 0, kTermWindow);
    std::memset(out.damage.term_3, 0, kTermWindow);
}

// Full payload (no leading type tag — that lives in the wire frame).
void serialise_full(const ClusterRecord& rec, std::vector<std::uint8_t>& out) {
    out.clear();
    push_pod(out, rec.cluster_id);
    push_pod(out, rec.flags);
    push_pod(out, rec.n_members);
    push_pod(out, rec.n_after_damage);
    push_pod(out, rec.parent_seq_len);
    push_bytes(out, rec.parent_seq.data(), rec.parent_seq.size());

    const std::uint32_t qual_len = (rec.flags & kFlagHasQuality)
        ? static_cast<std::uint32_t>(rec.parent_qual.size())
        : 0u;
    push_pod(out, qual_len);
    if (qual_len) push_bytes(out, rec.parent_qual.data(), qual_len);

    const std::uint32_t n_edges = static_cast<std::uint32_t>(rec.edges.size());
    push_pod(out, n_edges);
    if (n_edges) push_bytes(out, rec.edges.data(), n_edges * sizeof(Edge));

    push_bytes(out, rec.damage.term_5, kTermWindow);
    push_bytes(out, rec.damage.term_3, kTermWindow);

    std::uint32_t n_member_ids = 0;
    if (rec.flags & kFlagHasMemberIds)
        n_member_ids = static_cast<std::uint32_t>(rec.member_ids.size());
    push_pod(out, n_member_ids);
    for (std::uint32_t i = 0; i < n_member_ids; ++i) {
        const std::string& id = rec.member_ids[i];
        const std::uint32_t L = static_cast<std::uint32_t>(id.size());
        push_pod(out, L);
        push_bytes(out, id.data(), L);
    }
}

void parse_full(const std::uint8_t* p, const std::uint8_t* end, ClusterRecord& out) {
    auto need = [&](std::size_t n) {
        if (p + n > end) throw std::runtime_error("cluster_format: full block truncated");
    };
    need(8); std::memcpy(&out.cluster_id,     p, 8); p += 8;
    need(4); std::memcpy(&out.flags,          p, 4); p += 4;
    need(4); std::memcpy(&out.n_members,      p, 4); p += 4;
    need(4); std::memcpy(&out.n_after_damage, p, 4); p += 4;
    need(4); std::memcpy(&out.parent_seq_len, p, 4); p += 4;
    const std::size_t pseq_bytes = (out.parent_seq_len + 3) / 4;
    need(pseq_bytes);
    out.parent_seq.assign(p, p + pseq_bytes); p += pseq_bytes;

    need(4); std::uint32_t qual_len; std::memcpy(&qual_len, p, 4); p += 4;
    if (qual_len) {
        need(qual_len);
        out.parent_qual.assign(p, p + qual_len); p += qual_len;
    } else {
        out.parent_qual.clear();
    }

    need(4); std::uint32_t n_edges; std::memcpy(&n_edges, p, 4); p += 4;
    need(static_cast<std::size_t>(n_edges) * sizeof(Edge));
    out.edges.resize(n_edges);
    if (n_edges) std::memcpy(out.edges.data(), p, n_edges * sizeof(Edge));
    p += n_edges * sizeof(Edge);

    need(2 * kTermWindow);
    std::memcpy(out.damage.term_5, p, kTermWindow); p += kTermWindow;
    std::memcpy(out.damage.term_3, p, kTermWindow); p += kTermWindow;

    need(4); std::uint32_t n_ids; std::memcpy(&n_ids, p, 4); p += 4;
    out.member_ids.clear();
    out.member_ids.reserve(n_ids);
    for (std::uint32_t i = 0; i < n_ids; ++i) {
        need(4); std::uint32_t L; std::memcpy(&L, p, 4); p += 4;
        need(L);
        out.member_ids.emplace_back(reinterpret_cast<const char*>(p), L);
        p += L;
    }
}

std::uint32_t crc32_table(std::uint32_t i) {
    std::uint32_t c = i;
    for (int k = 0; k < 8; ++k) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
    return c;
}

std::uint32_t crc32(const void* data, std::size_t n) {
    static std::uint32_t T[256] = {};
    static bool init = false;
    if (!init) { for (std::uint32_t i = 0; i < 256; ++i) T[i] = crc32_table(i); init = true; }
    std::uint32_t c = 0xFFFFFFFFu;
    const auto* p = static_cast<const std::uint8_t*>(data);
    for (std::size_t i = 0; i < n; ++i) c = T[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

} // namespace

std::string serialise_meta_json(const WriterMetadata& m) {
    std::ostringstream o;
    o.precision(6);
    o << "{"
      << "\"tool\":\""         << m.tool         << "\","
      << "\"tool_version\":\"" << m.tool_version << "\","
      << "\"input_fastq\":\""  << m.input_fastq  << "\","
      << "\"n_input_reads\":"  << m.n_input_reads << ","
      << "\"n_clusters\":"     << m.n_clusters    << ","
      << "\"library_type\":\"" << m.library_type << "\","
      << "\"damage\":{"
      << "\"d_max_5\":"  << m.d_max_5  << ","
      << "\"d_max_3\":"  << m.d_max_3  << ","
      << "\"lambda_5\":" << m.lambda_5 << ","
      << "\"lambda_3\":" << m.lambda_3 << ","
      << "\"mask_pos_5\":" << m.mask_pos_5 << ","
      << "\"mask_pos_3\":" << m.mask_pos_3
      << "},"
      << "\"errcor\":{"
      << "\"snp_threshold\":" << m.snp_threshold << ","
      << "\"snp_min_count\":" << m.snp_min_count << ","
      << "\"bucket_cap\":"    << m.bucket_cap
      << "},"
      << "\"loss_counters\":{"
      << "\"bucket_overflow_drops\":"  << m.bucket_overflow_drops  << ","
      << "\"short_interior_skipped\":" << m.short_interior_skipped << ","
      << "\"short_brute_evaluated\":"  << m.short_brute_evaluated  << ","
      << "\"short_brute_found\":"      << m.short_brute_found      << ","
      << "\"short_too_small_skipped\":" << m.short_too_small_skipped
      << "},"
      << "\"block_compression\":\"" << m.block_compression << "\","
      << "\"n_singletons_tinyblock\":" << m.n_singletons_tinyblock
      << "}";
    return o.str();
}

Writer::Writer(const std::string& path, WriterMetadata meta)
    : path_(path), meta_(std::move(meta)) {
    ofs_.open(path_, std::ios::binary | std::ios::trunc);
    if (!ofs_) throw std::runtime_error("cluster_format: cannot open " + path_);
    write_magic_header_();
    write_meta_header_();
}

Writer::~Writer() {
    try { close(); } catch (...) {}
}

void Writer::write_magic_header_() {
    write_le<std::uint32_t>(ofs_, kMagic);
    write_le<std::uint32_t>(ofs_, kVersion);
    // meta_size placeholder; written for real in write_meta_header_().
    write_le<std::uint64_t>(ofs_, 0);
}

void Writer::write_meta_header_() {
    // Reserve a fixed-size meta region (8 KiB) so close() can rewrite the JSON
    // in place with the final n_clusters/n_singletons_tinyblock counts without
    // having to relocate the rest of the file.
    constexpr std::uint64_t kMetaReserve = 8192;
    const std::string json = serialise_meta_json(meta_);
    if (json.size() > kMetaReserve)
        throw std::runtime_error("cluster_format: meta JSON exceeds reserved region");
    const std::streampos cur = ofs_.tellp();
    ofs_.seekp(8, std::ios::beg);
    write_le<std::uint64_t>(ofs_, kMetaReserve);
    ofs_.seekp(cur);
    std::vector<char> region(kMetaReserve, ' ');
    std::memcpy(region.data(), json.data(), json.size());
    region.back() = '\n';
    write_bytes(ofs_, region.data(), region.size());
}

void Writer::write_cluster(const ClusterRecord& rec) {
    if (closed_) throw std::runtime_error("cluster_format: write_cluster after close");

    offsets_.push_back(static_cast<std::uint64_t>(ofs_.tellp()));

    std::vector<std::uint8_t> payload;
    if (is_singleton_tiny_eligible(rec)) {
        serialise_tiny(rec, payload);
        const std::uint32_t wire_size = static_cast<std::uint32_t>(1 + payload.size());
        write_le<std::uint32_t>(ofs_, wire_size);
        write_le<std::uint8_t >(ofs_, kBlockTypeTiny);
        write_bytes(ofs_, payload.data(), payload.size());
        ++meta_.n_singletons_tinyblock;
    } else {
        serialise_full(rec, payload);
        // Compress only when the gain outweighs the per-block frame overhead.
        if (meta_.block_compression == "zstd" && payload.size() >= kZstdMinBytes) {
            const std::size_t bound = ZSTD_compressBound(payload.size());
            std::vector<std::uint8_t> z(bound);
            const std::size_t zsz = ZSTD_compress(z.data(), bound,
                                                  payload.data(), payload.size(), 3);
            if (ZSTD_isError(zsz))
                throw std::runtime_error(std::string("cluster_format: zstd compress failed: ")
                                         + ZSTD_getErrorName(zsz));
            const std::uint32_t wire_size =
                static_cast<std::uint32_t>(1 + 4 + zsz);
            write_le<std::uint32_t>(ofs_, wire_size);
            write_le<std::uint8_t >(ofs_, kBlockTypeFullZstd);
            write_le<std::uint32_t>(ofs_, static_cast<std::uint32_t>(payload.size()));
            write_bytes(ofs_, z.data(), zsz);
        } else {
            const std::uint32_t wire_size = static_cast<std::uint32_t>(1 + payload.size());
            write_le<std::uint32_t>(ofs_, wire_size);
            write_le<std::uint8_t >(ofs_, kBlockTypeFull);
            write_bytes(ofs_, payload.data(), payload.size());
        }
    }
    ++n_clusters_;
}

void Writer::close() {
    if (closed_) return;
    closed_ = true;

    const std::uint64_t n = static_cast<std::uint64_t>(offsets_.size());
    write_bytes(ofs_, offsets_.data(), n * sizeof(std::uint64_t));
    write_le<std::uint64_t>(ofs_, n);
    const std::uint32_t crc = crc32(offsets_.data(), n * sizeof(std::uint64_t));
    write_le<std::uint32_t>(ofs_, crc);
    write_le<std::uint32_t>(ofs_, kMagic);

    // Rewrite meta JSON with final counts. The reserved region is fixed-size,
    // so we just overwrite within [8 + 8, 8 + 8 + kMetaReserve).
    meta_.n_clusters = n_clusters_;
    constexpr std::uint64_t kMetaReserve = 8192;
    const std::string json = serialise_meta_json(meta_);
    if (json.size() <= kMetaReserve) {
        ofs_.seekp(16, std::ios::beg);
        std::vector<char> region(kMetaReserve, ' ');
        std::memcpy(region.data(), json.data(), json.size());
        region.back() = '\n';
        write_bytes(ofs_, region.data(), region.size());
    }
    ofs_.flush();
    ofs_.close();
}

// ── Reader ──────────────────────────────────────────────────────────────────

namespace {

template <typename T>
T read_le_(std::ifstream& is) {
    T v{};
    is.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (!is) throw std::runtime_error("cluster_format: short read");
    return v;
}

} // anon

void decode_2bit(const std::uint8_t* packed, std::uint32_t L, char* dst) {
    static const char A[4] = {'A', 'C', 'G', 'T'};
    for (std::uint32_t i = 0; i < L; ++i) {
        std::uint8_t b = (packed[i >> 2] >> (6 - 2 * (i & 3))) & 0x3u;
        dst[i] = A[b];
    }
}

Reader::Reader(const std::string& path) : path_(path) {
    ifs_.open(path_, std::ios::binary);
    if (!ifs_) throw std::runtime_error("cluster_format: cannot open " + path_);

    auto magic   = read_le_<std::uint32_t>(ifs_);
    auto version = read_le_<std::uint32_t>(ifs_);
    auto msz     = read_le_<std::uint64_t>(ifs_);
    if (magic != kMagic)        throw std::runtime_error("cluster_format: bad magic in " + path_);
    if (version != kVersion)    throw std::runtime_error("cluster_format: unsupported version " + std::to_string(version));
    meta_json_.resize(msz);
    if (msz) ifs_.read(meta_json_.data(), static_cast<std::streamsize>(msz));
    if (!ifs_) throw std::runtime_error("cluster_format: short meta read");

    ifs_.seekg(-16, std::ios::end);
    n_clusters_ = read_le_<std::uint64_t>(ifs_);
    auto crc    = read_le_<std::uint32_t>(ifs_);
    auto magic2 = read_le_<std::uint32_t>(ifs_);
    if (magic2 != kMagic) throw std::runtime_error("cluster_format: bad trailing magic in " + path_);

    offsets_.resize(n_clusters_);
    if (n_clusters_) {
        ifs_.seekg(-(static_cast<std::streamoff>(16 + 8 * n_clusters_)), std::ios::end);
        ifs_.read(reinterpret_cast<char*>(offsets_.data()),
                  static_cast<std::streamsize>(8 * n_clusters_));
        if (!ifs_) throw std::runtime_error("cluster_format: short offsets read");
        std::uint32_t got = crc32(offsets_.data(), n_clusters_ * sizeof(std::uint64_t));
        if (got != crc) throw std::runtime_error("cluster_format: footer CRC mismatch");
    }
}

void Reader::read_cluster(std::uint64_t idx, ClusterRecord& out) {
    if (idx >= n_clusters_)
        throw std::out_of_range("cluster_format: cluster index out of range");
    ifs_.seekg(static_cast<std::streamoff>(offsets_[idx]), std::ios::beg);
    auto wire_size = read_le_<std::uint32_t>(ifs_);
    std::vector<std::uint8_t> wire(wire_size);
    ifs_.read(reinterpret_cast<char*>(wire.data()), wire_size);
    if (!ifs_) throw std::runtime_error("cluster_format: short block read");
    parse_block_(wire, out);
}

void Reader::parse_block_(const std::vector<std::uint8_t>& wire, ClusterRecord& out) {
    if (wire.empty()) throw std::runtime_error("cluster_format: empty block");
    const std::uint8_t type = wire[0];
    const std::uint8_t* p   = wire.data() + 1;
    const std::uint8_t* end = wire.data() + wire.size();

    switch (type) {
        case kBlockTypeFull:
            parse_full(p, end, out);
            break;
        case kBlockTypeTiny:
            parse_tiny(p, end, out);
            break;
        case kBlockTypeFullZstd: {
            if (p + 4 > end) throw std::runtime_error("cluster_format: zstd header truncated");
            std::uint32_t usz; std::memcpy(&usz, p, 4); p += 4;
            std::vector<std::uint8_t> raw(usz);
            const std::size_t got = ZSTD_decompress(raw.data(), usz,
                                                    p, static_cast<std::size_t>(end - p));
            if (ZSTD_isError(got))
                throw std::runtime_error(std::string("cluster_format: zstd decompress failed: ")
                                         + ZSTD_getErrorName(got));
            if (got != usz)
                throw std::runtime_error("cluster_format: zstd size mismatch");
            parse_full(raw.data(), raw.data() + raw.size(), out);
            break;
        }
        default:
            throw std::runtime_error("cluster_format: unknown block_type "
                                     + std::to_string(static_cast<int>(type)));
    }
}

} // namespace fqdup::clusterfmt
