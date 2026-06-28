#pragma once

// Damage-model serialization shared by `profile` (writer) and `split` (reader),
// so split can REUSE the model profile fits instead of re-deriving its own —
// keeping merge/profile/split in sync. Two artifacts:
//
//   JSON  (--damage-json / --damage-json-out): human-readable scalar model
//         (d_max, lambda, bg, pi + verdict, provenance). Small.
//   BIN   (--model-bin / --model-bin-out): the full DamageSplitModel (per-bin
//         LOD float tables + fallback) → split is ZERO-estimation on load.
//
// Both carry an `input_digest` (cheap hash of the first DIGEST_READS read
// sequences). split recomputes it from its own input and REFUSES a model whose
// digest differs — this is the guard against scoring one dataset with another's
// model. `n_reads` is the profile's full read count, carried for diagnostics.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

#include "fqdup/damage_profile.hpp"
#include "fqdup/fastq_common.hpp"

namespace fqdup_model_io {

inline constexpr int64_t DIGEST_READS = 200000;  // reads hashed for the provenance digest

// FNV-1a 64-bit over the first DIGEST_READS read sequences. Deterministic and
// cheap; identical for the same merged FASTQ across profile and split.
inline std::string compute_input_digest(const std::string& path) {
    uint64_t h = 1469598103934665603ULL;          // FNV offset basis
    auto reader = make_fastq_reader(path);
    FastqRecord rec;
    int64_t n = 0;
    while (n < DIGEST_READS && reader->read(rec)) {
        for (char c : rec.seq) { h ^= static_cast<uint8_t>(c); h *= 1099511628211ULL; }
        h ^= static_cast<uint64_t>(rec.seq.size()); h *= 1099511628211ULL;
        ++n;
    }
    h ^= static_cast<uint64_t>(n); h *= 1099511628211ULL;   // fold in count (catches truncation)
    char buf[20]; std::snprintf(buf, sizeof(buf), "%016llx", (unsigned long long)h);
    return std::string(buf);
}

// ───────────────────────── JSON scalar model (v1) ─────────────────────────

inline void write_damage_model_json(const DamageProfile& p, const std::string& path,
                                    int64_t n_reads, const std::string& digest) {
    std::ofstream o(path);
    o << "{\n"
      << "  \"schema\": \"fqdup_damage_model_v1\",\n"
      << "  \"n_reads\": " << n_reads << ",\n"
      << "  \"input_digest\": \"" << digest << "\",\n"
      << "  \"mode_3prime\": "   << static_cast<int>(p.mode_3prime) << ",\n"
      << "  \"d_max_5prime\": "  << p.d_max_5prime  << ",\n"
      << "  \"d_max_3prime\": "  << p.d_max_3prime  << ",\n"
      << "  \"lambda_5prime\": " << p.lambda_5prime << ",\n"
      << "  \"lambda_3prime\": " << p.lambda_3prime << ",\n"
      << "  \"bg_5_tc\": "       << p.bg_5_tc       << ",\n"
      << "  \"bg_3_channel\": "  << p.bg_3_channel  << ",\n"
      << "  \"mask_threshold\": "<< p.mask_threshold<< ",\n"
      << "  \"ss_mode\": "        << (p.ss_mode ? 1 : 0) << ",\n"
      << "  \"d_cpg_5prime\": "   << p.d_cpg_5prime    << ",\n"
      << "  \"d_noncpg_5prime\": "<< p.d_noncpg_5prime << ",\n"
      << "  \"mixture_converged\": " << (p.mixture_converged ? 1 : 0) << ",\n"
      << "  \"mixture_d_damaged\": "  << p.mixture_d_damaged  << ",\n"
      << "  \"mixture_pi_damaged\": " << p.mixture_pi_damaged << ",\n"
      << "  \"pi_detected\": " << (p.pi_detected ? 1 : 0) << ",\n"
      << "  \"pi_point\": " << p.pi_point << ",\n"
      << "  \"pi_lo\": "    << p.pi_lo    << ",\n"
      << "  \"pi_hi\": "    << p.pi_hi    << ",\n"
      << "  \"pi_state\": " << static_cast<int>(p.pi_state) << "\n"
      << "}\n";
}

inline double json_num_(const std::string& s, const char* key, double dflt) {
    std::string k = std::string("\"") + key + "\"";
    auto pos = s.find(k);
    if (pos == std::string::npos) return dflt;
    pos = s.find(':', pos + k.size());
    if (pos == std::string::npos) return dflt;
    return std::strtod(s.c_str() + pos + 1, nullptr);
}

inline std::string json_str_(const std::string& s, const char* key) {
    std::string k = std::string("\"") + key + "\"";
    auto pos = s.find(k);
    if (pos == std::string::npos) return "";
    pos = s.find(':', pos + k.size());
    if (pos == std::string::npos) return "";
    auto q1 = s.find('"', pos); if (q1 == std::string::npos) return "";
    auto q2 = s.find('"', q1 + 1); if (q2 == std::string::npos) return "";
    return s.substr(q1 + 1, q2 - q1 - 1);
}

// Returns false if unreadable. Fills p and out-params; recomputes mask.
inline bool load_damage_model_json(const std::string& path, DamageProfile& p,
                                   int64_t& n_reads, std::string& digest) {
    std::ifstream f(path);
    if (!f) return false;
    std::stringstream ss; ss << f.rdbuf();
    const std::string s = ss.str();
    p.enabled      = true;   // model came from a real fit; required by the LSD rebuild path
    n_reads        = static_cast<int64_t>(json_num_(s, "n_reads", -1));
    digest         = json_str_(s, "input_digest");
    p.mode_3prime  = static_cast<Damage3PrimeMode>((int)json_num_(s, "mode_3prime", 0));
    p.d_max_5prime = json_num_(s, "d_max_5prime", 0.0);
    p.d_max_3prime = json_num_(s, "d_max_3prime", 0.0);
    p.lambda_5prime= json_num_(s, "lambda_5prime", 0.5);
    p.lambda_3prime= json_num_(s, "lambda_3prime", 0.5);
    p.bg_5_tc      = json_num_(s, "bg_5_tc", 0.5);
    p.bg_3_channel = json_num_(s, "bg_3_channel", 0.5);
    p.mask_threshold = json_num_(s, "mask_threshold", 0.05);
    p.ss_mode        = json_num_(s, "ss_mode", 0.0) != 0.0;
    p.d_cpg_5prime   = json_num_(s, "d_cpg_5prime", 0.0);
    p.d_noncpg_5prime= json_num_(s, "d_noncpg_5prime", 0.0);
    p.mixture_converged  = json_num_(s, "mixture_converged", 0.0) != 0.0;
    p.mixture_d_damaged  = json_num_(s, "mixture_d_damaged", 0.0);
    p.mixture_pi_damaged = json_num_(s, "mixture_pi_damaged", 0.0);
    p.pi_detected  = json_num_(s, "pi_detected", 0.0) != 0.0;
    p.pi_point     = json_num_(s, "pi_point", -1.0);
    p.pi_lo        = json_num_(s, "pi_lo", -1.0);
    p.pi_hi        = json_num_(s, "pi_hi", -1.0);
    p.pi_state     = static_cast<taph::DamageConfidence>((int)json_num_(s, "pi_state", 0));
    p.populate_mask_from_model();
    return true;
}

// ──────────────────── binary full split model (v2) ────────────────────────

inline constexpr char     DMGMODEL_MAGIC[4] = {'F','Q','D','M'};
inline constexpr uint32_t DMGMODEL_VERSION  = 2;

inline void write_split_model_bin(const DamageSplitModel& m, const std::string& path,
                                  int64_t n_reads, const std::string& digest) {
    std::ofstream o(path, std::ios::binary);
    auto w = [&](const void* p, size_t n){ o.write(reinterpret_cast<const char*>(p), n); };
    w(DMGMODEL_MAGIC, 4);
    uint32_t ver = DMGMODEL_VERSION; w(&ver, 4);
    uint64_t nr = (uint64_t)n_reads; w(&nr, 8);
    uint32_t dl = (uint32_t)digest.size(); w(&dl, 4); w(digest.data(), dl);
    int32_t st = (int32_t)m.fallback.pi_state; w(&st, 4);
    double piv[3] = {m.fallback.pi_point, m.fallback.pi_lo, m.fallback.pi_hi}; w(piv, sizeof(piv));
    double sc[8] = {m.fallback.d_max_5prime, m.fallback.d_max_3prime,
                    m.fallback.lambda_5prime, m.fallback.lambda_3prime,
                    m.fallback.bg_5_tc, m.fallback.bg_3_channel,
                    m.fallback.mask_threshold, m.fallback.mixture_pi_damaged}; w(sc, sizeof(sc));
    int32_t mode = (int32_t)m.fallback.mode_3prime; w(&mode, 4);
    int32_t nb = (int32_t)m.bins.size(); w(&nb, 4);
    for (const auto& b : m.bins) {
        int32_t lo=b.lo, hi=b.hi; w(&lo,4); w(&hi,4);
        uint8_t ss = b.ss_mode ? 1 : 0; w(&ss,1);
        w(b.lod5_T.data(),     sizeof(b.lod5_T));
        w(b.lod5_C.data(),     sizeof(b.lod5_C));
        w(b.lod3_T.data(),     sizeof(b.lod3_T));
        w(b.lod3_C.data(),     sizeof(b.lod3_C));
        w(b.lod5_T_cpg.data(), sizeof(b.lod5_T_cpg));
        w(b.lod5_C_cpg.data(), sizeof(b.lod5_C_cpg));
    }
}

// Returns false on unreadable / bad magic / wrong version.
inline bool read_split_model_bin(const std::string& path, DamageSplitModel& m,
                                 int64_t& n_reads, std::string& digest) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    auto r = [&](void* p, size_t n){ f.read(reinterpret_cast<char*>(p), n); };
    char mg[4]; r(mg,4); if (std::memcmp(mg, DMGMODEL_MAGIC, 4) != 0) return false;
    uint32_t ver; r(&ver,4); if (ver != DMGMODEL_VERSION) return false;
    uint64_t nr; r(&nr,8); n_reads = (int64_t)nr;
    uint32_t dl; r(&dl,4); digest.resize(dl); if (dl) r(&digest[0], dl);
    int32_t st; r(&st,4); m.fallback.pi_state = (taph::DamageConfidence)st;
    double piv[3]; r(piv,sizeof(piv));
    m.fallback.pi_point=piv[0]; m.fallback.pi_lo=piv[1]; m.fallback.pi_hi=piv[2];
    double sc[8]; r(sc,sizeof(sc));
    m.fallback.d_max_5prime=sc[0]; m.fallback.d_max_3prime=sc[1];
    m.fallback.lambda_5prime=sc[2]; m.fallback.lambda_3prime=sc[3];
    m.fallback.bg_5_tc=sc[4]; m.fallback.bg_3_channel=sc[5];
    m.fallback.mask_threshold=sc[6]; m.fallback.mixture_pi_damaged=sc[7];
    int32_t mode; r(&mode,4); m.fallback.mode_3prime = (Damage3PrimeMode)mode;
    m.fallback.populate_mask_from_model();
    int32_t nb; r(&nb,4); if (nb < 0 || nb > 100000) return false;
    m.bins.resize(nb);
    for (auto& b : m.bins) {
        r(&b.lo,4); r(&b.hi,4); uint8_t ss; r(&ss,1); b.ss_mode = (ss != 0);
        r(b.lod5_T.data(),     sizeof(b.lod5_T));
        r(b.lod5_C.data(),     sizeof(b.lod5_C));
        r(b.lod3_T.data(),     sizeof(b.lod3_T));
        r(b.lod3_C.data(),     sizeof(b.lod3_C));
        r(b.lod5_T_cpg.data(), sizeof(b.lod5_T_cpg));
        r(b.lod5_C_cpg.data(), sizeof(b.lod5_C_cpg));
    }
    return static_cast<bool>(f);
}

}  // namespace fqdup_model_io
