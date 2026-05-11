// fqdup view — inspect .fqcl cluster genealogy files.
// Modes: header summary (default), single-cluster tree/staircase, sibling bundling.

#include "view_html_assets.hpp"
#include "fqdup/cluster_format.hpp"
#include "fqdup/version.hpp"
#include "fqdup/bundle_key.hpp"
#include "fqdup/fastq_types.hpp"

#include <xxhash.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cf = fqdup::clusterfmt;

namespace {

using fqdup::bundlekey::from_cluster;

// 5'-CT or 3'-GA terminal damage signature (MASK_POSITIONS=15 default).
constexpr std::uint32_t kDamageMaskPositions = 15;

// ---- JSON helpers -----------------------------------------------------------

void json_escape(std::ostream& os, std::string_view s) {
    os << '"';
    for (unsigned char c : s) {
        switch (c) {
            case '"':  os << "\\\""; break;
            case '\\': os << "\\\\"; break;
            case '\b': os << "\\b";  break;
            case '\f': os << "\\f";  break;
            case '\n': os << "\\n";  break;
            case '\r': os << "\\r";  break;
            case '\t': os << "\\t";  break;
            // Escape < and > so JSON embedded in <script> blocks can't inject closing tags
            case '<':  os << "\\u003c"; break;
            case '>':  os << "\\u003e"; break;
            case '&':  os << "\\u0026"; break;
            default:
                if (c < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    os << buf;
                } else {
                    os << static_cast<char>(c);
                }
        }
    }
    os << '"';
}

std::string json_escape_str(std::string_view s) {
    std::ostringstream o;
    json_escape(o, s);
    return o.str();
}

// Emit envelope opener up to and including: ,"data":
// Caller writes the data value, then the closing "}".
void emit_envelope_open(std::ostream& os,
                        const cf::Reader& rd,
                        const std::string& path,
                        std::string_view mode,
                        std::string_view args_json) {
    os << "{\"schema\":\"fqdup.view.v1\",\"fqcl\":{\"path\":";
    json_escape(os, path);
    os << ",\"n_clusters\":" << rd.n_clusters()
       << ",\"metadata\":" << (rd.meta_json().empty() ? std::string("{}") : rd.meta_json())
       << "},\"query\":{\"mode\":";
    json_escape(os, mode);
    os << ",\"args\":" << args_json << "},\"data\":";
}

void emit_envelope_close(std::ostream& os) { os << "}"; }

// ---- Edit walking -----------------------------------------------------------

// For a member node, walk up the edges to root, collect (pos -> base) edits.
// Edges store from_node→to_node; the parent of node N is found by scanning.
void apply_edits(const cf::ClusterRecord& r, std::uint32_t node, char* seq) {
    // parent_edge_of[to_node] = edge index
    std::unordered_map<std::uint32_t, std::uint32_t> parent_edge;
    parent_edge.reserve(r.edges.size());
    for (std::uint32_t i = 0; i < r.edges.size(); ++i)
        parent_edge[r.edges[i].to_node] = i;

    std::uint32_t cur = node;
    while (cur != 0) {
        auto it = parent_edge.find(cur);
        if (it == parent_edge.end()) break;
        const auto& e = r.edges[it->second];
        static const char A[4] = {'A','C','G','T'};
        seq[e.pos] = A[e.to_base & 0x3u];
        cur = e.from_node;
    }
}

// ---- Plaintext printers (unchanged) -----------------------------------------

void print_header(const cf::Reader& rd) {
    std::cout << "file              : " << rd.n_clusters() << " clusters\n";
    std::cout << "metadata          : " << rd.meta_json() << "\n";
}

void print_tree(const cf::ClusterRecord& r) {
    std::cout << "cluster " << r.cluster_id
              << "  n_members=" << r.n_members
              << "  parent_len=" << r.parent_seq_len
              << "  edges=" << r.edges.size() << "\n";

    std::vector<char> pbuf(r.parent_seq_len);
    cf::decode_2bit(r.parent_seq.data(), r.parent_seq_len, pbuf.data());
    std::cout << "  parent: " << std::string(pbuf.data(), pbuf.size()) << "\n";

    // children-by-parent index
    std::map<std::uint32_t, std::vector<std::uint32_t>> children;
    for (std::uint32_t i = 0; i < r.edges.size(); ++i)
        children[r.edges[i].from_node].push_back(i);

    static const char A[4] = {'A','C','G','T'};
    auto walk = [&](auto&& self, std::uint32_t node, int depth) -> void {
        auto it = children.find(node);
        if (it == children.end()) return;
        for (auto eidx : it->second) {
            const auto& e = r.edges[eidx];
            std::cout << "  " << std::string(depth * 2, ' ')
                      << "└─ node " << e.to_node
                      << " : pos " << e.pos
                      << " " << A[e.from_base & 3] << "→" << A[e.to_base & 3]
                      << " (n_reads=" << e.n_reads << ")\n";
            self(self, e.to_node, depth + 1);
        }
    };
    walk(walk, 0, 0);
}

void print_staircase(const cf::ClusterRecord& r) {
    if (r.parent_seq_len == 0) { std::cout << "(empty)\n"; return; }
    std::vector<char> parent(r.parent_seq_len);
    cf::decode_2bit(r.parent_seq.data(), r.parent_seq_len, parent.data());

    // Collect every node that appears (root + every to_node)
    std::vector<std::uint32_t> nodes{0};
    for (const auto& e : r.edges) nodes.push_back(e.to_node);
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());

    std::printf("# cluster=%lu  parent_len=%u  members=%u  nodes=%zu\n",
                static_cast<unsigned long>(r.cluster_id),
                r.parent_seq_len, r.n_members, nodes.size());
    std::printf("#                   ");
    for (std::uint32_t i = 0; i < r.parent_seq_len; ++i)
        std::putchar(parent[i]);
    std::printf("  (parent)\n");

    std::vector<char> rbuf(r.parent_seq_len);
    for (auto n : nodes) {
        std::memcpy(rbuf.data(), parent.data(), r.parent_seq_len);
        if (n != 0) apply_edits(r, n, rbuf.data());
        std::printf("node %6u : ", n);
        // dot grid: '.' where matches parent, base letter where differs
        for (std::uint32_t i = 0; i < r.parent_seq_len; ++i) {
            char c = (rbuf[i] == parent[i]) ? '.' : rbuf[i];
            std::putchar(c);
        }
        std::putchar('\n');
    }
}

void print_bundles(cf::Reader& rd, int end_k, std::size_t min_size) {
    std::map<std::uint64_t, std::vector<std::uint64_t>> by_key;
    cf::ClusterRecord r;
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        by_key[from_cluster(r, end_k)].push_back(i);
    }
    std::printf("bundle_key       \tn_clusters\tmember_clusters\n");
    std::size_t shown = 0, skipped = 0;
    for (auto& kv : by_key) {
        if (kv.second.size() < min_size) { skipped++; continue; }
        std::printf("%016lx\t%zu\t", static_cast<unsigned long>(kv.first), kv.second.size());
        for (std::size_t i = 0; i < kv.second.size(); ++i) {
            if (i) std::putchar(',');
            std::printf("%lu", static_cast<unsigned long>(kv.second[i]));
        }
        std::putchar('\n');
        shown++;
    }
    std::fprintf(stderr, "%zu bundles shown, %zu hidden (size<%zu)\n",
                 shown, skipped, min_size);
}

void print_bundle_staircase(cf::Reader& rd, std::uint64_t key, int end_k) {
    cf::ClusterRecord r;
    bool first = true;
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        if (from_cluster(r, end_k) != key) continue;
        if (!first) std::printf("\n");
        first = false;
        print_staircase(r);
    }
    if (first) std::fprintf(stderr, "no clusters in bundle %016lx\n",
                            static_cast<unsigned long>(key));
}

void print_top(cf::Reader& rd, std::size_t k, bool by_edges) {
    struct Hit { std::uint64_t id; std::uint32_t n_members; std::uint32_t n_edges; std::uint32_t parent_len; };
    std::vector<Hit> top;
    cf::ClusterRecord r;
    auto cmp = [&](const Hit& a, const Hit& b) {
        if (by_edges) return a.n_edges > b.n_edges;
        return a.n_members > b.n_members;
    };
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        Hit h{r.cluster_id,
              r.n_members,
              static_cast<std::uint32_t>(r.edges.size()),
              r.parent_seq_len};
        if (top.size() < k) {
            top.push_back(h);
            std::push_heap(top.begin(), top.end(), cmp);
        } else {
            if (cmp(h, top.front())) {
                std::pop_heap(top.begin(), top.end(), cmp);
                top.back() = h;
                std::push_heap(top.begin(), top.end(), cmp);
            }
        }
        if ((i % 5000000) == 0)
            std::fprintf(stderr, "scanned %lu / %lu\n",
                         static_cast<unsigned long>(i),
                         static_cast<unsigned long>(rd.n_clusters()));
    }
    std::sort(top.begin(), top.end(), cmp);
    std::printf("cluster_id\tn_members\tn_edges\tparent_len\n");
    for (const auto& h : top)
        std::printf("%lu\t%u\t%u\t%u\n",
                    static_cast<unsigned long>(h.id),
                    h.n_members, h.n_edges, h.parent_len);
}

// ---- JSON emitters ----------------------------------------------------------

// Emit the inner cluster object: {"cluster_id":..,"n_members":..,...,"nodes":[...],"edges":[...]}
void emit_cluster_object(std::ostream& os, const cf::ClusterRecord& r) {
    std::vector<char> pbuf(r.parent_seq_len);
    if (r.parent_seq_len > 0)
        cf::decode_2bit(r.parent_seq.data(), r.parent_seq_len, pbuf.data());

    os << "{\"cluster_id\":" << r.cluster_id
       << ",\"n_members\":" << r.n_members
       << ",\"n_edges\":" << r.edges.size()
       << ",\"parent_len\":" << r.parent_seq_len
       << ",\"parent_seq\":";
    json_escape(os, std::string_view(pbuf.data(), pbuf.size()));

    // Build children index for BFS with depth/parent tracking
    std::map<std::uint32_t, std::vector<std::uint32_t>> children;
    for (std::uint32_t i = 0; i < r.edges.size(); ++i)
        children[r.edges[i].from_node].push_back(i);

    // BFS over the DAG, visiting each to_node exactly once under its first
    // encountered parent. The .fqcl genealogy is a DAG (a node may have
    // multiple candidate incoming edges); for tree-style viz we pick the
    // first parent and let the full edge list expose the rest.
    struct NodeRec { std::uint32_t id; std::int64_t parent; int depth; bool is_parent; };
    std::vector<NodeRec> bfs;
    std::unordered_set<std::uint32_t> visited;
    bfs.push_back({0, -1, 0, true});
    visited.insert(0);
    for (std::size_t head = 0; head < bfs.size(); ++head) {
        auto cur = bfs[head];
        auto it = children.find(cur.id);
        if (it == children.end()) continue;
        for (auto eidx : it->second) {
            const auto& e = r.edges[eidx];
            if (!visited.insert(e.to_node).second) continue;
            bfs.push_back({e.to_node, static_cast<std::int64_t>(cur.id), cur.depth + 1, false});
        }
    }

    os << ",\"nodes\":[";
    for (std::size_t i = 0; i < bfs.size(); ++i) {
        if (i) os << ",\n    ";
        const auto& n = bfs[i];
        if (n.is_parent) {
            os << "{\"id\":" << n.id << ",\"is_parent\":true}";
        } else {
            os << "{\"id\":" << n.id
               << ",\"parent\":" << n.parent
               << ",\"depth\":" << n.depth << "}";
        }
    }
    os << "]";

    static const char A[4] = {'A','C','G','T'};
    os << ",\"edges\":[";
    for (std::size_t i = 0; i < r.edges.size(); ++i) {
        if (i) os << ",\n    ";
        const auto& e = r.edges[i];
        char ref = A[e.from_base & 3];
        char alt = A[e.to_base   & 3];
        bool damage_like =
            (ref == 'C' && alt == 'T' && e.pos < kDamageMaskPositions) ||
            (ref == 'G' && alt == 'A' && r.parent_seq_len >= kDamageMaskPositions
                && e.pos >= r.parent_seq_len - kDamageMaskPositions);
        os << "{\"from\":" << e.from_node
           << ",\"to\":" << e.to_node
           << ",\"pos\":" << e.pos
           << ",\"ref\":\"" << ref << "\""
           << ",\"alt\":\"" << alt << "\""
           << ",\"n_reads\":" << e.n_reads
           << ",\"damage_like\":" << (damage_like ? "true" : "false")
           << (std::isnan(e.score)
                   ? ",\"score\":null,\"score_evaluated\":false"
                   : ",\"score\":" + std::to_string(e.score) + ",\"score_evaluated\":true")
           << "}";
    }
    os << "]}";
}

// Emit staircase mutation rows for a cluster.
void emit_staircase_rows(std::ostream& os, const cf::ClusterRecord& r) {
    os << "{\"rows\":[";
    if (r.parent_seq_len == 0) { os << "]}"; return; }

    std::vector<char> parent(r.parent_seq_len);
    cf::decode_2bit(r.parent_seq.data(), r.parent_seq_len, parent.data());

    std::vector<std::uint32_t> nodes{0};
    for (const auto& e : r.edges) nodes.push_back(e.to_node);
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());

    std::vector<char> rbuf(r.parent_seq_len);
    for (std::size_t ni = 0; ni < nodes.size(); ++ni) {
        if (ni) os << ",\n    ";
        std::uint32_t n = nodes[ni];
        std::memcpy(rbuf.data(), parent.data(), r.parent_seq_len);
        if (n != 0) apply_edits(r, n, rbuf.data());

        os << "{\"node\":" << n << ",\"mutations\":[";
        bool first_mut = true;
        for (std::uint32_t i = 0; i < r.parent_seq_len; ++i) {
            if (rbuf[i] == parent[i]) continue;
            if (!first_mut) os << ",";
            first_mut = false;
            os << "{\"pos\":" << i << ",\"base\":\"" << rbuf[i] << "\"}";
        }
        os << "]}";
    }
    os << "]}";
}

void emit_header_json(std::ostream& os, const cf::Reader& rd, const std::string& path) {
    emit_envelope_open(os, rd, path, "header", "{}");
    os << "{}";
    emit_envelope_close(os);
    os << "\n";
}

void emit_cluster_json(std::ostream& os, const cf::Reader& rd,
                       const std::string& path, std::uint64_t cluster_id,
                       const cf::ClusterRecord& r) {
    std::ostringstream args;
    args << "{\"cluster_id\":" << cluster_id << "}";
    emit_envelope_open(os, rd, path, "cluster", args.str());
    os << "{\"cluster\":";
    emit_cluster_object(os, r);
    os << "}";
    emit_envelope_close(os);
    os << "\n";
}

void emit_staircase_json(std::ostream& os, const cf::Reader& rd,
                         const std::string& path, std::uint64_t cluster_id,
                         const cf::ClusterRecord& r) {
    std::ostringstream args;
    args << "{\"cluster_id\":" << cluster_id << "}";
    emit_envelope_open(os, rd, path, "staircase", args.str());
    os << "{\"cluster\":";
    emit_cluster_object(os, r);
    os << ",\"staircase\":";
    emit_staircase_rows(os, r);
    os << "}";
    emit_envelope_close(os);
    os << "\n";
}

void emit_bundles_json(std::ostream& os, cf::Reader& rd,
                       const std::string& path, int end_k, std::size_t min_size) {
    std::map<std::uint64_t, std::vector<std::uint64_t>> by_key;
    cf::ClusterRecord r;
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        by_key[from_cluster(r, end_k)].push_back(i);
    }

    std::ostringstream args;
    args << "{\"end_k\":" << end_k << ",\"min_size\":" << min_size << "}";
    emit_envelope_open(os, rd, path, "bundle", args.str());

    os << "{\"bundles\":[";
    bool first = true;
    std::size_t skipped = 0;
    for (auto& kv : by_key) {
        if (kv.second.size() < min_size) { skipped++; continue; }
        if (!first) os << ",\n    ";
        first = false;
        char hex[32];
        std::snprintf(hex, sizeof(hex), "%016lx", static_cast<unsigned long>(kv.first));
        os << "{\"bundle_key\":\"" << hex << "\""
           << ",\"n_clusters\":" << kv.second.size()
           << ",\"clusters\":[";
        for (std::size_t i = 0; i < kv.second.size(); ++i) {
            if (i) os << ",";
            os << kv.second[i];
        }
        os << "]}";
    }
    os << "],\"skipped_below_min_size\":" << skipped << "}";
    emit_envelope_close(os);
    os << "\n";
}

void emit_bundle_staircase_json(std::ostream& os, cf::Reader& rd,
                                const std::string& path,
                                std::uint64_t key, int end_k) {
    char hex[32];
    std::snprintf(hex, sizeof(hex), "%016lx", static_cast<unsigned long>(key));

    std::ostringstream args;
    args << "{\"bundle_key\":\"" << hex << "\",\"end_k\":" << end_k << "}";
    emit_envelope_open(os, rd, path, "bundle_staircase", args.str());

    os << "{\"bundle_key\":\"" << hex << "\",\"clusters\":[";
    cf::ClusterRecord r;
    bool first = true;
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        if (from_cluster(r, end_k) != key) continue;
        if (!first) os << ",\n    ";
        first = false;
        os << "{\"cluster\":";
        emit_cluster_object(os, r);
        os << ",\"staircase\":";
        emit_staircase_rows(os, r);
        os << "}";
    }
    os << "]}";
    emit_envelope_close(os);
    os << "\n";
}

void emit_top_json(std::ostream& os, cf::Reader& rd, const std::string& path,
                   std::size_t k, bool by_edges) {
    struct Hit {
        std::uint64_t id;
        std::uint32_t n_members;
        std::uint32_t n_edges;
        std::uint32_t parent_len;
        std::uint64_t bundle_key;
    };
    std::vector<Hit> top;
    cf::ClusterRecord r;
    auto cmp = [&](const Hit& a, const Hit& b) {
        if (by_edges) return a.n_edges > b.n_edges;
        return a.n_members > b.n_members;
    };
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        rd.read_cluster(i, r);
        Hit h{r.cluster_id,
              r.n_members,
              static_cast<std::uint32_t>(r.edges.size()),
              r.parent_seq_len,
              from_cluster(r, 16)};
        if (top.size() < k) {
            top.push_back(h);
            std::push_heap(top.begin(), top.end(), cmp);
        } else {
            if (cmp(h, top.front())) {
                std::pop_heap(top.begin(), top.end(), cmp);
                top.back() = h;
                std::push_heap(top.begin(), top.end(), cmp);
            }
        }
        if ((i % 5000000) == 0)
            std::fprintf(stderr, "scanned %lu / %lu\n",
                         static_cast<unsigned long>(i),
                         static_cast<unsigned long>(rd.n_clusters()));
    }
    std::sort(top.begin(), top.end(), cmp);

    std::ostringstream args;
    args << "{\"k\":" << k << "}";
    const char* mode = by_edges ? "top_edges" : "top_members";
    emit_envelope_open(os, rd, path, mode, args.str());
    os << "{\"by\":\"" << (by_edges ? "edges" : "members") << "\""
       << ",\"k\":" << k << ",\"clusters\":[";
    for (std::size_t i = 0; i < top.size(); ++i) {
        if (i) os << ",\n    ";
        const auto& h = top[i];
        char hex[32];
        std::snprintf(hex, sizeof(hex), "%016lx", static_cast<unsigned long>(h.bundle_key));
        os << "{\"cluster_id\":" << h.id
           << ",\"n_members\":" << h.n_members
           << ",\"n_edges\":" << h.n_edges
           << ",\"parent_len\":" << h.parent_len
           << ",\"bundle_key\":\"" << hex << "\"}";
    }
    os << "]}";
    emit_envelope_close(os);
    os << "\n";
}

// ── static HTML export ────────────────────────────────────────────────────────
// Produces a single self-contained HTML file with embedded CSS, JS (viz4), and
// all cluster data inlined as window.FQCL_DATA.  No server needed to view it.
void emit_html_static(cf::Reader& rd, const std::string& path,
                      const std::string& out_path, std::size_t n_top,
                      bool by_edges = false) {
    // 1. header JSON
    std::ostringstream hdr;
    emit_header_json(hdr, rd, path);

    // 2. top-N clusters list + per-cluster JSON
    //    Collect top N by member count (reuse same heap logic as emit_top_json)
    struct Hit {
        std::uint64_t id;
        std::uint32_t n_members;
        std::uint32_t n_edges;
        std::uint32_t parent_len;
        std::uint64_t bundle_key;
    };
    std::vector<Hit> top;
    cf::ClusterRecord r;
    auto cmp = [by_edges](const Hit& a, const Hit& b){
        return by_edges ? a.n_edges > b.n_edges : a.n_members > b.n_members;
    };
    for (std::uint64_t i = 0; i < rd.n_clusters(); ++i) {
        cf::ClusterRecord rec;
        rd.read_cluster(i, rec);
        Hit h{ rec.cluster_id,
               rec.n_members,
               static_cast<std::uint32_t>(rec.edges.size()),
               rec.parent_seq_len,
               from_cluster(rec, 16) };
        if (top.size() < n_top) {
            top.push_back(h);
            std::push_heap(top.begin(), top.end(), cmp);
        } else if (cmp(h, top.front())) {
            std::pop_heap(top.begin(), top.end(), cmp);
            top.back() = h;
            std::push_heap(top.begin(), top.end(), cmp);
        }
        if ((i % 5000000) == 0)
            std::fprintf(stderr, "scanned %lu / %lu\r",
                         static_cast<unsigned long>(i),
                         static_cast<unsigned long>(rd.n_clusters()));
    }
    std::sort(top.begin(), top.end(), [](const Hit& a, const Hit& b){
        return a.n_members > b.n_members;
    });

    // top_members JSON (reuse envelope format viz4 expects)
    std::ostringstream tm;
    tm << "{\"schema\":\"fqdup.view.v1\",\"fqcl\":{\"path\":";
    json_escape(tm, path);
    tm << "},\"query\":{\"mode\":\"top_members\",\"args\":{\"k\":" << n_top << "}},"
       << "\"data\":{\"by\":\"members\",\"k\":" << n_top << ",\"clusters\":[";
    for (std::size_t i = 0; i < top.size(); ++i) {
        if (i) tm << ",";
        char hex[32];
        std::snprintf(hex, sizeof(hex), "%016lx",
                      static_cast<unsigned long>(top[i].bundle_key));
        tm << "{\"cluster_id\":" << top[i].id
           << ",\"n_members\":" << top[i].n_members
           << ",\"n_edges\":" << top[i].n_edges
           << ",\"parent_len\":" << top[i].parent_len
           << ",\"bundle_key\":\"" << hex << "\"}";
    }
    tm << "]}}";

    // per-cluster JSON map: second pass reading only the top clusters
    std::map<std::uint64_t, std::string> cluster_jsons;
    for (const auto& h : top) {
        cf::ClusterRecord rec;
        rd.read_cluster(h.id, rec);
        std::ostringstream cj;
        cj << "{\"schema\":\"fqdup.view.v1\",\"fqcl\":{\"path\":";
        json_escape(cj, path);
        cj << ",\"n_clusters\":" << rd.n_clusters();
        cj << ",\"metadata\":" << (rd.meta_json().empty() ? "{}" : rd.meta_json());
        cj << "},\"query\":{\"mode\":\"cluster\",\"args\":{\"cluster_id\":"
           << rec.cluster_id << "}},\"data\":{\"cluster\":";
        emit_cluster_object(cj, rec);
        cj << "}}";
        cluster_jsons[rec.cluster_id] = cj.str();
    }

    // 3. write HTML
    std::ofstream out(out_path);
    if (!out) throw std::runtime_error("cannot open output: " + out_path);

    out << R"(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>fqcl · cluster genealogy</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
)";
    out << fqdup_html::CSS;
    out << R"(
</style>
</head>
<body data-theme="light">
  <header id="topbar">
    <div class="brand"><span class="b1">fqcl</span><span class="b2">cluster genealogy · v)" FQDUP_VERSION R"(</span></div>
    <div id="ds-meta"></div>
    <div class="topbar-actions">
      <button id="theme-toggle" title="theme">◐</button>
      <button id="export-svg" title="export panel as SVG">↓ svg</button>
    </div>
  </header>
  <main>
    <aside id="nav">
      <div class="nav-head">
        <input id="nav-filter" placeholder="filter cluster id…">
        <span id="nav-count"></span>
      </div>
      <div id="size-hist" class="size-hist"></div>
      <div id="nav-list"></div>
    </aside>
    <section id="sheet">
      <div id="cluster-head">
        <div class="ch-id">
          <span class="ch-label">cluster</span>
          <span class="ch-num" id="cluster-id">—</span>
          <div class="ch-stats" id="cluster-stats"></div>
        </div>
        <div id="qc-chips"></div>
      </div>
      <div class="grid">
        <div class="panel wide">
          <div class="panel-head"><span class="ph-l1">parent sequence</span></div>
          <div id="parent-tracks"></div>
        </div>
        <div class="panel wide" id="heatmap-panel">
          <div class="panel-head"><span class="ph-l1">damage fingerprint</span>
            <span class="ph-l2 muted" style="font-size:11px"> · rows=depth cols=position color=mutation class</span></div>
          <div id="heatmap-canvas"></div>
        </div>
        <div class="panel wide">
          <div class="panel-head">
            <span class="ph-l1">genealogy tree</span>
            <span class="ph-l2 mono" id="gn-rowcount"></span>
            <div class="panel-controls">
              <div class="ctrl-group"><span>color</span>
                <button data-color="class" class="active">class</button>
                <button data-color="depth">depth</button>
              </div>
            </div>
          </div>
          <div id="tree-canvas"></div>
        </div>
        <div class="panel">
          <div class="panel-head"><span class="ph-l1">mutation class</span></div>
          <div id="class-bar"></div>
          <div id="class-legend" style="display:flex;flex-wrap:wrap;gap:10px;margin-top:8px;font-size:11px"></div>
        </div>
        <div class="panel">
          <div class="panel-head"><span class="ph-l1">evidence</span></div>
          <div id="evidence"></div>
        </div>
        <div class="panel wide">
          <div class="panel-head"><span class="ph-l1">mutation ledger</span></div>
          <div class="ledger-row head">
            <div></div><div>class</div><div>pos</div><div>change</div>
            <div>zone</div><div>reason</div><div>reads</div><div>fraction</div>
          </div>
          <div id="ledger"></div>
        </div>
      </div>
    </section>
  </main>
  <div id="tooltip"></div>
)";

    // embed data
    out << "<script>\nwindow.FQCL_DATA = {\n";
    out << "  header: " << hdr.str() << ",\n";
    out << "  top_members: " << tm.str() << ",\n";
    out << "  clusters: {\n";
    bool first = true;
    for (const auto& [cid, json] : cluster_jsons) {
        if (!first) out << ",\n";
        out << "    \"" << cid << "\": " << json;
        first = false;
    }
    out << "\n  }\n};\n</script>\n";

    // d3 + viz
    out << R"(<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
)";
    out << fqdup_html::JS;
    out << "\n</script>\n</body>\n</html>\n";
}

void emit_dump_members_ndjson(std::ostream& os, cf::Reader& rd) {
    rd.for_each([&os](std::uint64_t /*idx*/, const cf::ClusterRecord& r) {
        os << "{\"cluster_id\":" << r.cluster_id << ",\"members\":[";
        for (std::size_t i = 0; i < r.member_ids.size(); ++i) {
            if (i) os << ",";
            json_escape(os, r.member_ids[i]);
        }
        os << "]}\n";
    });
}

// ---- --member-of: map raw FASTQ reads → cluster_id by canonical sequence ----
//
// Walks the .fqcl, reconstructs every node sequence (parent + edits along path),
// builds a SequenceFingerprint(canonical_unmasked) → cluster_id map, then streams
// the input FASTQ and emits TSV: read_name<TAB>cluster_id (-1 = unmapped).
//
// Limitation: matches by EXACT unmasked canonical sequence. Reads merged via
// Phase-1 damage-mask hash collision (terminal-base-altered twins of a rep)
// will report -1 even though fqdup absorbed them. For interior-variant audits
// this is acceptable; a damage-aware variant would need per-position mask_pos
// preserved in .fqcl meta (currently only the count is stored).

struct MoFingerprint {
    std::uint64_t lo;
    std::uint64_t hi;
    std::uint32_t len;
    bool operator==(const MoFingerprint& o) const noexcept {
        return lo == o.lo && hi == o.hi && len == o.len;
    }
};

struct MoFingerprintHash {
    std::size_t operator()(const MoFingerprint& fp) const noexcept {
        std::uint64_t mixed = fp.lo
                            ^ (fp.hi  * 0x9e3779b97f4a7c15ULL)
                            ^ (static_cast<std::uint64_t>(fp.len) * 0x517cc1b727220a95ULL);
        return static_cast<std::size_t>(mixed ^ (mixed >> 32));
    }
};

inline void mo_revcomp_upper(const char* in, std::size_t n, char* out) {
    for (std::size_t i = 0; i < n; ++i) {
        unsigned char c = static_cast<unsigned char>(in[n - 1 - i]);
        switch (c) {
            case 'A': case 'a': out[i] = 'T'; break;
            case 'C': case 'c': out[i] = 'G'; break;
            case 'G': case 'g': out[i] = 'C'; break;
            case 'T': case 't': out[i] = 'A'; break;
            default:            out[i] = 'N'; break;
        }
    }
}

inline MoFingerprint mo_canonical_fp(const char* seq, std::size_t n,
                                     std::vector<char>& s_up,
                                     std::vector<char>& s_rc) {
    if (s_up.size() < n) s_up.resize(n);
    if (s_rc.size() < n) s_rc.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        unsigned char c = static_cast<unsigned char>(seq[i]);
        s_up[i] = static_cast<char>(std::toupper(c));
    }
    XXH128_hash_t h1 = XXH3_128bits(s_up.data(), n);
    mo_revcomp_upper(s_up.data(), n, s_rc.data());
    XXH128_hash_t h2 = XXH3_128bits(s_rc.data(), n);
    XXH128_hash_t hc;
    if (h1.high64 < h2.high64 || (h1.high64 == h2.high64 && h1.low64 <= h2.low64))
        hc = h1;
    else
        hc = h2;
    return MoFingerprint{hc.low64, hc.high64, static_cast<std::uint32_t>(n)};
}

void emit_member_of(std::ostream& os, cf::Reader& rd, const std::string& fastq_path) {
    std::unordered_map<MoFingerprint, std::uint64_t, MoFingerprintHash> fp_to_cluster;
    fp_to_cluster.reserve(rd.n_clusters() * 2);

    std::vector<char> scratch_up, scratch_rc;
    std::vector<std::vector<char>> node_seq;  // per-cluster scratch, reused

    std::uint64_t n_seqs_indexed = 0;
    std::uint64_t n_collisions   = 0;

    rd.for_each([&](std::uint64_t /*idx*/, const cf::ClusterRecord& rec) {
        if (rec.parent_seq_len == 0) return;

        // Decode parent (node 0).
        std::uint32_t max_node = 0;
        for (const auto& e : rec.edges)
            if (e.to_node > max_node) max_node = e.to_node;

        if (node_seq.size() < static_cast<std::size_t>(max_node) + 1)
            node_seq.resize(static_cast<std::size_t>(max_node) + 1);
        for (std::uint32_t n = 0; n <= max_node; ++n)
            node_seq[n].assign(rec.parent_seq_len, '\0');

        cf::decode_2bit(rec.parent_seq.data(), rec.parent_seq_len, node_seq[0].data());

        // Apply edges in storage order: each edge sets node_seq[to_node] =
        // node_seq[from_node] with one base swap at edge.pos. derep emits edges
        // in tree-construction order (child after its parent), so from_node is
        // always already populated.
        static const char A4[4] = {'A','C','G','T'};
        for (const auto& e : rec.edges) {
            std::vector<char>& src = node_seq[e.from_node];
            std::vector<char>& dst = node_seq[e.to_node];
            // Copy from parent if not yet built (first time we see this to_node).
            if (dst.size() != src.size()
                || std::memcmp(dst.data(), src.data(), src.size()) != 0)
                dst.assign(src.begin(), src.end());
            if (e.pos < dst.size())
                dst[e.pos] = A4[e.to_base & 0x3u];
        }

        // Hash every node and insert.
        for (std::uint32_t n = 0; n <= max_node; ++n) {
            const auto& s = node_seq[n];
            MoFingerprint fp = mo_canonical_fp(s.data(), s.size(), scratch_up, scratch_rc);
            auto ins = fp_to_cluster.emplace(fp, rec.cluster_id);
            if (ins.second) {
                ++n_seqs_indexed;
            } else if (ins.first->second != rec.cluster_id) {
                ++n_collisions;
            }
        }
    });

    std::fprintf(stderr,
        "[member-of] indexed %lu canonical sequences from %lu clusters (cross-cluster collisions: %lu)\n",
        static_cast<unsigned long>(n_seqs_indexed),
        static_cast<unsigned long>(rd.n_clusters()),
        static_cast<unsigned long>(n_collisions));

    auto reader = make_fastq_reader(fastq_path);
    FastqRecord fr;
    std::uint64_t n_in = 0, n_hit = 0;
    while (reader->read(fr)) {
        ++n_in;
        MoFingerprint fp = mo_canonical_fp(fr.seq.data(), fr.seq.size(),
                                           scratch_up, scratch_rc);
        auto it = fp_to_cluster.find(fp);
        std::int64_t cid = (it == fp_to_cluster.end())
                         ? -1
                         : static_cast<std::int64_t>(it->second);
        if (cid >= 0) ++n_hit;

        // Strip leading '@' and trim at first whitespace for the read name.
        const std::string& h = fr.header;
        std::size_t start = (!h.empty() && (h[0] == '@' || h[0] == '>')) ? 1 : 0;
        std::size_t sp    = h.find_first_of(" \t", start);
        std::string name  = (sp == std::string::npos)
                          ? h.substr(start)
                          : h.substr(start, sp - start);
        os << name << '\t' << cid << '\n';

        if ((n_in % 5000000) == 0)
            std::fprintf(stderr,
                "[member-of] %lu reads scanned, %lu mapped (%.2f%%)\n",
                static_cast<unsigned long>(n_in),
                static_cast<unsigned long>(n_hit),
                100.0 * static_cast<double>(n_hit) / static_cast<double>(n_in));
    }
    std::fprintf(stderr,
        "[member-of] DONE: %lu reads scanned, %lu mapped (%.2f%%)\n",
        static_cast<unsigned long>(n_in),
        static_cast<unsigned long>(n_hit),
        n_in ? 100.0 * static_cast<double>(n_hit) / static_cast<double>(n_in) : 0.0);
}

void usage(const char* prog) {
    std::cerr <<
      "Usage: " << prog << " view <file.fqcl> [options]\n"
      "  (no opts)                summary header\n"
      "  --cluster N              ASCII tree of cluster N\n"
      "  --staircase N            per-node mismatch grid for cluster N\n"
      "  --bundle [--end-k K]     group clusters by start+end k-mer (default K=16)\n"
      "  --bundle-staircase HEX   render staircase across one bundle\n"
      "  --min-bundle-size N      skip bundles smaller than N (default 2)\n"
      "  --top-members N          top N clusters by member count\n"
      "  --top-edges N            top N clusters by edge count\n"
      "  --dump-members           emit TSV: cluster_id\\tmember_id (all clusters)\n"
      "  --member-of FASTQ        emit TSV: read_name\\tcluster_id for every\n"
      "                           record in FASTQ; cluster_id=-1 if unmapped\n"
      "                           (matches by EXACT unmasked canonical sequence;\n"
      "                            damage-mask-merged twins report -1)\n"
      "  --json                   emit structured JSON (schema fqdup.view.v1)\n"
      "                           NOTE: --dump-members --json emits NDJSON,\n"
      "  --html PATH              write self-contained HTML viz (top 50 clusters)\n"
      "                           combine with --top-members N to change count\n"
      "                                 one cluster object per line (no envelope).\n";
}

} // namespace

int view_main(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); return 1; }

    std::string path;
    long cluster_n = -1;
    long staircase_n = -1;
    bool list_bundles = false;
    std::uint64_t bundle_hex = 0;
    bool bundle_staircase = false;
    int end_k = 16;
    std::size_t min_bundle = 2;
    long top_members = -1;
    long top_edges   = -1;
    bool dump_members = false;
    bool emit_json = false;
    std::string member_of_fastq;
    std::string html_out;

    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (a == "--cluster" && i + 1 < argc)        cluster_n   = std::stol(argv[++i]);
        else if (a == "--staircase" && i + 1 < argc)      staircase_n = std::stol(argv[++i]);
        else if (a == "--bundle")                          list_bundles = true;
        else if (a == "--bundle-staircase" && i + 1 < argc) {
            bundle_staircase = true;
            bundle_hex = std::stoull(argv[++i], nullptr, 16);
        }
        else if (a == "--end-k" && i + 1 < argc)         end_k = std::stoi(argv[++i]);
        else if (a == "--min-bundle-size" && i + 1 < argc) min_bundle = std::stoul(argv[++i]);
        else if (a == "--top-members" && i + 1 < argc)    top_members = std::stol(argv[++i]);
        else if (a == "--top-edges" && i + 1 < argc)      top_edges   = std::stol(argv[++i]);
        else if (a == "--dump-members")                    dump_members = true;
        else if (a == "--member-of" && i + 1 < argc)      member_of_fastq = argv[++i];
        else if (a == "--json")                            emit_json = true;
        else if (a == "--html" && i + 1 < argc)          html_out = argv[++i];
        else if (path.empty())                            path = a;
        else { std::cerr << "unknown arg: " << a << "\n"; usage(argv[0]); return 1; }
    }
    if (path.empty()) { usage(argv[0]); return 1; }

    try {
        cf::Reader rd(path);
        if (!html_out.empty()) {
            bool by_edges_html = (top_edges > 0 && top_members <= 0);
            std::size_t n = by_edges_html ? static_cast<std::size_t>(top_edges)
                          : (top_members > 0 ? static_cast<std::size_t>(top_members) : 50);
            emit_html_static(rd, path, html_out, n, by_edges_html);
            std::cerr << "wrote " << html_out << "\n";
            return 0;
        }
        if (cluster_n >= 0) {
            cf::ClusterRecord r;
            rd.read_cluster(static_cast<std::uint64_t>(cluster_n), r);
            if (emit_json) emit_cluster_json(std::cout, rd, path,
                                             static_cast<std::uint64_t>(cluster_n), r);
            else           print_tree(r);
        } else if (staircase_n >= 0) {
            cf::ClusterRecord r;
            rd.read_cluster(static_cast<std::uint64_t>(staircase_n), r);
            if (emit_json) emit_staircase_json(std::cout, rd, path,
                                               static_cast<std::uint64_t>(staircase_n), r);
            else           print_staircase(r);
        } else if (dump_members) {
            if (emit_json) {
                emit_dump_members_ndjson(std::cout, rd);
            } else {
                rd.for_each([](std::uint64_t /*idx*/, const cf::ClusterRecord& r) {
                    for (const auto& mid : r.member_ids)
                        std::printf("%lu\t%s\n",
                                    static_cast<unsigned long>(r.cluster_id),
                                    mid.c_str());
                });
            }
        } else if (!member_of_fastq.empty()) {
            emit_member_of(std::cout, rd, member_of_fastq);
        } else if (top_members > 0) {
            if (emit_json) emit_top_json(std::cout, rd, path,
                                         static_cast<std::size_t>(top_members), false);
            else           print_top(rd, static_cast<std::size_t>(top_members), false);
        } else if (top_edges > 0) {
            if (emit_json) emit_top_json(std::cout, rd, path,
                                         static_cast<std::size_t>(top_edges), true);
            else           print_top(rd, static_cast<std::size_t>(top_edges), true);
        } else if (bundle_staircase) {
            if (emit_json) emit_bundle_staircase_json(std::cout, rd, path, bundle_hex, end_k);
            else           print_bundle_staircase(rd, bundle_hex, end_k);
        } else if (list_bundles) {
            if (emit_json) emit_bundles_json(std::cout, rd, path, end_k, min_bundle);
            else           print_bundles(rd, end_k, min_bundle);
        } else {
            if (emit_json) emit_header_json(std::cout, rd, path);
            else           print_header(rd);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
