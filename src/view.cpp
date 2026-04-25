// fqdup view — inspect .fqcl cluster genealogy files.
// Modes: header summary (default), single-cluster tree/staircase, sibling bundling.

#include "fqdup/cluster_format.hpp"
#include "fqdup/bundle_key.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cf = fqdup::clusterfmt;

namespace {

using fqdup::bundlekey::from_cluster;

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
      "  --dump-members           emit TSV: cluster_id\\tmember_id (all clusters)\n";
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
        else if (path.empty())                            path = a;
        else { std::cerr << "unknown arg: " << a << "\n"; usage(argv[0]); return 1; }
    }
    if (path.empty()) { usage(argv[0]); return 1; }

    try {
        cf::Reader rd(path);
        if (cluster_n >= 0) {
            cf::ClusterRecord r;
            rd.read_cluster(static_cast<std::uint64_t>(cluster_n), r);
            print_tree(r);
        } else if (staircase_n >= 0) {
            cf::ClusterRecord r;
            rd.read_cluster(static_cast<std::uint64_t>(staircase_n), r);
            print_staircase(r);
        } else if (dump_members) {
            rd.for_each([](std::uint64_t /*idx*/, const cf::ClusterRecord& r) {
                for (const auto& mid : r.member_ids)
                    std::printf("%lu\t%s\n",
                                static_cast<unsigned long>(r.cluster_id),
                                mid.c_str());
            });
        } else if (top_members > 0) {
            print_top(rd, static_cast<std::size_t>(top_members), /*by_edges=*/false);
        } else if (top_edges > 0) {
            print_top(rd, static_cast<std::size_t>(top_edges), /*by_edges=*/true);
        } else if (bundle_staircase) {
            print_bundle_staircase(rd, bundle_hex, end_k);
        } else if (list_bundles) {
            print_bundles(rd, end_k, min_bundle);
        } else {
            print_header(rd);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
