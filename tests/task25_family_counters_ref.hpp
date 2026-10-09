#pragma once
// Frozen from 8f65564; only namespace and the inline hiq wrapper changed.
#include "task25_columns_ref.hpp"
namespace frozen_family_counters {
using frozen_columns::columns;
using frozen_columns::Column;
bool second_molecule(const std::vector<Pair>& P, const std::vector<int>& g) {
    const std::array<std::vector<Column>,2> all_columns{columns(P,g,0),columns(P,g,1)};
    std::map<std::pair<int, int>, int> shared; long multi = 0, cover = 0;
    for (int side = 0; side < 2; ++side) for (const auto& col : all_columns[side]) if (col.n >= 2) { ++cover; if (col.votes.size() > 1) ++multi; }
    if (cover && double(multi) / cover > MAXDIS) return true;
    for (int side = 0; side < 2; ++side) for (const auto& col : all_columns[side]) { if (col.top < 0) continue;
        for (size_t v = 0; v < col.votes.size(); ++v) { if (int(v) == col.top || col.votes[v].second.size() < 2) continue; const auto& ids = col.votes[v].second;
            for (size_t x = 0; x < ids.size(); ++x) for (size_t y = x + 1; y < ids.size(); ++y) ++shared[{ids[x], ids[y]}]; } }
    for (auto& [k, v] : shared) if (v >= 2) return true;
    return false;
}
std::array<int,MAXCOMP> high_quality_counts(const std::vector<Pair>& P,const std::vector<int>& g0) {
            std::map<int, int> hiq;
            for (int side = 0; side < 2; ++side) { auto cols = columns(P, g0, side);
                for (size_t c = 0; c < cols.size(); ++c) { auto& col = cols[c]; if (col.top < 0) continue; int topq = -1; for (int i : col.votes[col.top].second) topq = std::max(topq, qv(P[i].r[side], c));
                    if (topq < 30) continue; for (size_t v = 0; v < col.votes.size(); ++v) if (int(v) != col.top) for (int i : col.votes[v].second) if (qv(P[i].r[side], c) >= 30) ++hiq[i]; } }

    std::array<int,MAXCOMP> out{};for(std::size_t i=0;i<g0.size();++i)out[i]=hiq[g0[i]];return out;
}
}
