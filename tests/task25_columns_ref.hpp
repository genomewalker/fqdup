#pragma once
// Frozen from 9fe83c1; namespace only changed.
namespace frozen_columns {
struct Column { int n = 0; std::vector<std::pair<char, std::vector<int>>> votes; int top = -1; }; // votes in first-seen order
std::vector<Column> columns(const std::vector<Pair>& P, const std::vector<int>& g, int side) {
    size_t L = 0; for (int i : g) L = std::max(L, P[i].r[side].s.size());
    std::vector<Column> cols(L);
    for (size_t c = 0; c < L; ++c) { auto& col = cols[c];
        for (int i : g) { const auto& s = P[i].r[side].s; if (c >= s.size() || s[c] == 'N') continue; ++col.n; char b = s[c];
            auto it = std::find_if(col.votes.begin(), col.votes.end(), [&](auto& v) { return v.first == b; });
            if (it == col.votes.end()) col.votes.push_back({b, {i}}); else it->second.push_back(i); }
        int top = -1; for (size_t v = 0; v < col.votes.size(); ++v) if (top < 0 || col.votes[v].second.size() > col.votes[top].second.size()) top = int(v);
        bool tie = false; if (top >= 0) for (size_t v = 0; v < col.votes.size(); ++v) if (int(v) != top && col.votes[v].second.size() == col.votes[top].second.size()) tie = true;
        col.top = tie ? -1 : top; }
    return cols;
}
}
