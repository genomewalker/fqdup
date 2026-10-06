#pragma once
#include <array>
#include <algorithm>
namespace fqalign { namespace full_complexity {
static constexpr int MAX_READ_LEN=512,DUST_W=64;
static constexpr float DUST_MAX=.10f;
static bool overlap_is_low_complexity(const char* seq, int len) {
    if (len < 5) return false;
    auto b2 = [](char c) -> int { switch (c) { case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3; } return -1; };
    std::array<int, MAX_READ_LEN> tri{};
    const int nt = len - 2;
    for (int i = 0; i < nt; ++i) {
        const int a = b2(seq[i]), b = b2(seq[i + 1]), c = b2(seq[i + 2]);
        tri[i] = (a < 0 || b < 0 || c < 0) ? -1 : (a * 16 + b * 4 + c);
    }
    const int wt = std::min(DUST_W - 2, nt);
    if (wt < 2) return false;
    int cnt[64] = {}, valid = 0;
    long score = 0;
    auto add = [&](int t) { if (t >= 0) { score += cnt[t]; ++cnt[t]; ++valid; } };
    auto del = [&](int t) { if (t >= 0) { --cnt[t]; score -= cnt[t]; --valid; } };
    auto exceeds = [&]() {
        if (valid < 2) return false;
        const float maximum = static_cast<float>(valid) * (valid - 1) / 2.f;
        return static_cast<float>(score) / maximum > DUST_MAX;
    };
    for (int i = 0; i < wt; ++i) add(tri[i]);
    if (exceeds()) return true;
    for (int i = wt; i < nt; ++i) {
        del(tri[i - wt]); add(tri[i]);
        if (exceeds()) return true;
    }
    return false;
}

}}

