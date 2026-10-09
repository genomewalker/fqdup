#pragma once
// Frozen key extraction before Step38 (746d472).
std::vector<std::string> frozen_keys_of(const std::vector<Pair>& P) {
    std::vector<std::string> keys(P.size());
    for (size_t i = 0; i < P.size(); ++i) {
        const auto &a = P[i].r[0].s, &b = P[i].r[1].s; if (int(a.size()) < KEYLEN || int(b.size()) < KEYLEN) continue;
        std::string k = a.substr(0, KEYLEN) + b.substr(0, KEYLEN);
        if (std::count(k.begin(), k.end(), 'N') > MAXN) continue;
        if (!g_rev8) {
            if (low_complexity(a) || low_complexity(b)) continue;
            bool adapter = false; for (int s = 0; s < 2; ++s) { auto at = P[i].r[s].s.find(ADAPTERS[s].substr(0, ADSEED)); if (at != std::string::npos && int(at) < ADZONE) adapter = true; }
            if (!adapter) keys[i] = k;
        } else { // Revision 8: complexity on the biological prefix (before adapter seed or terminal G>=10 run); adapter must start at >= KEYLEN
            bool ok = true;
            for (int s = 0; s < 2 && ok; ++s) { const auto& r = P[i].r[s].s; auto at = r.find(ADAPTERS[s].substr(0, ADSEED)); if (at != std::string::npos && int(at) < KEYLEN) { ok = false; break; }
                size_t e = r.size(); if (at != std::string::npos) e = std::min(e, at); size_t g = r.size(); while (g > 0 && r[g - 1] == 'G') --g; if (r.size() - g >= 10) e = std::min(e, g);
                if (low_complexity(r.substr(0, e))) ok = false; }
            if (ok) keys[i] = k;
        }
    }
    return keys;
}
