/*
 * biclique_dds_peeling_eldp.cpp
 * Peeling-ELDP: Port of pqBicliqeDensest peeling to noisy graph G'.
 * All operations use G' only. ELDP compliant.
 */

#include "biclique_dds_peeling_eldp.h"
#include "biclique_dds_eldp.h"
#include "biclique.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <unordered_map>
#include <vector>

using namespace std;

// Rank vertices by noisy degree (from G' only)
static vector<pair<int, int>> rank_by_noisy_degree_g(BiGraph& g, BiGraph& g_noisy,
    long double p_flip, mt19937& rng) {
    int n = g.num_nodes();
    vector<pair<int, int>> rank;
    rank.reserve(n);

    for (int u = 0; u < g.num_v1; u++) {
        int cnt = 0;
        for (int v = g.num_v1; v < n; v++) {
            if (query_noisy_edge(g, g_noisy, u, v, p_flip, rng)) cnt++;
        }
        rank.push_back({cnt, u});
    }
    for (int lv = 0; lv < g.num_v2; lv++) {
        int v = g.num_v1 + lv;
        int cnt = 0;
        for (int u = 0; u < g.num_v1; u++) {
            if (query_noisy_edge(g, g_noisy, u, v, p_flip, rng)) cnt++;
        }
        rank.push_back({cnt, v});
    }
    sort(rank.rbegin(), rank.rend());
    return rank;
}

vector<int> peeling_eldp_core(BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1, mt19937& rng) {
    int n = g.num_nodes();
    auto deg_sorted = rank_by_noisy_degree_g(g, g_noisy, p_flip, rng);
    int K = min(n, 500);

    vector<int> S;
    vector<int> S_upper, S_lower;
    int uc = 0, lc = 0;
    for (int i = 0; i < K && i < (int)deg_sorted.size(); i++) {
        int v = deg_sorted[i].second;
        S.push_back(v);
        if (g.is_upper(v)) { S_upper.push_back(v); uc++; }
        else { S_lower.push_back(v); lc++; }
    }

    if (uc < P || lc < Q) {
        return {};
    }

    long double total_est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    if (!isfinite((double)total_est) || total_est < 0) total_est = 0.0L;

    long double best_density = total_est / (long double)S.size();
    vector<int> best_set = S;
    int mostLoops = (int)S.size();

    vector<bool> vis(n, false);
    for (int v : S) vis[v] = true;

    for (int i = 0; i < mostLoops - P - Q; i++) {
        int best_rem = -1;
        long double min_contrib = 1e30L;

        for (int v : S) {
            vector<int> tu, tl;
            for (int u : S_upper) if (u != v) tu.push_back(u);
            for (int u : S_lower) if (u != v) tl.push_back(u);
            if ((int)tu.size() < P || (int)tl.size() < Q) continue;

            long double contrib = estimate_biclique_gain(tu, tl, v, g.is_upper(v), P, Q, g_noisy, g, p_flip, eps1, rng);
            if (contrib < min_contrib) {
                min_contrib = contrib;
                best_rem = v;
            }
        }

        if (best_rem < 0) break;

        total_est -= min_contrib;
        if (total_est < 0) total_est = 0.0L;

        S.erase(remove(S.begin(), S.end(), best_rem), S.end());
        if (g.is_upper(best_rem)) {
            S_upper.erase(remove(S_upper.begin(), S_upper.end(), best_rem), S_upper.end());
        } else {
            S_lower.erase(remove(S_lower.begin(), S_lower.end(), best_rem), S_lower.end());
        }
        vis[best_rem] = false;

        long double cur_density = (S.size() > 0) ? total_est / (long double)S.size() : 0.0L;
        if (cur_density > best_density) {
            best_density = cur_density;
            best_set = S;
        }

        if (total_est <= 0 || total_est < (P + Q) * best_density) break;
    }

    return best_set;
}

DDSResult peeling_eldp(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    DDSResult result;
    mt19937 rng(seed);
    int n = g.num_nodes();

    long double p_flip = 1.0L / (expl(epsilon) + 1.0L);
    long double eps1 = epsilon;

    BiGraph g_noisy(g);
    g_noisy.edge_vector.resize(n);
    for (int i = 0; i < n; i++) g_noisy.edge_vector[i].clear();

    vector<int> best_set = peeling_eldp_core(g, g_noisy, P, Q, p_flip, eps1, rng);

    result.vertex_set = best_set;
    if (best_set.empty()) {
        result.estimated_density = 0.0L;
        result.real_density = 0.0L;
        cout << "  [peeling_eldp] Final |S|=0" << endl;
        return result;
    }

    vector<int> su, sl;
    for (int v : best_set) {
        if (g.is_upper(v)) su.push_back(v);
        else sl.push_back(v);
    }
    result.estimated_density = estimate_biclique_count_dp(su, sl, P, Q, g_noisy, g, p_flip, eps1, rng)
        / (long double)best_set.size();
    result.real_density = count_bicliques_in_set(g, best_set, P, Q) / (long double)best_set.size();

    cout << "  [peeling_eldp] Final |S|=" << best_set.size()
         << " est_d=" << fixed << setprecision(6) << result.estimated_density
         << " real_d=" << result.real_density << endl;

    return result;
}
