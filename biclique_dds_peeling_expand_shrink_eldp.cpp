/*
 * biclique_dds_peeling_expand_shrink_eldp.cpp
 * Peeling + Expand + Shrink on noisy graph G'. ELDP compliant.
 */

#include "biclique_dds_peeling_expand_shrink_eldp.h"
#include "biclique_dds_peeling_eldp.h"
#include "biclique_dds_eldp.h"
#include "biclique.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

using namespace std;

static pair<vector<int>, long double> greedy_add_eldp(
    BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1,
    const vector<int>& seed_set, int max_size, mt19937& rng) {

    int n = g.num_nodes();
    vector<int> S = seed_set;
    vector<bool> in_S(n, false);
    vector<int> S_upper, S_lower;

    for (int v : S) {
        in_S[v] = true;
        if (g.is_upper(v)) S_upper.push_back(v);
        else S_lower.push_back(v);
    }

    if ((int)S_upper.size() < P || (int)S_lower.size() < Q)
        return {S, 0.0L};

    long double current_est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double best_density = current_est / (long double)S.size();
    vector<int> best_set = S;
    int no_improve = 0;

    while ((int)S.size() < max_size && no_improve < 5) {
        unordered_set<int> candidate_set;
        for (int v : S) {
            for (auto& [key, val] : g_noisy.edge_vector[v]) {
                if (val && !in_S[key]) candidate_set.insert(key);
            }
        }
        if ((int)candidate_set.size() < 20) {
            for (int v : S) {
                for (int u = 0; u < n; u++) {
                    if (in_S[u]) continue;
                    if ((g.is_upper(v) && !g.is_upper(u)) || (!g.is_upper(v) && g.is_upper(u))) {
                        if (query_noisy_edge(g, g_noisy, v, u, p_flip, rng))
                            candidate_set.insert(u);
                    }
                }
            }
        }

        vector<pair<int, int>> candidates_scored;
        for (int v : candidate_set) {
            bool is_up = g.is_upper(v);
            const vector<int>& opp = is_up ? S_lower : S_upper;
            if (opp.empty()) continue;
            int links = 0;
            for (int u : opp) {
                if (query_noisy_edge(g, g_noisy, v, u, p_flip, rng)) links++;
            }
            candidates_scored.push_back({links, v});
        }
        sort(candidates_scored.rbegin(), candidates_scored.rend());

        int K_eval = min((int)candidates_scored.size(), 15);
        if (K_eval == 0) break;

        vector<pair<long double, int>> gain_scored;
        for (int i = 0; i < K_eval; i++) {
            int v = candidates_scored[i].second;
            bool is_up = g.is_upper(v);
            long double gain = estimate_biclique_gain(S_upper, S_lower, v, is_up, P, Q, g_noisy, g, p_flip, eps1, rng);
            gain_scored.push_back({gain, v});
        }
        sort(gain_scored.rbegin(), gain_scored.rend());

        int pick_idx = 0;
        if (gain_scored.size() > 1) {
            uniform_int_distribution<int> pick_dist(0, min((int)gain_scored.size(), 2) - 1);
            pick_idx = pick_dist(rng);
        }

        int best_v = gain_scored[pick_idx].second;
        long double gain = gain_scored[pick_idx].first;

        S.push_back(best_v);
        in_S[best_v] = true;
        if (g.is_upper(best_v)) S_upper.push_back(best_v);
        else S_lower.push_back(best_v);

        current_est += gain;
        long double cur_d = current_est / (long double)S.size();

        if (cur_d > best_density) {
            best_density = cur_d;
            best_set = S;
            no_improve = 0;
        } else {
            no_improve++;
        }
    }
    return {best_set, best_density};
}

static pair<vector<int>, long double> greedy_shrink_eldp(
    const vector<int>& initial_set,
    BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1, mt19937& rng) {

    vector<int> S = initial_set;
    vector<int> S_upper, S_lower;
    for (int u : S) {
        if (g.is_upper(u)) S_upper.push_back(u);
        else S_lower.push_back(u);
    }

    if (S.empty()) return {{}, 0.0L};
    if ((int)S_upper.size() < P || (int)S_lower.size() < Q) return {S, 0.0L};

    long double est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double cur_d = est / (long double)S.size();
    long double best_d = cur_d;
    vector<int> best_set = S;

    bool improved = true;
    while (improved && (int)S.size() > P + Q) {
        improved = false;
        long double best_rem_d = cur_d;
        int best_rem_idx = -1;

        for (int idx = 0; idx < (int)S.size(); idx++) {
            int v = S[idx];
            vector<int> tu, tl;
            for (int u : S_upper) if (u != v) tu.push_back(u);
            for (int u : S_lower) if (u != v) tl.push_back(u);
            if ((int)tu.size() < P || (int)tl.size() < Q) continue;

            long double ne = estimate_biclique_count_dp(tu, tl, P, Q, g_noisy, g, p_flip, eps1, rng);
            long double nd = ne / (long double)(S.size() - 1);
            if (nd > best_rem_d) {
                best_rem_d = nd;
                best_rem_idx = idx;
                improved = true;
            }
        }

        if (best_rem_idx >= 0) {
            S.erase(S.begin() + best_rem_idx);
            S_upper.clear();
            S_lower.clear();
            for (int u : S) {
                if (g.is_upper(u)) S_upper.push_back(u);
                else S_lower.push_back(u);
            }
            cur_d = best_rem_d;
            if (cur_d > best_d) {
                best_d = cur_d;
                best_set = S;
            }
        }
    }
    return {best_set, best_d};
}

DDSResult peeling_expand_shrink_eldp(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    DDSResult result;
    mt19937 rng(seed);
    int n = g.num_nodes();

    long double p_flip = 1.0L / (expl(epsilon) + 1.0L);
    long double eps1 = epsilon;

    BiGraph g_noisy(g);
    g_noisy.edge_vector.resize(n);
    for (int i = 0; i < n; i++) g_noisy.edge_vector[i].clear();

    vector<int> seed_set = peeling_eldp_core(g, g_noisy, P, Q, p_flip, eps1, rng);
    if (seed_set.empty()) {
        result.vertex_set.clear();
        result.estimated_density = 0.0L;
        result.real_density = 0.0L;
        cout << "  [peeling_expand_shrink_eldp] Peeling returned empty" << endl;
        return result;
    }

    auto [add_set, add_density] = greedy_add_eldp(g, g_noisy, P, Q, p_flip, eps1, seed_set, min(n, 1000), rng);
    auto [final_set, final_density] = greedy_shrink_eldp(add_set, g, g_noisy, P, Q, p_flip, eps1, rng);

    result.vertex_set = final_set;
    result.estimated_density = final_density;
    if (!final_set.empty()) {
        result.real_density = count_bicliques_in_set(g, final_set, P, Q) / (long double)final_set.size();
    } else {
        result.real_density = 0.0L;
    }

    cout << "  [peeling_expand_shrink_eldp] Final |S|=" << final_set.size()
         << " est_d=" << fixed << setprecision(6) << result.estimated_density
         << " real_d=" << result.real_density << endl;

    return result;
}
