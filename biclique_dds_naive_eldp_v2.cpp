#include "biclique_dds_naive_eldp_v2.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

// Reuse existing unbiased biclique estimator and true counter.
extern long double estimate_biclique_count_dp(
    const vector<int>& S_upper,
    const vector<int>& S_lower,
    int P,
    int Q,
    BiGraph& g_noisy,
    BiGraph& g_orig,
    long double p_flip,
    long double eps1,
    mt19937& rng);

extern long double count_bicliques_in_set(BiGraph& g, const vector<int>& S, int P, int Q);

static vector<pair<long double, int>> rank_by_noisy_degree(
    BiGraph& g,
    long double p_flip,
    mt19937& rng) {

    vector<pair<long double, int>> rank;
    rank.reserve(g.num_nodes());

    for (int u = 0; u < g.num_v1; ++u) {
        int true_deg = (int)g.neighbor[u].size();
        int n_opp = g.num_v2;
        binomial_distribution<int> keep_true(true_deg, (double)(1.0L - p_flip));
        binomial_distribution<int> add_false(max(0, n_opp - true_deg), (double)p_flip);
        int noisy_deg = keep_true(rng) + add_false(rng);
        rank.push_back({(long double)noisy_deg, u});
    }

    for (int lv = 0; lv < g.num_v2; ++lv) {
        int v = g.num_v1 + lv;
        int true_deg = (int)g.neighbor[v].size();
        int n_opp = g.num_v1;
        binomial_distribution<int> keep_true(true_deg, (double)(1.0L - p_flip));
        binomial_distribution<int> add_false(max(0, n_opp - true_deg), (double)p_flip);
        int noisy_deg = keep_true(rng) + add_false(rng);
        rank.push_back({(long double)noisy_deg, v});
    }

    sort(rank.rbegin(), rank.rend());
    return rank;
}

static long double estimate_density_on_set(
    BiGraph& g,
    const vector<int>& S,
    int P,
    int Q,
    BiGraph& g_noisy,
    long double p_flip,
    long double eps1,
    mt19937& rng) {

    vector<int> su, sl;
    su.reserve(S.size());
    sl.reserve(S.size());
    for (int v : S) {
        if (g.is_upper(v)) su.push_back(v);
        else sl.push_back(v);
    }

    if ((int)su.size() < P || (int)sl.size() < Q || S.empty()) return 0.0L;

    long double bc_est = estimate_biclique_count_dp(su, sl, P, Q, g_noisy, g, p_flip, eps1, rng);
    if (!isfinite((double)bc_est) || bc_est < 0) bc_est = 0.0L;
    return bc_est / (long double)S.size();
}

NaiveV2Result naive_eldp_v2(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    NaiveV2Result result;

    int n = g.num_nodes();
    mt19937 rng(seed);

    // No extra budget split: build ranking from noisy graph degree only.
    long double p_flip = 1.0L / (expl(epsilon) + 1.0L);
    long double eps1 = epsilon;

    BiGraph g_noisy(g);
    g_noisy.edge_vector.resize(n);
    for (int i = 0; i < n; ++i) g_noisy.edge_vector[i].clear();

    vector<pair<long double, int>> noisy_rank = rank_by_noisy_degree(g, p_flip, rng);

    // More aggressive init: larger K and candidate sizes (fix |S|=0 on big graphs)
    int K = min(n, 500);
    vector<int> candidate_sizes;
    candidate_sizes.push_back(P + Q);
    for (int s : {50, 100, 150, 200, 250, 300, 400}) {
        if (s >= P + Q && s <= K) candidate_sizes.push_back(s);
    }
    sort(candidate_sizes.begin(), candidate_sizes.end());
    candidate_sizes.erase(unique(candidate_sizes.begin(), candidate_sizes.end()), candidate_sizes.end());

    long double best_est_density = -1.0L;
    vector<int> best_set;

    for (int keep : candidate_sizes) {
        vector<int> cur;
        cur.reserve(keep);
        int up = 0, low = 0;
        for (int i = 0; i < keep; ++i) {
            int v = noisy_rank[i].second;
            cur.push_back(v);
            if (g.is_upper(v)) up++;
            else low++;
        }
        if (up < P || low < Q) continue;

        long double est_d = estimate_density_on_set(g, cur, P, Q, g_noisy, p_flip, eps1, rng);
        if (est_d > best_est_density) {
            best_est_density = est_d;
            best_set = cur;
        }
    }

    if (best_set.empty()) {
        result.vertex_set.clear();
        result.estimated_density = 0.0L;
        result.real_density = 0.0L;
        cout << "  [naive_eldp_v2] Final |S|=0 est_d=0.000000 p_flip="
             << fixed << setprecision(4) << p_flip << endl;
        return result;
    }

    // IMPORTANT: not used in optimization, only for final evaluation/reporting.
    long double true_bc = count_bicliques_in_set(g, best_set, P, Q);
    long double real_density = true_bc / (long double)best_set.size();

    result.vertex_set = best_set;
    result.estimated_density = max((long double)0.0L, best_est_density);
    result.real_density = max((long double)0.0L, real_density);

    cout << "  [naive_eldp_v2] Final |S|=" << best_set.size()
         << " est_d=" << fixed << setprecision(6) << result.estimated_density
         << " real_d=" << setprecision(6) << result.real_density
         << " p_flip=" << setprecision(4) << p_flip << endl;

    return result;
}
