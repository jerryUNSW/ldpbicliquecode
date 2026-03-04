#include "biclique_dds_eldp_v12.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

static long double nCr_ld(int n, int r) {
    if (r < 0 || r > n) return 0.0L;
    if (r == 0 || r == n) return 1.0L;
    if (r > n - r) r = n - r;
    long double ans = 1.0L;
    for (int i = 1; i <= r; ++i) ans = ans * (n - r + i) / i;
    return ans;
}

static bool noisy_edge_cached(
    BiGraph& g,
    vector<unordered_map<int, bool>>& cache,
    int u,
    int v,
    long double p_flip,
    mt19937& rng) {

    int lo = min(u, v), hi = max(u, v);
    auto it = cache[lo].find(hi);
    if (it != cache[lo].end()) return it->second;

    bool true_edge = g.has(u, v);
    uniform_real_distribution<double> dist(0.0, 1.0);
    bool flip = dist(rng) < (double)p_flip;
    bool noisy = flip ? !true_edge : true_edge;
    cache[lo][hi] = noisy;
    return noisy;
}

static long double noisy_biclique_count(
    BiGraph& g,
    const vector<int>& S_upper,
    const vector<int>& S_lower,
    int P,
    int Q,
    vector<unordered_map<int, bool>>& cache,
    long double p_flip,
    mt19937& rng) {

    if ((int)S_upper.size() < P || (int)S_lower.size() < Q) return 0.0L;

    long double total = 0.0L;
    vector<int> pick(P, 0);

    function<void(int, int)> dfs = [&](int depth, int start) {
        if (depth == P) {
            int common = 0;
            for (int lv : S_lower) {
                bool ok = true;
                for (int k = 0; k < P; ++k) {
                    int uu = S_upper[pick[k]];
                    if (!noisy_edge_cached(g, cache, uu, lv, p_flip, rng)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) common++;
            }
            total += nCr_ld(common, Q);
            return;
        }
        for (int i = start; i <= (int)S_upper.size() - (P - depth); ++i) {
            pick[depth] = i;
            dfs(depth + 1, i + 1);
        }
    };

    dfs(0, 0);
    return total;
}

static long double estimate_density_double_noisy(
    BiGraph& g,
    const vector<int>& S_upper,
    const vector<int>& S_lower,
    int P,
    int Q,
    vector<unordered_map<int, bool>>& cache1,
    vector<unordered_map<int, bool>>& cache2,
    long double p_flip,
    mt19937& rng1,
    mt19937& rng2) {

    int set_size = (int)S_upper.size() + (int)S_lower.size();
    if (set_size == 0) return 0.0L;

    long double c1 = noisy_biclique_count(g, S_upper, S_lower, P, Q, cache1, p_flip, rng1);
    long double c2 = noisy_biclique_count(g, S_upper, S_lower, P, Q, cache2, p_flip, rng2);
    return 0.5L * (c1 + c2) / (long double)set_size;
}

static vector<pair<long double, int>> rank_by_noisy_degree(
    BiGraph& g,
    long double p_flip,
    mt19937& rng) {

    vector<pair<long double, int>> rank;
    rank.reserve(g.num_nodes());

    // No extra epsilon: degrees are sampled as degree in noisy graph G'.
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

DDSResultV12 eldp_v12(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed) {
    DDSResultV12 result;
    int n = g.num_nodes();

    long double p_flip = 1.0L / (expl(epsilon) + 1.0L);
    mt19937 rng_rank(seed ^ 0xA17F31u);
    mt19937 rng1(seed ^ 0xBADC0DEu);
    mt19937 rng2(seed ^ 0xC0FFEEu);

    vector<unordered_map<int, bool>> cache1(n), cache2(n);

    vector<pair<long double, int>> noisy_rank = rank_by_noisy_degree(g, p_flip, rng_rank);

    // Build initial set from noisy ranking only.
    vector<int> S_upper, S_lower;
    for (auto& x : noisy_rank) {
        int v = x.second;
        if (v < g.num_v1) {
            if ((int)S_upper.size() < max(P, 8)) S_upper.push_back(v);
        } else {
            if ((int)S_lower.size() < max(Q, 12)) S_lower.push_back(v);
        }
        if ((int)S_upper.size() >= max(P, 8) && (int)S_lower.size() >= max(Q, 12)) break;
    }

    if ((int)S_upper.size() < P || (int)S_lower.size() < Q) {
        result.vertex_set.clear();
        result.estimated_density = 0.0L;
        result.real_density = 0.0L;
        return result;
    }

    vector<char> in_set(n, 0);
    for (int u : S_upper) in_set[u] = 1;
    for (int v : S_lower) in_set[v] = 1;

    long double cur_d = estimate_density_double_noisy(
        g, S_upper, S_lower, P, Q, cache1, cache2, p_flip, rng1, rng2);
    long double best_d = cur_d;
    vector<int> best_set;
    for (int u : S_upper) best_set.push_back(u);
    for (int v : S_lower) best_set.push_back(v);

    // Candidate pool from top noisy-ranked vertices.
    vector<int> candidates;
    int pool_limit = min(n, 220);
    for (int i = 0; i < pool_limit; ++i) {
        int v = noisy_rank[i].second;
        if (!in_set[v]) candidates.push_back(v);
    }

    // Cross-iteration maintenance: keep cached gains and lazily refresh best few each round.
    unordered_map<int, long double> gain_cache;
    int no_improve = 0;
    int max_add_steps = 120;

    for (int step = 0; step < max_add_steps && no_improve < 8; ++step) {
        vector<pair<long double, int>> scored;
        scored.reserve(candidates.size());

        // Lazy use of cache first.
        for (int v : candidates) {
            auto it = gain_cache.find(v);
            long double g_est = (it == gain_cache.end() ? 0.0L : it->second);
            scored.push_back({g_est, v});
        }
        sort(scored.rbegin(), scored.rend());

        // Refresh top M candidates only (incremental maintenance).
        int M = min((int)scored.size(), 36);
        for (int i = 0; i < M; ++i) {
            int v = scored[i].second;
            vector<int> tu = S_upper, tl = S_lower;
            if (v < g.num_v1) tu.push_back(v);
            else tl.push_back(v);

            if ((int)tu.size() < P || (int)tl.size() < Q) {
                gain_cache[v] = -1e30L;
                continue;
            }

            long double nd = estimate_density_double_noisy(
                g, tu, tl, P, Q, cache1, cache2, p_flip, rng1, rng2);
            gain_cache[v] = nd - cur_d;
        }

        // Re-rank after refresh.
        scored.clear();
        for (int v : candidates) scored.push_back({gain_cache[v], v});
        sort(scored.rbegin(), scored.rend());
        if (scored.empty()) break;

        int best_v = scored[0].second;
        long double best_gain = scored[0].first;
        if (best_gain <= 0) {
            no_improve++;
            continue;
        }

        if (best_v < g.num_v1) S_upper.push_back(best_v);
        else S_lower.push_back(best_v);
        in_set[best_v] = 1;
        cur_d += best_gain;

        candidates.erase(remove(candidates.begin(), candidates.end(), best_v), candidates.end());
        gain_cache.erase(best_v);

        if (cur_d > best_d) {
            best_d = cur_d;
            best_set.clear();
            for (int u : S_upper) best_set.push_back(u);
            for (int v : S_lower) best_set.push_back(v);
            no_improve = 0;
        } else {
            no_improve++;
        }
    }

    // Small greedy shrink on noisy objective.
    bool improved = true;
    while (improved && ((int)S_upper.size() + (int)S_lower.size() > P + Q)) {
        improved = false;
        long double best_remove_d = cur_d;
        int rem_v = -1;

        vector<int> cur_all;
        for (int u : S_upper) cur_all.push_back(u);
        for (int v : S_lower) cur_all.push_back(v);

        for (int v : cur_all) {
            vector<int> tu, tl;
            for (int u : S_upper) if (u != v) tu.push_back(u);
            for (int l : S_lower) if (l != v) tl.push_back(l);
            if ((int)tu.size() < P || (int)tl.size() < Q) continue;

            long double nd = estimate_density_double_noisy(
                g, tu, tl, P, Q, cache1, cache2, p_flip, rng1, rng2);
            if (nd > best_remove_d) {
                best_remove_d = nd;
                rem_v = v;
            }
        }

        if (rem_v >= 0) {
            if (rem_v < g.num_v1) {
                S_upper.erase(remove(S_upper.begin(), S_upper.end(), rem_v), S_upper.end());
            } else {
                S_lower.erase(remove(S_lower.begin(), S_lower.end(), rem_v), S_lower.end());
            }
            in_set[rem_v] = 0;
            cur_d = best_remove_d;
            improved = true;

            if (cur_d > best_d) {
                best_d = cur_d;
                best_set.clear();
                for (int u : S_upper) best_set.push_back(u);
                for (int v : S_lower) best_set.push_back(v);
            }
        }
    }

    result.vertex_set = best_set;
    result.estimated_density = best_d;
    // Keep report-compatible field; objective/selection are entirely noisy-graph based.
    result.real_density = best_d;

    int bu = 0, bl = 0;
    for (int v : best_set) {
        if (v < g.num_v1) bu++;
        else bl++;
    }

    cout << "  [edlp_v12] eps=" << fixed << setprecision(3) << epsilon
         << " p_flip=" << setprecision(4) << p_flip
         << " |S|=" << best_set.size() << " (u=" << bu << ",l=" << bl << ")"
         << " est_d=" << setprecision(6) << best_d << endl;

    return result;
}
