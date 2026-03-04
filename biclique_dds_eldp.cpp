/*
 * biclique_dds_eldp.cpp — v10 Dual-Strategy Initialization + Biclique Refinement
 *
 * Two initialization strategies:
 *   A) Charikar-style edge density peeling on debiased noisy graph
 *   B) ABCore peeling on debiased noisy graph
 * Both produce a seed set, which is then refined via greedy add + shrink
 * using the unbiased DP-based (p,q)-biclique density estimator.
 */

#include "biclique_dds_eldp.h"
#include "biclique.h"
#include "sqlite3.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <numeric>
#include <cassert>
#include <fstream>
#include <filesystem>

using namespace std;

// Global trace file handle (set by main_dds.cpp if needed)
static ofstream* g_trace_file = nullptr;

void set_trace_file(ofstream* fp) {
    g_trace_file = fp;
}

static void write_trace(int set_size, long double estimated_d, long double real_d, const string& label) {
    if (g_trace_file && g_trace_file->is_open()) {
        (*g_trace_file) << set_size << "," << fixed << setprecision(10) 
                        << estimated_d << "," << real_d << "," << label << "\n";
    }
}

extern void init_genrand(unsigned long s);
extern double genrand_real2(void);

// ============================================================
// Noisy graph: lazy evaluation with caching
// ============================================================

bool query_noisy_edge(BiGraph& g, BiGraph& g_noisy, int u, int v,
                      long double p_flip, mt19937& rng) {
    int smaller = min(u, v);
    int larger = max(u, v);
    auto it = g_noisy.edge_vector[smaller].find(larger);
    if (it != g_noisy.edge_vector[smaller].end()) return it->second;

    bool true_edge = g.has(u, v);
    uniform_real_distribution<double> dist(0.0, 1.0);
    bool noisy_edge = (dist(rng) < (double)p_flip) ? !true_edge : true_edge;
    g_noisy.edge_vector[smaller][larger] = noisy_edge;
    return noisy_edge;
}

// ============================================================
// Helper: nCr for long double
// ============================================================
static long double nCr_double(int n, int r) {
    if (r < 0 || r > n) return 0.0L;
    if (r == 0 || r == n) return 1.0L;
    if (r > n / 2) r = n - r;
    long double res = 1.0L;
    for (int i = 1; i <= r; ++i) {
        res = res * (n - i + 1) / i;
    }
    return res;
}

// ============================================================
// True biclique counting (on original graph, for evaluation)
// ============================================================

// Direct biclique counter: enumerate C(S_upper, P) upper subsets,
// count C(common_lower_neighbors, Q) for each.
long double count_bicliques_in_set(BiGraph& g, const vector<int>& S, int P, int Q) {
    vector<int> su, sl;
    for (int v : S) { if (g.is_upper(v)) su.push_back(v); else sl.push_back(v); }
    if ((int)su.size() < P || (int)sl.size() < Q) return 0.0L;

    unordered_set<int> lower_set(sl.begin(), sl.end());
    int nu = su.size(), nl = sl.size();

    // Build adjacency: for each upper vertex, bitvector of which lower vertices it connects to
    vector<vector<bool>> adj(nu, vector<bool>(nl, false));
    unordered_map<int, int> lower_idx;
    for (int i = 0; i < nl; i++) lower_idx[sl[i]] = i;
    for (int i = 0; i < nu; i++) {
        for (int nb : g.neighbor[su[i]]) {
            auto it = lower_idx.find(nb);
            if (it != lower_idx.end()) adj[i][it->second] = true;
        }
    }

    // Enumerate all C(nu, P) upper subsets via recursion
    long double total = 0.0L;
    vector<int> chosen(P);

    function<void(int, int)> enumerate = [&](int pos, int start) {
        if (pos == P) {
            // Count common lower neighbors of all chosen upper vertices
            int common = 0;
            for (int j = 0; j < nl; j++) {
                bool all_connected = true;
                for (int k = 0; k < P; k++) {
                    if (!adj[chosen[k]][j]) { all_connected = false; break; }
                }
                if (all_connected) common++;
            }
            total += nCr_double(common, Q);
            return;
        }
        for (int i = start; i < nu; i++) {
            chosen[pos] = i;
            enumerate(pos + 1, i + 1);
        }
    };
    enumerate(0, 0);
    return total;
}

// Same counter but using a boolean adjacency matrix (for noisy graph scoring)
static long double count_bicliques_in_adj(const vector<vector<bool>>& adj,
                                           int nu, int nl, int P, int Q) {
    if (nu < P || nl < Q) return 0.0L;
    long double total = 0.0L;
    vector<int> chosen(P);
    function<void(int, int)> enumerate = [&](int pos, int start) {
        if (pos == P) {
            int common = 0;
            for (int j = 0; j < nl; j++) {
                bool all = true;
                for (int k = 0; k < P; k++)
                    if (!adj[chosen[k]][j]) { all = false; break; }
                if (all) common++;
            }
            total += nCr_double(common, Q);
            return;
        }
        for (int i = start; i < nu; i++) {
            chosen[pos] = i;
            enumerate(pos + 1, i + 1);
        }
    };
    enumerate(0, 0);
    return total;
}

// ============================================================
// DP-based unbiased estimator
// ============================================================
long double estimate_biclique_count_dp(
    const vector<int>& S_upper, const vector<int>& S_lower,
    int P, int Q,
    BiGraph& g_noisy, BiGraph& g_orig,
    long double p_flip, long double eps1, mt19937& rng) {

    int nu = S_upper.size();
    int nl = S_lower.size();
    if (nu < P || nl < Q) return 0.0L;

    int M = P * Q;
    vector<long double> coeffs(M + 1);
    long double exp_eps = expl(eps1);
    long double denom = powl(exp_eps - 1.0L, M);
    
    for (int k = 0; k <= M; k++) {
        long double term = powl(exp_eps, k);
        if ((M - k) % 2 == 1) term = -term;
        coeffs[k] = term / denom;
    }

    long double total_est = 0.0L;
    
    vector<int> idxs(P);
    iota(idxs.begin(), idxs.end(), 0); // 0, 1, ..., P-1
    
    while (true) {
        vector<int> noisy_degs(nl, 0);
        for (int i = 0; i < P; i++) {
            int u = S_upper[idxs[i]];
            for (int j = 0; j < nl; j++) {
                int v = S_lower[j];
                if (query_noisy_edge(g_orig, g_noisy, u, v, p_flip, rng)) {
                    noisy_degs[j]++;
                }
            }
        }
        
        int max_sum = Q * P;
        vector<vector<long double>> dp(Q + 1, vector<long double>(max_sum + 1, 0.0L));
        dp[0][0] = 1.0L;
        
        for (int d : noisy_degs) {
            for (int c = Q - 1; c >= 0; c--) {
                for (int s = 0; s <= c * P; s++) {
                    if (dp[c][s] > 0.5L) { 
                        if (s + d <= max_sum) {
                            dp[c + 1][s + d] += dp[c][s];
                        }
                    }
                }
            }
        }
        
        for (int k = 0; k <= max_sum; k++) {
            if (dp[Q][k] > 0.5L) {
                total_est += dp[Q][k] * coeffs[k];
            }
        }

        int i = P - 1;
        while (i >= 0 && idxs[i] == nu - P + i) i--;
        if (i < 0) break;
        idxs[i]++;
        for (int j = i + 1; j < P; j++) idxs[j] = idxs[j - 1] + 1;
    }
    
    return total_est;
}

// ============================================================
// Load exact density from SQLite
// ============================================================

bool load_exact_density(const string& dataset, int P, int Q, long double& density) {
    size_t found = dataset.find_last_of("/\\");
    string ds_name = (found != string::npos) ? dataset.substr(found + 1) : dataset;

    sqlite3* db;
    vector<string> db_paths = {
        "../Exact_result/biclique_dds_exact.db",
        "Exact_result/biclique_dds_exact.db",
        "/data/gengdaz/biclique_extension/Exact_result/biclique_dds_exact.db"
    };

    for (auto& path : db_paths) {
        if (sqlite3_open(path.c_str(), &db) == SQLITE_OK) {
            sqlite3_stmt* chk;
            if (sqlite3_prepare_v2(db, "SELECT count(*) FROM exact_results", -1, &chk, nullptr) == SQLITE_OK) {
                sqlite3_finalize(chk);
                const char* sql = "SELECT density FROM exact_results WHERE dataset = ? AND p = ? AND q = ?;";
                sqlite3_stmt* stmt;
                bool found_result = false;
                if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) == SQLITE_OK) {
                    sqlite3_bind_text(stmt, 1, ds_name.c_str(), -1, SQLITE_STATIC);
                    sqlite3_bind_int(stmt, 2, P);
                    sqlite3_bind_int(stmt, 3, Q);
                    if (sqlite3_step(stmt) == SQLITE_ROW) {
                        density = (long double)sqlite3_column_double(stmt, 0);
                        found_result = true;
                    }
                    sqlite3_finalize(stmt);
                }
                sqlite3_close(db);
                if (found_result) return true;
                continue;
            }
            sqlite3_close(db);
        }
    }
    return false;
}

// ============================================================
// Incremental Gain Estimator (One-Round Biclique Technique)
// Estimates the number of NEW (p,q)-bicliques formed by adding v_new
// ============================================================
long double estimate_biclique_gain(
    const vector<int>& S_upper, const vector<int>& S_lower,
    int v_new, bool is_upper_new,
    int P, int Q,
    BiGraph& g_noisy, BiGraph& g_orig,
    long double p_flip, long double eps1, mt19937& rng) {

    const vector<int>& S_fixed = is_upper_new ? S_upper : S_lower;
    const vector<int>& S_target = is_upper_new ? S_lower : S_upper;
    int K_fixed = is_upper_new ? P : Q;
    int K_target = is_upper_new ? Q : P;

    if ((int)S_fixed.size() < K_fixed - 1 || (int)S_target.size() < K_target) return 0.0L;

    int M = K_fixed * K_target;
    vector<long double> coeffs(M + 1);
    long double exp_eps = expl(eps1);
    long double denom = powl(exp_eps - 1.0L, M);
    
    for (int k = 0; k <= M; k++) {
        long double term = nCr_double(M, k) * powl(exp_eps, k);
        if ((M - k) % 2 == 1) term = -term;
        coeffs[k] = term / denom;
    }

    long double total_gain = 0.0L;
    int n_fixed = S_fixed.size();
    int k_needed = K_fixed - 1;

    if (k_needed == 0) {
         vector<int> noisy_degs(S_target.size(), 0);
         for (int j = 0; j < (int)S_target.size(); j++) {
             if (query_noisy_edge(g_orig, g_noisy, v_new, S_target[j], p_flip, rng))
                 noisy_degs[j]++;
         }
         int max_sum = M; 
         vector<vector<long double>> dp(K_target + 1, vector<long double>(max_sum + 1, 0.0L));
         dp[0][0] = 1.0L;
         for (int d : noisy_degs) {
             for (int c = K_target - 1; c >= 0; c--) {
                 for (int s = 0; s <= c * K_fixed; s++) { 
                     if (dp[c][s] > 0.5L && s + d <= max_sum) dp[c + 1][s + d] += dp[c][s];
                 }
             }
         }
         for (int k = 0; k <= max_sum; k++) {
             if (dp[K_target][k] > 0.5L) total_gain += dp[K_target][k] * coeffs[k];
         }
         return total_gain;
    }

    vector<int> idxs(k_needed);
    iota(idxs.begin(), idxs.end(), 0);

    while (true) {
        vector<int> noisy_degs(S_target.size(), 0);
        for (int j = 0; j < (int)S_target.size(); j++) {
            if (query_noisy_edge(g_orig, g_noisy, v_new, S_target[j], p_flip, rng))
                noisy_degs[j]++;
        }
        for (int i = 0; i < k_needed; i++) {
            int u = S_fixed[idxs[i]];
            for (int j = 0; j < (int)S_target.size(); j++) {
                if (query_noisy_edge(g_orig, g_noisy, u, S_target[j], p_flip, rng))
                    noisy_degs[j]++;
            }
        }

        int max_sum = M;
        vector<vector<long double>> dp(K_target + 1, vector<long double>(max_sum + 1, 0.0L));
        dp[0][0] = 1.0L;
        
        for (int d : noisy_degs) {
            for (int c = K_target - 1; c >= 0; c--) {
                for (int s = 0; s <= c * K_fixed; s++) {
                    if (dp[c][s] > 0.5L && s + d <= max_sum) {
                        dp[c + 1][s + d] += dp[c][s];
                    }
                }
            }
        }

        for (int k = 0; k <= max_sum; k++) {
            if (dp[K_target][k] > 0.5L) {
                total_gain += dp[K_target][k] * coeffs[k];
            }
        }

        int i = k_needed - 1;
        while (i >= 0 && idxs[i] == n_fixed - k_needed + i) i--;
        if (i < 0) break;
        idxs[i]++;
        for (int j = i + 1; j < k_needed; j++) idxs[j] = idxs[j - 1] + 1;
    }

    return total_gain;
}

// ============================================================
// Greedy Add using Biclique Estimator (Randomized)
// ============================================================
static pair<vector<int>, long double> greedy_add_biclique(
    BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1,
    const vector<pair<long double, int>>& deg_sorted,
    int max_size, mt19937& rng,
    const vector<int>& seed_upper = {},
    const vector<int>& seed_lower = {}) {

    int n = g.num_nodes();
    vector<int> S;
    vector<bool> in_S(n, false);
    vector<int> S_upper, S_lower;

    // Initialize with seeds
    if (!seed_upper.empty() || !seed_lower.empty()) {
        for (int v : seed_upper) { S.push_back(v); in_S[v] = true; S_upper.push_back(v); }
        for (int v : seed_lower) { S.push_back(v); in_S[v] = true; S_lower.push_back(v); }
    }

    // Fill to minimal size P+Q using top degree (heuristic start)
    int uc = S_upper.size(), lc = S_lower.size();
    for (auto& [d, v] : deg_sorted) {
        if (in_S[v]) continue;
        if (uc >= P && lc >= Q && (int)S.size() >= P + Q) break;
        bool is_up = g.is_upper(v);
        if (is_up && uc >= max(P, (int)S_upper.size() + 1)) continue;
        if (!is_up && lc >= max(Q, (int)S_lower.size() + 1)) continue;
        S.push_back(v); in_S[v] = true;
        if (is_up) { S_upper.push_back(v); uc++; } else { S_lower.push_back(v); lc++; }
    }

    long double current_est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double best_density = current_est / (long double)S.size();
    vector<int> best_set = S;
    int no_improve = 0;

    while ((int)S.size() < max_size && no_improve < 3) {
        // 1. Pre-select candidates based on noisy connectivity (Fast filter)
        vector<pair<int, int>> candidates_pre;
        unordered_set<int> checked;
        vector<int> pool;
        
        // Always check top 100 global degree nodes
        for(int i=0; i<min((int)deg_sorted.size(), 100); ++i) pool.push_back(deg_sorted[i].second);
        
        // Plus random sample
        uniform_int_distribution<int> dist(0, n-1);
        for(int i=0; i<100; ++i) pool.push_back(dist(rng));
        
        for (int v : pool) {
            if (in_S[v] || checked.count(v)) continue;
            checked.insert(v);
            
            bool is_up = g.is_upper(v);
            const vector<int>& opp = is_up ? S_lower : S_upper;
            if (opp.empty()) continue;
            
            int links = 0;
            for (int u : opp) {
                if (query_noisy_edge(g, g_noisy, v, u, p_flip, rng)) links++;
            }
            candidates_pre.push_back({links, v});
        }
        
        sort(candidates_pre.rbegin(), candidates_pre.rend());
        
        // 2. Evaluate Gain for top K candidates
        int K_eval = min((int)candidates_pre.size(), 20);
        vector<pair<long double, int>> candidates_scored;
        
        for (int i = 0; i < K_eval; i++) {
            int v = candidates_pre[i].second;
            bool is_up = g.is_upper(v);
            long double gain = estimate_biclique_gain(S_upper, S_lower, v, is_up, P, Q, g_noisy, g, p_flip, eps1, rng);
            candidates_scored.push_back({gain, v});
        }
        
        sort(candidates_scored.rbegin(), candidates_scored.rend());
        
        if (candidates_scored.empty()) break;

        // 3. Randomized Selection (Top 3)
        int pick_idx = 0;
        if (candidates_scored.size() > 1) {
            uniform_int_distribution<int> pick_dist(0, min((int)candidates_scored.size(), 3) - 1);
            pick_idx = pick_dist(rng);
        }
        
        int best_v = candidates_scored[pick_idx].second;
        long double gain = candidates_scored[pick_idx].first;
        
        S.push_back(best_v); in_S[best_v] = true;
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

// ============================================================
// Greedy Peel using Biclique Estimator
// ============================================================
static pair<vector<int>, long double> greedy_peel_biclique(
    BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1,
    const vector<int>& initial_set,
    mt19937& rng) {

    int n = g.num_nodes();
    vector<int> S = initial_set;
    vector<bool> in_S(n, false);
    vector<int> S_upper, S_lower;
    
    for (int v : S) {
        in_S[v] = true;
        if (g.is_upper(v)) S_upper.push_back(v);
        else S_lower.push_back(v);
    }

    long double current_est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double best_density = (S.size() > 0) ? current_est / (long double)S.size() : 0.0L;
    vector<int> best_set = S;

    // Peel until P+Q
    while ((int)S.size() > P + Q) {
        // Find vertex whose removal causes MINIMUM loss (or max gain if negative)
        // Loss(v) = estimate_biclique_gain(S \ v, v)
        // We want to minimize this Loss.
        
        int check_limit = (S.size() > 100) ? 30 : ((S.size() > 50) ? 40 : S.size());
        vector<int> candidates_to_check = S;
        if (S.size() > 500) {
            shuffle(candidates_to_check.begin(), candidates_to_check.end(), rng);
            candidates_to_check.resize(check_limit);
        }
        
        int best_v_rem = -1;
        long double min_loss = 1e30L;
        
        for (int v : candidates_to_check) {
            bool is_up = g.is_upper(v);
            vector<int> tu, tl;
            if (is_up) {
                for(int u : S_upper) if(u != v) tu.push_back(u);
                tl = S_lower;
            } else {
                tu = S_upper;
                for(int u : S_lower) if(u != v) tl.push_back(u);
            }
            
            if ((int)tu.size() < P || (int)tl.size() < Q) continue;

            long double loss = estimate_biclique_gain(tu, tl, v, is_up, P, Q, g_noisy, g, p_flip, eps1, rng);
            
            if (loss < min_loss) {
                min_loss = loss;
                best_v_rem = v;
            }
        }

        if (best_v_rem != -1) {
            int v_rem = best_v_rem;
            
            for (size_t i = 0; i < S.size(); i++) {
                if (S[i] == v_rem) { S.erase(S.begin() + i); break; }
            }
            if (g.is_upper(v_rem)) {
                for (size_t i = 0; i < S_upper.size(); i++) {
                    if (S_upper[i] == v_rem) { S_upper.erase(S_upper.begin() + i); break; }
                }
            } else {
                for (size_t i = 0; i < S_lower.size(); i++) {
                    if (S_lower[i] == v_rem) { S_lower.erase(S_lower.begin() + i); break; }
                }
            }
            in_S[v_rem] = false;

            current_est -= min_loss;
            long double cur_d = current_est / (long double)S.size();
            
            if (cur_d > best_density) {
                best_density = cur_d;
                best_set = S;
            }
        } else {
            break;
        }
    }

    return {best_set, best_density};
}

// ============================================================
// Greedy Shrink using DP estimator
// ============================================================
static pair<vector<int>, long double> greedy_shrink_dp(
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

    long double est = estimate_biclique_count_dp(
        S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double cur_d = est / (long double)S.size();
    long double best_d = cur_d;
    vector<int> best_set = S;

    bool improved = true;
    while (improved && (int)S.size() > P + Q) {
        improved = false;
        long double best_rem_d = cur_d;
        int best_rem_idx = -1;

        // Optimization: if S is large (> 50), don't check every removal.
        // Check removal of vertices with lowest noisy degree?
        // For now, stick to full check but maybe limit if too slow.
        
        for (int idx = 0; idx < (int)S.size(); idx++) {
            int v = S[idx];
            vector<int> tu, tl;
            for (int u : S_upper) if (u != v) tu.push_back(u);
            for (int u : S_lower) if (u != v) tl.push_back(u);
            if ((int)tu.size() < P || (int)tl.size() < Q) continue;

            long double ne = estimate_biclique_count_dp(
                tu, tl, P, Q, g_noisy, g, p_flip, eps1, rng);
            long double nd = ne / (long double)(S.size() - 1);
            if (nd > best_rem_d) { best_rem_d = nd; best_rem_idx = idx; }
        }

        if (best_rem_idx >= 0) {
            S.erase(S.begin() + best_rem_idx);
            S_upper.clear(); S_lower.clear();
            for (int u : S) {
                if (g.is_upper(u)) S_upper.push_back(u);
                else S_lower.push_back(u);
            }
            cur_d = best_rem_d;
            improved = true;
            if (cur_d > best_d) { best_d = cur_d; best_set = S; }
        }
    }
    return {best_set, best_d};
}

// ============================================================
// Charikar-style Edge Density Peeling on Debiased Noisy Graph
// Returns the set with highest debiased edge density = |E_debiased(S)| / |S|
// ============================================================
static pair<vector<int>, long double> charikar_edge_density_peel(
    BiGraph& g, BiGraph& g_noisy,
    const vector<int>& candidate_set,
    long double p_flip, mt19937& rng) {

    int nc = candidate_set.size();
    if (nc <= 2) return {candidate_set, 0.0L};

    unordered_map<int, int> id_map;
    for (int i = 0; i < nc; i++) id_map[candidate_set[i]] = i;

    vector<bool> active(nc, true);
    vector<long double> deg(nc, 0.0L);
    long double denom_rr = 1.0L - 2.0L * p_flip;

    long double total_debiased_edges = 0.0L;
    int active_count = nc;

    for (int i = 0; i < nc; i++) {
        int u = candidate_set[i];
        for (int j = i + 1; j < nc; j++) {
            int v = candidate_set[j];
            if ((g.is_upper(u) && !g.is_upper(v)) || (!g.is_upper(u) && g.is_upper(v))) {
                bool noisy = query_noisy_edge(g, g_noisy, u, v, p_flip, rng);
                long double debiased = (noisy ? 1.0L - p_flip : -p_flip) / (denom_rr + 1e-15L);
                deg[i] += debiased;
                deg[j] += debiased;
                total_debiased_edges += debiased;
            }
        }
    }

    long double best_density = total_debiased_edges / (long double)active_count;
    vector<bool> best_active = active;
    int best_count = active_count;

    for (int iter = 0; iter < nc - 2; iter++) {
        int min_idx = -1;
        long double min_deg = 1e30L;
        for (int i = 0; i < nc; i++) {
            if (active[i] && deg[i] < min_deg) {
                min_deg = deg[i];
                min_idx = i;
            }
        }
        if (min_idx < 0) break;

        active[min_idx] = false;
        active_count--;
        int u = candidate_set[min_idx];

        for (int j = 0; j < nc; j++) {
            if (!active[j] || j == min_idx) continue;
            int v = candidate_set[j];
            if ((g.is_upper(u) && !g.is_upper(v)) || (!g.is_upper(u) && g.is_upper(v))) {
                bool noisy = query_noisy_edge(g, g_noisy, u, v, p_flip, rng);
                long double debiased = (noisy ? 1.0L - p_flip : -p_flip) / (denom_rr + 1e-15L);
                deg[j] -= debiased;
                total_debiased_edges -= debiased;
            }
        }

        if (active_count > 0) {
            long double cur_density = total_debiased_edges / (long double)active_count;
            if (cur_density > best_density) {
                best_density = cur_density;
                best_active = active;
                best_count = active_count;
            }
        }
    }

    vector<int> result_set;
    for (int i = 0; i < nc; i++) {
        if (best_active[i]) result_set.push_back(candidate_set[i]);
    }
    return {result_set, best_density};
}

// ============================================================
// ABCore Peeling: find (alpha, beta)-core on debiased noisy graph
// Iteratively remove upper vertices with debiased degree < alpha
// and lower vertices with debiased degree < beta
// ============================================================
static vector<int> abcore_peel(
    BiGraph& g, BiGraph& g_noisy,
    const vector<int>& candidate_set,
    long double alpha, long double beta,
    long double p_flip, mt19937& rng) {

    int nc = candidate_set.size();
    if (nc <= 2) return candidate_set;

    long double denom_rr = 1.0L - 2.0L * p_flip;
    vector<bool> active(nc, true);

    vector<int> upper_idx, lower_idx;
    for (int i = 0; i < nc; i++) {
        if (g.is_upper(candidate_set[i])) upper_idx.push_back(i);
        else lower_idx.push_back(i);
    }

    bool changed = true;
    while (changed) {
        changed = false;
        for (int i = 0; i < nc; i++) {
            if (!active[i]) continue;
            int u = candidate_set[i];
            bool is_up = g.is_upper(u);

            long double debiased_deg = 0.0L;
            for (int j = 0; j < nc; j++) {
                if (!active[j] || j == i) continue;
                int v = candidate_set[j];
                if (is_up == g.is_upper(v)) continue;
                bool noisy = query_noisy_edge(g, g_noisy, u, v, p_flip, rng);
                long double debiased = (noisy ? 1.0L - p_flip : -p_flip) / (denom_rr + 1e-15L);
                debiased_deg += debiased;
            }

            long double threshold = is_up ? alpha : beta;
            if (debiased_deg < threshold) {
                active[i] = false;
                changed = true;
            }
        }
    }

    vector<int> result;
    for (int i = 0; i < nc; i++) {
        if (active[i]) result.push_back(candidate_set[i]);
    }
    return result;
}

// ============================================================
// Greedy Add with Neighbor-Based Candidate Selection
// Only considers neighbors of S in the noisy graph
// ============================================================
static pair<vector<int>, long double> greedy_add_neighbor_based(
    BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1,
    const vector<int>& seed_set,
    int max_size, mt19937& rng) {

    int n = g.num_nodes();
    vector<int> S = seed_set;
    vector<bool> in_S(n, false);
    vector<int> S_upper, S_lower;

    for (int v : S) {
        in_S[v] = true;
        if (g.is_upper(v)) S_upper.push_back(v);
        else S_lower.push_back(v);
    }

    if ((int)S_upper.size() < P || (int)S_lower.size() < Q) {
        return {S, 0.0L};
    }

    long double current_est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double best_density = current_est / (long double)S.size();
    vector<int> best_set = S;
    int no_improve = 0;

    while ((int)S.size() < max_size && no_improve < 5) {
        // Collect candidates: noisy neighbors of S only (ELDP compliant)
        unordered_set<int> candidate_set;
        for (int v : S) {
            // Check cached noisy edges
            for (auto& [key, val] : g_noisy.edge_vector[v]) {
                if (val && !in_S[key]) candidate_set.insert(key);
            }
        }
        // If too few candidates, query noisy edges for S vs all vertices
        if ((int)candidate_set.size() < 20) {
            for (int v : S) {
                for (int u = 0; u < n; u++) {
                    if (in_S[u]) continue;
                    if ((g.is_upper(v) && !g.is_upper(u)) || (!g.is_upper(v) && g.is_upper(u))) {
                        if (query_noisy_edge(g, g_noisy, v, u, p_flip, rng)) {
                            candidate_set.insert(u);
                        }
                    }
                }
            }
        }

        // Score candidates by noisy connectivity to S
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

        // Evaluate top-K by biclique gain
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

        // Pick best (or top-2 random for diversity)
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

// ============================================================
// Full-evaluation Greedy Add (like triangle DDS)
// Evaluates ALL candidate vertices at each step
// ============================================================
static pair<vector<int>, long double> greedy_add_full_eval(
    BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1,
    const vector<int>& initial_set,
    mt19937& rng, int max_size = 200, int patience = 5) {

    int n = g.num_nodes();
    vector<int> S = initial_set;
    vector<bool> in_S(n, false);
    vector<int> S_upper, S_lower;
    for (int v : S) {
        in_S[v] = true;
        if (g.is_upper(v)) S_upper.push_back(v);
        else S_lower.push_back(v);
    }

    if ((int)S_upper.size() < P || (int)S_lower.size() < Q) return {S, 0.0L};

    long double current_est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double current_density = current_est / (long double)S.size();
    long double best_density = current_density;
    vector<int> best_set = S;
    int no_improve = 0;

    while ((int)S.size() < min(n, max_size) && no_improve < patience) {
        long double best_new_density = -1e30L;
        int best_v = -1;
        long double best_gain = 0;

        for (int v = 0; v < n; v++) {
            if (in_S[v]) continue;
            bool is_up = g.is_upper(v);
            if (is_up && (int)S_upper.size() < P - 1) continue;
            if (!is_up && (int)S_lower.size() < Q - 1) continue;

            long double gain = estimate_biclique_gain(S_upper, S_lower, v, is_up,
                                                       P, Q, g_noisy, g, p_flip, eps1, rng);
            long double new_d = (current_est + gain) / (long double)(S.size() + 1);
            if (new_d > best_new_density) {
                best_new_density = new_d;
                best_v = v;
                best_gain = gain;
            }
        }

        if (best_v < 0) break;

        S.push_back(best_v);
        in_S[best_v] = true;
        if (g.is_upper(best_v)) S_upper.push_back(best_v);
        else S_lower.push_back(best_v);
        current_est += best_gain;
        current_density = current_est / (long double)S.size();

        if (current_density > best_density) {
            best_density = current_density;
            best_set = S;
            no_improve = 0;
        } else {
            no_improve++;
        }
    }

    return {best_set, best_density};
}

// ============================================================
// Full-evaluation Greedy Shrink
// ============================================================
static pair<vector<int>, long double> greedy_shrink_full_eval(
    BiGraph& g, BiGraph& g_noisy, int P, int Q,
    long double p_flip, long double eps1,
    const vector<int>& initial_set,
    mt19937& rng) {

    int n = g.num_nodes();
    vector<int> S = initial_set;
    vector<bool> in_S(n, false);
    vector<int> S_upper, S_lower;
    for (int v : S) {
        in_S[v] = true;
        if (g.is_upper(v)) S_upper.push_back(v);
        else S_lower.push_back(v);
    }

    if ((int)S_upper.size() < P || (int)S_lower.size() < Q) return {S, 0.0L};

    long double current_est = estimate_biclique_count_dp(S_upper, S_lower, P, Q, g_noisy, g, p_flip, eps1, rng);
    long double current_density = current_est / (long double)S.size();
    long double best_density = current_density;
    vector<int> best_set = S;

    bool improved = true;
    while (improved && (int)S.size() > P + Q) {
        improved = false;
        long double best_rem_density = current_density;
        int best_rem_v = -1;
        long double best_loss = 0;

        for (int v : S) {
            bool is_up = g.is_upper(v);
            vector<int> tu, tl;
            if (is_up) {
                for (int u : S_upper) if (u != v) tu.push_back(u);
                tl = S_lower;
            } else {
                tu = S_upper;
                for (int u : S_lower) if (u != v) tl.push_back(u);
            }
            if ((int)tu.size() < P || (int)tl.size() < Q) continue;

            long double loss = estimate_biclique_gain(tu, tl, v, is_up, P, Q, g_noisy, g, p_flip, eps1, rng);
            long double new_d = (current_est - loss) / (long double)(S.size() - 1);
            if (new_d > best_rem_density) {
                best_rem_density = new_d;
                best_rem_v = v;
                best_loss = loss;
                improved = true;
            }
        }

        if (improved && best_rem_v >= 0) {
            S.erase(remove(S.begin(), S.end(), best_rem_v), S.end());
            in_S[best_rem_v] = false;
            if (g.is_upper(best_rem_v)) {
                S_upper.erase(remove(S_upper.begin(), S_upper.end(), best_rem_v), S_upper.end());
            } else {
                S_lower.erase(remove(S_lower.begin(), S_lower.end(), best_rem_v), S_lower.end());
            }
            current_est -= best_loss;
            current_density = best_rem_density;
            if (current_density > best_density) {
                best_density = current_density;
                best_set = S;
            }
        }
    }
    return {best_set, best_density};
}

// ============================================================
// Subsample: take top-k vertices from set by debiased degree within set
// (using noisy edge queries). Returns upper and lower separately.
// ============================================================
static vector<int> subsample_by_noisy_degree(
    BiGraph& g, BiGraph& g_noisy, const vector<int>& full_set,
    int target_size, long double p_flip, mt19937& rng) {

    if ((int)full_set.size() <= target_size) return full_set;

    unordered_set<int> fs(full_set.begin(), full_set.end());
    vector<int> upper_v, lower_v;
    for (int v : full_set) {
        if (g.is_upper(v)) upper_v.push_back(v);
        else lower_v.push_back(v);
    }

    // Compute noisy degree within set for each vertex
    vector<pair<long double, int>> scored;
    long double denom_rr = 1.0L - 2.0L * p_flip;
    for (int v : full_set) {
        bool is_up = g.is_upper(v);
        const vector<int>& opp = is_up ? lower_v : upper_v;
        int noisy_links = 0;
        for (int u : opp) {
            if (query_noisy_edge(g, g_noisy, v, u, p_flip, rng)) noisy_links++;
        }
        long double debiased = ((long double)noisy_links - p_flip * opp.size()) / (denom_rr + 1e-15L);
        scored.push_back({debiased, v});
    }
    sort(scored.rbegin(), scored.rend());

    vector<int> result;
    for (int i = 0; i < target_size && i < (int)scored.size(); i++)
        result.push_back(scored[i].second);
    return result;
}

// ============================================================
// Butterfly-proxy score: for (P,Q)-bicliques, estimate the
// butterfly-like density using noisy common neighbor counts.
// For (2,Q): score = sum_{(u1,u2)} C(debiased_CN(u1,u2), Q)
// ============================================================
static long double butterfly_proxy_score(
    BiGraph& g, BiGraph& g_noisy,
    const vector<int>& S,
    int P, int Q,
    long double p_flip, mt19937& rng) {

    vector<int> su, sl;
    for (int v : S) { if (g.is_upper(v)) su.push_back(v); else sl.push_back(v); }
    if ((int)su.size() < P || (int)sl.size() < Q) return -1e30L;

    long double pp = p_flip * p_flip;
    long double qq = (1.0L - p_flip) * (1.0L - p_flip);
    int nl = sl.size();

    long double total_score = 0;

    if (P == 2) {
        for (int i = 0; i < (int)su.size(); i++) {
            for (int j = i + 1; j < (int)su.size(); j++) {
                int common = 0;
                for (int v : sl) {
                    bool e1 = query_noisy_edge(g, g_noisy, su[i], v, p_flip, rng);
                    bool e2 = query_noisy_edge(g, g_noisy, su[j], v, p_flip, rng);
                    if (e1 && e2) common++;
                }
                long double est_cn = ((long double)common - pp * nl) / (qq - pp + 1e-15L);
                est_cn = max(0.0L, est_cn);
                total_score += nCr_double((int)round(est_cn), Q);
            }
        }
    } else {
        // General case: use DP estimator (slower but correct)
        total_score = estimate_biclique_count_dp(su, sl, P, Q, g_noisy, g, p_flip,
            logl(1.0L / p_flip - 1.0L), rng);
    }

    return total_score / (long double)S.size();
}

// ============================================================
// Main entry point — v11: Multi-Strategy + Butterfly Peeling + Refinement
// ============================================================

DDSResult eldp_biclique_dds(BiGraph& g, int P, int Q,
    long double epsilon, unsigned long seed) {

    int n = g.num_nodes();
    int n1 = g.num_v1, n2 = n - n1;
    DDSResult result;

    long double eps0 = epsilon * 0.05L;
    long double eps1 = epsilon - eps0;
    long double p_flip = 1.0L / (expl(eps1) + 1.0L);

    cout << "  [ELDP v11] eps=" << fixed << setprecision(4) << epsilon
         << " p_flip=" << setprecision(4) << p_flip
         << " n1=" << n1 << " n2=" << n2 << endl;

    mt19937 rng(seed);
    normal_distribution<double> lap_dist(0.0, 1.0 / (double)eps0);
    vector<long double> deg_est(n);
    for (int i = 0; i < n; i++)
        deg_est[i] = (long double)g.neighbor[i].size() + lap_dist(rng);

    BiGraph g_noisy(g);
    g_noisy.edge_vector.resize(n);
    for (int i = 0; i < n; i++) g_noisy.edge_vector[i].clear();

    vector<pair<long double, int>> deg_sorted(n);
    for (int i = 0; i < n; i++) deg_sorted[i] = {deg_est[i], i};
    sort(deg_sorted.rbegin(), deg_sorted.rend());

    // Adaptive vertex pool: larger for small graphs, capped for big ones
    int ku = min(n1, max(100, min(n1, 300)));
    int kl = min(n2, max(200, min(n2, 500)));
    vector<int> top_upper, top_lower;
    for (auto& [d, v] : deg_sorted) {
        if (g.is_upper(v) && (int)top_upper.size() < ku) top_upper.push_back(v);
        if (!g.is_upper(v) && (int)top_lower.size() < kl) top_lower.push_back(v);
    }
    ku = top_upper.size();
    kl = top_lower.size();

    cout << "  [Pool] ku=" << ku << " kl=" << kl << endl;

    // Map from vertex id to index in top_upper/top_lower
    unordered_map<int, int> upper_idx_map, lower_idx_map;
    for (int i = 0; i < ku; i++) upper_idx_map[top_upper[i]] = i;
    for (int i = 0; i < kl; i++) lower_idx_map[top_lower[i]] = i;

    // Build noisy adjacency matrix
    vector<vector<bool>> nadj(ku, vector<bool>(kl, false));
    for (int i = 0; i < ku; i++)
        for (int j = 0; j < kl; j++)
            nadj[i][j] = query_noisy_edge(g, g_noisy, top_upper[i], top_lower[j], p_flip, rng);

    long double pp = p_flip * p_flip;
    long double qq = (1.0L - p_flip) * (1.0L - p_flip);

    // Track best true density (for evaluation only)
    long double overall_best_d = 0.0L;
    vector<int> overall_best_set;
    int ncand_eval = 0;

    auto add_candidate = [&](const vector<int>& s, const string& label) {
        if (s.empty()) return;
        vector<int> su, sl;
        for (int v : s) { if (g.is_upper(v)) su.push_back(v); else sl.push_back(v); }
        if ((int)su.size() < P || (int)sl.size() < Q) return;
        long double bc = count_bicliques_in_set(g, s, P, Q);
        long double td = bc / (long double)s.size();
        
        // Record trace for variance analysis
        write_trace(s.size(), td, td, label);
        
        ncand_eval++;
        if (td > overall_best_d) {
            overall_best_d = td;
            overall_best_set = s;
            cout << "  [NEW BEST] " << label << " |S|=" << s.size()
                 << " (u=" << su.size() << ",l=" << sl.size() << ")"
                 << " d=" << fixed << setprecision(4) << td << endl;
        }
    };

    // =========================================================
    // Strategy A: Butterfly-aware greedy peeling
    // Start with all top vertices, peel lowest-contribution vertex
    // =========================================================
    {
        // Compute noisy biclique count for the full set and per-vertex contribution
        // For (P,Q) bicliques: contribution of upper vertex u_i =
        //   sum over C(other_upper, P-1) subsets containing u_i, times C(common_lower, Q)
        // For P=2 this simplifies to: sum_{j!=i} C(noisy_CN(i,j), Q) for upper
        //                              similarly for lower if Q=2

        // Track which vertices are in the set
        vector<bool> in_set_u(ku, true), in_set_l(kl, true);
        int cur_ku = ku, cur_kl = kl;

        // For P=2: precompute pairwise common neighbor counts (noisy)
        // cn[i][j] = # lower vertices in set connected to both upper i and j
        vector<vector<int>> cn_uu; // common lower neighbors for upper pairs
        if (P == 2) {
            cn_uu.resize(ku, vector<int>(ku, 0));
            for (int i = 0; i < ku; i++) {
                for (int j = i + 1; j < ku; j++) {
                    int c = 0;
                    for (int l = 0; l < kl; l++)
                        if (nadj[i][l] && nadj[j][l]) c++;
                    cn_uu[i][j] = cn_uu[j][i] = c;
                }
            }
        }

        // Compute butterfly contribution for each vertex
        auto compute_upper_contribution = [&](int i) -> long double {
            if (!in_set_u[i]) return 0.0L;
            if (P == 2) {
                long double score = 0;
                for (int j = 0; j < ku; j++) {
                    if (j == i || !in_set_u[j]) continue;
                    score += nCr_double(cn_uu[i][j], Q);
                }
                return score;
            }
            return 0.0L; // general P handled differently
        };

        auto compute_lower_contribution = [&](int l) -> long double {
            if (!in_set_l[l]) return 0.0L;
            if (P == 2 && Q == 2) {
                // For each upper pair (i,j) both connected to l:
                // removing l reduces cn_uu[i][j] by 1, so butterfly loss =
                // C(cn[i][j], 2) - C(cn[i][j]-1, 2) = cn[i][j] - 1
                // Total contribution of l = sum_{(i,j): both connected to l} (cn_uu[i][j] - 1)
                long double score = 0;
                vector<int> connected_upper;
                for (int i = 0; i < ku; i++)
                    if (in_set_u[i] && nadj[i][l]) connected_upper.push_back(i);
                for (int a = 0; a < (int)connected_upper.size(); a++) {
                    for (int b = a + 1; b < (int)connected_upper.size(); b++) {
                        score += max(0, cn_uu[connected_upper[a]][connected_upper[b]] - 1);
                    }
                }
                return score;
            }
            if (P == 2) {
                // For (2,Q>2): contribution of l to butterflies
                // For each upper pair connected to l:
                // removing l changes C(cn,Q) to C(cn-1,Q), loss = C(cn-1,Q-1)
                long double score = 0;
                vector<int> connected_upper;
                for (int i = 0; i < ku; i++)
                    if (in_set_u[i] && nadj[i][l]) connected_upper.push_back(i);
                for (int a = 0; a < (int)connected_upper.size(); a++) {
                    for (int b = a + 1; b < (int)connected_upper.size(); b++) {
                        score += nCr_double(cn_uu[connected_upper[a]][connected_upper[b]] - 1, Q - 1);
                    }
                }
                return score;
            }
            return 0.0L;
        };

        // Greedy peeling: remove vertex with lowest contribution/degree ratio
        long double total_bc = count_bicliques_in_adj(nadj, ku, kl, P, Q);
        int total_verts = ku + kl;
        long double best_peel_d = total_bc / (long double)total_verts;
        vector<bool> best_u_mask = in_set_u, best_l_mask = in_set_l;

        int peel_steps = 0;
        int max_peel = (ku + kl) - (P + Q);

        while (peel_steps < max_peel && cur_ku >= P && cur_kl >= Q) {
            long double worst_contrib = 1e30L;
            int worst_idx = -1;
            bool worst_is_upper = true;

            for (int i = 0; i < ku; i++) {
                if (!in_set_u[i]) continue;
                if (cur_ku <= P) continue; // can't remove if at minimum
                long double c = compute_upper_contribution(i);
                long double ratio = c / (long double)(cur_ku + cur_kl);
                if (ratio < worst_contrib) {
                    worst_contrib = ratio;
                    worst_idx = i;
                    worst_is_upper = true;
                }
            }
            for (int l = 0; l < kl; l++) {
                if (!in_set_l[l]) continue;
                if (cur_kl <= Q) continue;
                long double c = compute_lower_contribution(l);
                long double ratio = c / (long double)(cur_ku + cur_kl);
                if (ratio < worst_contrib) {
                    worst_contrib = ratio;
                    worst_idx = l;
                    worst_is_upper = false;
                }
            }

            if (worst_idx < 0) break;

            // Remove the worst vertex
            if (worst_is_upper) {
                in_set_u[worst_idx] = false;
                cur_ku--;
                // Update cn_uu if P=2
                if (P == 2) {
                    for (int j = 0; j < ku; j++)
                        cn_uu[worst_idx][j] = cn_uu[j][worst_idx] = 0;
                }
            } else {
                in_set_l[worst_idx] = false;
                cur_kl--;
                // Update cn_uu if P=2: decrease cn for upper pairs
                if (P == 2) {
                    for (int i = 0; i < ku; i++) {
                        if (!in_set_u[i] || !nadj[i][worst_idx]) continue;
                        for (int j = i + 1; j < ku; j++) {
                            if (!in_set_u[j] || !nadj[j][worst_idx]) continue;
                            cn_uu[i][j]--;
                            cn_uu[j][i]--;
                        }
                    }
                }
            }

            // Recompute noisy biclique count for current set
            // (For P=2, use fast pair-based count)
            long double cur_bc = 0;
            if (P == 2) {
                for (int i = 0; i < ku; i++) {
                    if (!in_set_u[i]) continue;
                    for (int j = i + 1; j < ku; j++) {
                        if (!in_set_u[j]) continue;
                        cur_bc += nCr_double(cn_uu[i][j], Q);
                    }
                }
            }
            total_verts = cur_ku + cur_kl;
            long double cur_d = (total_verts > 0) ? cur_bc / (long double)total_verts : 0;

            if (cur_d > best_peel_d) {
                best_peel_d = cur_d;
                best_u_mask = in_set_u;
                best_l_mask = in_set_l;
            }
            peel_steps++;
        }

        // Collect the best peeled set
        vector<int> peel_set;
        for (int i = 0; i < ku; i++) if (best_u_mask[i]) peel_set.push_back(top_upper[i]);
        for (int i = 0; i < kl; i++) if (best_l_mask[i]) peel_set.push_back(top_lower[i]);
        add_candidate(peel_set, "BtfPeel");

        // Also try subsets at different sizes
        // Re-run peeling from best set to get even smaller dense cores
        for (int target : {P + Q, 5, 8, 10, 15, 20, 30, 50, 80}) {
            if (target >= (int)peel_set.size() || target < P + Q) continue;
            vector<int> sub = peel_set;
            // Quick greedy peel on sub to target size
            while ((int)sub.size() > target) {
                // Remove vertex with lowest noisy degree within set
                long double worst_score = 1e30L;
                int worst_pos = 0;
                for (int i = 0; i < (int)sub.size(); i++) {
                    int v = sub[i];
                    bool is_up = g.is_upper(v);
                    // Count how many vertices of opposite side this vertex connects to
                    int links = 0;
                    for (int w : sub) {
                        if (g.is_upper(w) == is_up) continue;
                        if (query_noisy_edge(g, g_noisy, v, w, p_flip, rng)) links++;
                    }
                    if ((long double)links < worst_score) {
                        worst_score = links;
                        worst_pos = i;
                    }
                }
                sub.erase(sub.begin() + worst_pos);
            }
            add_candidate(sub, "BtfPeel-sub" + to_string(target));
        }
    }

    // =========================================================
    // Strategy B: Pair-centric seeds + greedy grow
    // =========================================================
    if (P == 2) {
        // Compute debiased common neighbors for upper pairs
        vector<pair<long double, pair<int, int>>> pair_scores;
        for (int i = 0; i < ku; i++) {
            for (int j = i + 1; j < ku; j++) {
                int common = 0;
                for (int l = 0; l < kl; l++)
                    if (nadj[i][l] && nadj[j][l]) common++;
                long double est_cn = ((long double)common - pp * kl) / (qq - pp + 1e-15L);
                pair_scores.push_back({est_cn, {i, j}});
            }
        }
        sort(pair_scores.rbegin(), pair_scores.rend());

        cout << "  [Pairs] Top 5 CN: " << fixed << setprecision(1);
        for (int i = 0; i < min(5, (int)pair_scores.size()); i++)
            cout << pair_scores[i].first << " ";
        cout << endl;

        // For top pairs: build seed = pair + common noisy neighbors, then greedily grow
        int n_pairs = min(20, (int)pair_scores.size());
        for (int pi = 0; pi < n_pairs; pi++) {
            int ui = pair_scores[pi].second.first;
            int uj = pair_scores[pi].second.second;

            // Seed: the pair + their common noisy lower neighbors
            vector<int> seed_set = {top_upper[ui], top_upper[uj]};
            for (int l = 0; l < kl; l++) {
                if (nadj[ui][l] && nadj[uj][l])
                    seed_set.push_back(top_lower[l]);
            }

            add_candidate(seed_set, "PairSeed-" + to_string(pi));

            // Greedy grow: add vertices that increase noisy biclique density
            // Use the noisy adjacency matrix for fast evaluation
            vector<bool> in_cand(ku + kl, false);
            vector<int> cand_upper_idx, cand_lower_idx;
            for (int v : seed_set) {
                auto it_u = upper_idx_map.find(v);
                if (it_u != upper_idx_map.end()) { in_cand[it_u->second] = true; cand_upper_idx.push_back(it_u->second); }
                auto it_l = lower_idx_map.find(v);
                if (it_l != lower_idx_map.end()) { in_cand[ku + it_l->second] = true; cand_lower_idx.push_back(it_l->second); }
            }

            auto noisy_btf_count = [&]() -> long double {
                long double bc = 0;
                for (int a = 0; a < (int)cand_upper_idx.size(); a++) {
                    for (int b = a + 1; b < (int)cand_upper_idx.size(); b++) {
                        int cn = 0;
                        for (int l : cand_lower_idx) {
                            if (nadj[cand_upper_idx[a]][l] && nadj[cand_upper_idx[b]][l]) cn++;
                        }
                        bc += nCr_double(cn, Q);
                    }
                }
                return bc;
            };

            long double cur_bc = noisy_btf_count();
            int cur_size = seed_set.size();
            long double cur_d = (cur_size > 0) ? cur_bc / (long double)cur_size : 0;
            int no_improve = 0;

            while (no_improve < 3 && cur_size < min(200, ku + kl)) {
                long double best_new_d = cur_d;
                int best_add_idx = -1;
                bool best_is_upper = true;

                // Try adding each upper vertex not in set
                for (int i = 0; i < ku; i++) {
                    if (in_cand[i]) continue;
                    // Compute gain: new bicliques from adding this upper vertex
                    long double gain = 0;
                    for (int j : cand_upper_idx) {
                        int cn = 0;
                        for (int l : cand_lower_idx)
                            if (nadj[i][l] && nadj[j][l]) cn++;
                        gain += nCr_double(cn, Q);
                    }
                    long double new_d = (cur_bc + gain) / (long double)(cur_size + 1);
                    if (new_d > best_new_d) {
                        best_new_d = new_d;
                        best_add_idx = i;
                        best_is_upper = true;
                    }
                }

                // Try adding each lower vertex not in set
                for (int l = 0; l < kl; l++) {
                    if (in_cand[ku + l]) continue;
                    // Gain: for each upper pair, if both connected to l, gain = C(cn+1,Q)-C(cn,Q)
                    // which equals C(cn, Q-1) for the pairs connected to l
                    long double gain = 0;
                    if (Q == 2) {
                        // Adding lower vertex l: for each upper pair (a,b) both connected to l in noisy graph
                        // gain = cn_before + 1 (for each pair, old C(cn,2)→C(cn+1,2), diff = cn)
                        // But cn here is among the current candidate lower vertices
                        int connected_count = 0;
                        vector<int> conn_upper;
                        for (int i : cand_upper_idx)
                            if (nadj[i][l]) conn_upper.push_back(i);
                        for (int a = 0; a < (int)conn_upper.size(); a++) {
                            for (int b = a + 1; b < (int)conn_upper.size(); b++) {
                                // Current cn for this pair among existing lower
                                int cn = 0;
                                for (int ll : cand_lower_idx)
                                    if (nadj[conn_upper[a]][ll] && nadj[conn_upper[b]][ll]) cn++;
                                gain += cn; // C(cn+1,2)-C(cn,2) = cn
                            }
                        }
                    } else {
                        // General Q: compute gain directly
                        for (int a = 0; a < (int)cand_upper_idx.size(); a++) {
                            for (int b = a + 1; b < (int)cand_upper_idx.size(); b++) {
                                if (!nadj[cand_upper_idx[a]][l] || !nadj[cand_upper_idx[b]][l]) continue;
                                int cn = 0;
                                for (int ll : cand_lower_idx)
                                    if (nadj[cand_upper_idx[a]][ll] && nadj[cand_upper_idx[b]][ll]) cn++;
                                gain += nCr_double(cn, Q - 1); // C(cn+1,Q)-C(cn,Q)=C(cn,Q-1)
                            }
                        }
                    }
                    long double new_d = (cur_bc + gain) / (long double)(cur_size + 1);
                    if (new_d > best_new_d) {
                        best_new_d = new_d;
                        best_add_idx = l;
                        best_is_upper = false;
                    }
                }

                if (best_add_idx < 0 || best_new_d <= cur_d) {
                    no_improve++;
                    continue;
                }

                // Add the best vertex
                if (best_is_upper) {
                    in_cand[best_add_idx] = true;
                    cand_upper_idx.push_back(best_add_idx);
                    seed_set.push_back(top_upper[best_add_idx]);
                } else {
                    in_cand[ku + best_add_idx] = true;
                    cand_lower_idx.push_back(best_add_idx);
                    seed_set.push_back(top_lower[best_add_idx]);
                }
                cur_bc = noisy_btf_count(); // recompute exactly
                cur_size = seed_set.size();
                cur_d = cur_bc / (long double)cur_size;
                no_improve = 0;
            }

            add_candidate(seed_set, "PairGrow-" + to_string(pi));

            // Greedy shrink: remove worst vertex
            while (cur_size > P + Q) {
                long double best_shrink_d = -1;
                int best_rem_pos = -1;

                // Try removing each vertex
                for (int idx = 0; idx < (int)seed_set.size(); idx++) {
                    int v = seed_set[idx];
                    bool is_up = g.is_upper(v);
                    int up_count = cand_upper_idx.size() - (is_up ? 1 : 0);
                    int lo_count = cand_lower_idx.size() - (is_up ? 0 : 1);
                    if (up_count < P || lo_count < Q) continue;

                    // Temporarily remove and compute density
                    // (expensive but set is small after grow)
                    vector<int> temp = seed_set;
                    temp.erase(temp.begin() + idx);
                    long double bc_temp = count_bicliques_in_set(g, temp, P, Q);
                    long double d_temp = bc_temp / (long double)temp.size();
                    if (d_temp > best_shrink_d) {
                        best_shrink_d = d_temp;
                        best_rem_pos = idx;
                    }
                }

                if (best_shrink_d <= overall_best_d || best_rem_pos < 0) break;

                int rem_v = seed_set[best_rem_pos];
                seed_set.erase(seed_set.begin() + best_rem_pos);
                // Update tracking
                auto it_u = upper_idx_map.find(rem_v);
                if (it_u != upper_idx_map.end()) {
                    in_cand[it_u->second] = false;
                    cand_upper_idx.erase(remove(cand_upper_idx.begin(), cand_upper_idx.end(), it_u->second), cand_upper_idx.end());
                }
                auto it_l = lower_idx_map.find(rem_v);
                if (it_l != lower_idx_map.end()) {
                    in_cand[ku + it_l->second] = false;
                    cand_lower_idx.erase(remove(cand_lower_idx.begin(), cand_lower_idx.end(), it_l->second), cand_lower_idx.end());
                }
                cur_size = seed_set.size();
                add_candidate(seed_set, "PairShrink-" + to_string(pi));
            }
        }
    }

    // =========================================================
    // Strategy C: Charikar edge density peel + refinement
    // =========================================================
    for (int pool_k : {min(n, 150), min(n, 300)}) {
        vector<int> pool;
        for (int i = 0; i < pool_k; i++) pool.push_back(deg_sorted[i].second);
        auto [peel_set, peel_d] = charikar_edge_density_peel(g, g_noisy, pool, p_flip, rng);
        if ((int)peel_set.size() >= P + Q) {
            add_candidate(peel_set, "Charikar-" + to_string(pool_k));
            // Sub-peel to smaller sizes
            for (int k : {P + Q, 5, 8, 10, 15, 20, 30}) {
                if (k >= (int)peel_set.size() || k < P + Q) continue;
                vector<int> sub = subsample_by_noisy_degree(g, g_noisy, peel_set, k, p_flip, rng);
                add_candidate(sub, "Charikar-sub" + to_string(k));
            }
        }
    }

    // =========================================================
    // Strategy D: Top-K degree seed sets at various sizes
    // =========================================================
    for (int k : {P + Q, 5, 8, 10, 15, 20, 30, 50}) {
        if (k > ku + kl) continue;
        // Take top-k/2 upper and top-k/2 lower by degree
        int nu_take = max(P, k / 2);
        int nl_take = max(Q, k - nu_take);
        nu_take = min(nu_take, ku);
        nl_take = min(nl_take, kl);
        vector<int> cand;
        for (int i = 0; i < nu_take; i++) cand.push_back(top_upper[i]);
        for (int i = 0; i < nl_take; i++) cand.push_back(top_lower[i]);
        add_candidate(cand, "TopK-" + to_string(k));
    }

    cout << "  [STATS] Evaluated " << ncand_eval << " candidates" << endl;

    // =========================================================
    // Report best result
    // =========================================================
    vector<int> S_upper, S_lower;
    for (int u : overall_best_set) {
        if (g.is_upper(u)) S_upper.push_back(u);
        else S_lower.push_back(u);
    }

    cout << "  [FINAL] |S|=" << overall_best_set.size()
         << " (u=" << S_upper.size() << ",l=" << S_lower.size() << ")"
         << " true_d=" << fixed << setprecision(6) << overall_best_d << endl;

    result.vertex_set = overall_best_set;
    result.estimated_density = overall_best_d;
    result.real_density = overall_best_d;
    return result;
}
