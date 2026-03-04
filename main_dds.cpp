/*
 * main_dds.cpp — Entry point for (p,q)-Biclique DDS under ELDP experiments.
 *
 * Usage: ./biclique_dds <epsilon> <dataset> <num_rounds> <P> <Q>
 */

#include "biclique_dds_eldp.h"
#include "biclique_dds_eldp_v12.h"
#include "biclique_dds_eldp_v13.h"
#include "biclique_dds_eldp_v14.h"
#include "biclique_dds_eldp_v16.h"
#include "biclique_dds_eldp_v17.h"
#include "biclique_dds_naive_eldp_v2.h"
#include "biclique_dds_peeling_eldp.h"
#include "biclique_dds_peeling_expand_shrink_eldp.h"
#include "biclique.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
#include <cmath>
#include <unordered_set>
#include <algorithm>
#include <omp.h>

using namespace std;

// Globals required by biclique.cpp (define here, biclique.cpp externs them)
long double Eps = 0, Eps0 = 0, Eps1 = 0, Eps2 = 0, p = 0;
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0;
long double communication_cost = 0;
long double alpha0 = 0.05, alpha1 = 0.95, alpha2 = 0.0;
vector<int> priv_deg;
vector<long double> naive_estis;
int priv_dmax_1 = 0, priv_dmax_2 = 0;
vector<vector<int>> up_options, lo_options;
int iteration = 0;
int K___, P___;
long double real = 0.0, real_ld = 0.0;
long double avg_estimated_variance = 0;
long double RR_time = 0, server_side_time = 0, naive_server_side = 0,
            local_count_time = 0, deg_esti_time = 0;
int num_rounds = 1;
vector<long double> estis;

// Globals defined in biclique.cpp — extern here
extern long double _cate, _wedge, _btf;
extern bool one_round, count_cate, sampling_noisy_graph, eva_comm;
extern bool edge_clipping, count_cc;
extern double p____;
extern int alpha;
extern stats::rand_engine_t engine;
extern bool two_noisy_graph_switch, multi_estimator_switch;
extern long double gamma__;

// These are extern in biclique.h, defined here
bool use_probability_filtering = false;
int post_processing_mode = 0;
long double winsorize_percentile = 95.0;
long double truncation_bound = 0.0;

struct EngagementSummary {
    long double mean = 0.0L;
    long double median = 0.0L;
    long double p90 = 0.0L;
    long double min_v = 0.0L;
    long double max_v = 0.0L;
};

static bool load_exact_vertex_set(const string& dataset, int P, int Q, vector<int>& out_set) {
    size_t found = dataset.find_last_of("/\\");
    string ds_name = (found != string::npos) ? dataset.substr(found + 1) : dataset;
    vector<string> candidates = {
        "/data/gengdaz/biclique_extension/Exact_result/" + ds_name + "_" + to_string(P) + "_" + to_string(Q) + "_exact.txt",
        "../Exact_result/" + ds_name + "_" + to_string(P) + "_" + to_string(Q) + "_exact.txt"
    };
    for (const auto& path : candidates) {
        ifstream fin(path);
        if (!fin.good()) continue;
        string line;
        while (getline(fin, line)) {
            if (line.rfind("vertex_set:", 0) == 0) {
                out_set.clear();
                istringstream iss(line.substr(string("vertex_set:").size()));
                int v;
                while (iss >> v) out_set.push_back(v);
                return !out_set.empty();
            }
        }
    }
    return false;
}

static vector<long double> compute_biclique_engagement(BiGraph& g, const vector<int>& S, int P, int Q) {
    vector<long double> e(S.size(), 0.0L);
    if (S.empty()) return e;
    long double total = count_bicliques_in_set(g, S, P, Q);
    vector<int> sub;
    sub.reserve(S.size() - 1);
    for (size_t i = 0; i < S.size(); i++) {
        sub.clear();
        for (size_t j = 0; j < S.size(); j++) {
            if (j == i) continue;
            sub.push_back(S[j]);
        }
        long double without = count_bicliques_in_set(g, sub, P, Q);
        long double part = total - without;
        e[i] = (part < 0 ? 0.0L : part);
    }
    return e;
}

static EngagementSummary summarize_engagement(const vector<long double>& e) {
    EngagementSummary s;
    if (e.empty()) return s;
    vector<long double> v = e;
    sort(v.begin(), v.end());
    long double sum = 0.0L;
    for (auto x : v) sum += x;
    s.mean = sum / (long double)v.size();
    s.median = v[v.size() / 2];
    s.p90 = v[(size_t)floor((long double)(v.size() - 1) * 0.90L)];
    s.min_v = v.front();
    s.max_v = v.back();
    return s;
}

static void dump_engagement_csv(
    const string& dataset,
    const string& mode_name,
    int P, int Q,
    int round_idx,
    const vector<int>& S,
    const vector<long double>& e,
    const string& tag) {

    namespace fs = std::filesystem;
    fs::path out_dir("/data/gengdaz/biclique_extension/experiment_results/engagement");
    fs::create_directories(out_dir);
    size_t pos = dataset.find_last_of("/\\");
    string ds_name = (pos == string::npos) ? dataset : dataset.substr(pos + 1);
    fs::path out = out_dir / (mode_name + "_" + ds_name + "_p" + to_string(P) + "_q" + to_string(Q) +
                              "_r" + to_string(round_idx) + "_" + tag + ".csv");
    ofstream fout(out.string());
    fout << "vertex_id,engagement\n";
    for (size_t i = 0; i < S.size() && i < e.size(); i++) {
        fout << S[i] << "," << fixed << setprecision(10) << e[i] << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        cerr << "Usage: ./biclique_dds <epsilon> <dataset> <num_rounds> <P> <Q> <mode>" << endl;
        cerr << "Mode: 0=edlp_v12, 1=naive_eldp_v2, 2=edlp_v13, 3=edlp_v14, 4=edlp_v16, 5=edlp_v17, 6=peeling_eldp, 7=peeling_expand_shrink_eldp" << endl;
        return 1;
    }

    long double epsilon = atof(argv[1]);
    string dataset = argv[2];
    int rounds = atoi(argv[3]);
    P___ = atoi(argv[4]);
    K___ = atoi(argv[5]);
    int mode = atoi(argv[6]);

    string mode_name;
    if (mode == 0) mode_name = "edlp_v12";
    else if (mode == 1) mode_name = "naive_eldp_v2";
    else if (mode == 2) mode_name = "edlp_v13";
    else if (mode == 3) mode_name = "edlp_v14";
    else if (mode == 4) mode_name = "edlp_v16";
    else if (mode == 5) mode_name = "edlp_v17";
    else if (mode == 6) mode_name = "peeling_eldp";
    else if (mode == 7) mode_name = "peeling_expand_shrink_eldp";
    else {
        cerr << "Unsupported mode: " << mode << endl;
        cerr << "Mode: 0=edlp_v12, 1=naive_eldp_v2, 2=edlp_v13, 3=edlp_v14, 4=edlp_v16, 5=edlp_v17, 6=peeling_eldp, 7=peeling_expand_shrink_eldp" << endl;
        return 1;
    }

    cout << "=== (p,q)-Biclique DDS under ELDP ===" << endl;
    cout << "Dataset: " << dataset << " | Mode: " << mode_name << endl;
    cout << "P=" << P___ << ", Q=" << K___ << ", Epsilon=" << epsilon << " | Rounds=" << rounds << endl;

    // Load graph
    BiGraph g(dataset);
    cout << "Graph: n1=" << g.num_v1 << ", n2=" << g.num_v2 << ", m=" << g.num_edges << endl;

    // Load exact optimal density
    long double exact_density = 0;
    bool has_exact = load_exact_density(dataset, P___, K___, exact_density);
    if (has_exact) {
        cout << "Exact optimal density: " << fixed << setprecision(10) << exact_density << endl;
    } else {
        cout << "WARNING: No exact density found for this dataset/P/Q" << endl;
    }

    // Run multiple rounds
    random_device rd;
    uint64_t tnow = (uint64_t)chrono::high_resolution_clock::now().time_since_epoch().count();
    seed_seq seq{
        rd(), rd(), rd(), rd(),
        (uint32_t)(tnow & 0xffffffffULL), (uint32_t)((tnow >> 32) & 0xffffffffULL),
        (uint32_t)hash<string>{}(dataset), (uint32_t)P___, (uint32_t)K___, (uint32_t)mode
    };
    mt19937_64 seed_gen(seq);
    unordered_set<unsigned long> used_seeds;

    vector<long double> real_densities;
    vector<long double> density_ratios;
    vector<int> exact_set;
    bool has_exact_set = load_exact_vertex_set(dataset, P___, K___, exact_set);
    if (has_exact_set) {
        cout << "Exact set loaded: |S*|=" << exact_set.size() << endl;
    } else {
        cout << "WARNING: Exact set file not found, skip S vs S* and engagement compare." << endl;
    }

    double t0 = omp_get_wtime();

    for (int round = 0; round < rounds; round++) {
        cout << "\n--- Round " << (round + 1) << "/" << rounds << " ---" << endl;
        unsigned long seed = seed_gen();
        while (used_seeds.find(seed) != used_seeds.end()) seed = seed_gen();
        used_seeds.insert(seed);
        cout << "Seed: " << seed << endl;

        // Setup trace file for this round
        namespace fs = std::filesystem;
        fs::path trace_dir("/data/gengdaz/biclique_extension/experiment_results/Special_test");
        fs::create_directories(trace_dir);
        size_t pos = dataset.find_last_of("/\\");
        string ds_name = (pos == string::npos) ? dataset : dataset.substr(pos + 1);
        fs::path trace_path = trace_dir / ("greedy_trace_" + mode_name + "_" + ds_name + 
                                           "_p" + to_string(P___) + "_q" + to_string(K___) + 
                                           "_r" + to_string(round + 1) + ".csv");
        ofstream trace_file(trace_path.string());
        trace_file << "set_size,estimated_density,real_density,label\n";
        set_trace_file(&trace_file);

        DDSResult result;
        if (mode == 0) {
            DDSResultV12 r12 = eldp_v12(g, P___, K___, epsilon, seed);
            result.vertex_set = r12.vertex_set;
            result.real_density = r12.real_density;
            result.estimated_density = r12.estimated_density;
        } else if (mode == 1) {
            NaiveV2Result nr = naive_eldp_v2(g, P___, K___, epsilon, seed);
            result.vertex_set = nr.vertex_set;
            result.real_density = nr.real_density;
            result.estimated_density = nr.estimated_density;
        } else if (mode == 2) {
            DDSResultV13 r13 = eldp_v13(g, P___, K___, epsilon, seed);
            result.vertex_set = r13.vertex_set;
            result.real_density = r13.real_density;
            result.estimated_density = r13.estimated_density;
        } else if (mode == 3) {
            DDSResultV14 r14 = eldp_v14(g, P___, K___, epsilon, seed);
            result.vertex_set = r14.vertex_set;
            result.real_density = r14.real_density;
            result.estimated_density = r14.estimated_density;
        } else if (mode == 4) {
            // v16 with default parameters (can be tuned)
            int size_threshold = 100;
            long double variance_threshold = 0.5L;
            long double precision_threshold = 0.5L;
            int max_size = 500;
            int top_k = 3;
            
            // Read from environment variables if set
            const char* env_size_th = getenv("V16_SIZE_THRESHOLD");
            const char* env_var_th = getenv("V16_VAR_THRESHOLD");
            const char* env_prec_th = getenv("V16_PREC_THRESHOLD");
            const char* env_max_sz = getenv("V16_MAX_SIZE");
            const char* env_top_k = getenv("V16_TOP_K");
            
            if (env_size_th) size_threshold = atoi(env_size_th);
            if (env_var_th) variance_threshold = atof(env_var_th);
            if (env_prec_th) precision_threshold = atof(env_prec_th);
            if (env_max_sz) max_size = atoi(env_max_sz);
            if (env_top_k) top_k = atoi(env_top_k);
            
            DDSResultV16 r16 = eldp_v16(g, P___, K___, epsilon, seed,
                                        size_threshold, variance_threshold,
                                        precision_threshold, max_size, top_k);
            result.vertex_set = r16.vertex_set;
            result.real_density = r16.real_density;
            result.estimated_density = r16.estimated_density;
        } else if (mode == 5) {
            // v17 data collection
            size_t pos = dataset.find_last_of("/\\");
            string ds_name = (pos == string::npos) ? dataset : dataset.substr(pos + 1);
            
            DDSResultV17 r17 = eldp_v17(g, P___, K___, epsilon, seed, ds_name, round + 1);
            result.vertex_set = r17.vertex_set;
            result.real_density = r17.real_density;
            result.estimated_density = r17.estimated_density;
        } else if (mode == 6) {
            DDSResult r6 = peeling_eldp(g, P___, K___, epsilon, seed);
            result.vertex_set = r6.vertex_set;
            result.real_density = r6.real_density;
            result.estimated_density = r6.estimated_density;
        } else {
            // mode == 7: peeling_expand_shrink_eldp
            DDSResult r7 = peeling_expand_shrink_eldp(g, P___, K___, epsilon, seed);
            result.vertex_set = r7.vertex_set;
            result.real_density = r7.real_density;
            result.estimated_density = r7.estimated_density;
        }

        // Close trace file
        set_trace_file(nullptr);
        trace_file.close();

        real_densities.push_back(result.real_density);
        cout << "  Real density: " << fixed << setprecision(6) << result.real_density << endl;
        cout << "  Set size: " << result.vertex_set.size() << endl;
        if (has_exact && exact_density > 0) {
            long double ratio = result.real_density / exact_density;
            density_ratios.push_back(ratio);
            cout << "  Density ratio: " << fixed << setprecision(6) << ratio << endl;
        }

        if (has_exact_set) {
            unordered_set<int> exact_idx(exact_set.begin(), exact_set.end());
            int inter = 0;
            for (int v : result.vertex_set) if (exact_idx.count(v)) inter++;
            long double recall_s_star = exact_set.empty() ? 0.0L : (long double)inter / (long double)exact_set.size();
            long double precision_on_s = result.vertex_set.empty() ? 0.0L : (long double)inter / (long double)result.vertex_set.size();
            cout << "  Overlap with S*: |S∩S*|=" << inter
                 << " recall(|S*|)=" << fixed << setprecision(6) << recall_s_star
                 << " precision(|S|)=" << precision_on_s << endl;

            vector<long double> e_s = compute_biclique_engagement(g, result.vertex_set, P___, K___);
            vector<long double> e_star = compute_biclique_engagement(g, exact_set, P___, K___);
            EngagementSummary ss = summarize_engagement(e_s);
            EngagementSummary st = summarize_engagement(e_star);

            cout << "  Engagement S    mean=" << fixed << setprecision(3) << ss.mean
                 << " median=" << ss.median << " p90=" << ss.p90
                 << " min=" << ss.min_v << " max=" << ss.max_v << endl;
            cout << "  Engagement S*   mean=" << fixed << setprecision(3) << st.mean
                 << " median=" << st.median << " p90=" << st.p90
                 << " min=" << st.min_v << " max=" << st.max_v << endl;

            dump_engagement_csv(dataset, mode_name, P___, K___, round + 1, result.vertex_set, e_s, "S");
            dump_engagement_csv(dataset, mode_name, P___, K___, round + 1, exact_set, e_star, "Sstar");
        }
    }

    double t1 = omp_get_wtime();

    cout << "\n=== Summary ===" << endl;
    cout << "Time: " << fixed << setprecision(2) << (t1 - t0) << "s" << endl;

    if (!real_densities.empty()) {
        long double mean_d = 0;
        for (auto d : real_densities) mean_d += d;
        mean_d /= real_densities.size();
        cout << "Mean real density: " << fixed << setprecision(6) << mean_d << endl;
    }
    if (!density_ratios.empty()) {
        long double mean_r = 0;
        for (auto r : density_ratios) mean_r += r;
        mean_r /= density_ratios.size();
        long double var_r = 0;
        for (auto r : density_ratios) var_r += (r - mean_r) * (r - mean_r);
        var_r /= density_ratios.size();
        cout << "Mean density ratio: " << fixed << setprecision(6) << mean_r << endl;
        cout << "Std density ratio: " << fixed << setprecision(6) << sqrt((double)var_r) << endl;
        cout << "dds_ratio = " << fixed << setprecision(6) << mean_r << endl;
    }

    return 0;
}
