/*
 * biclique_dds_eldp_v17.cpp — v17 数据收集版本
 * 目标：收集详细数据来验证 OPTIMIZATION_ROADMAP 中的假设
 */

#include "biclique_dds_eldp_v17.h"
#include "biclique_dds_eldp.h"
#include "biclique.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <filesystem>

using namespace std;

extern DDSResult eldp_biclique_dds(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed);

// 计算 engagement
static unordered_map<int, long double> compute_engagement(BiGraph& g, const vector<int>& S, int P, int Q) {
    unordered_map<int, long double> engagement;
    if (S.empty()) return engagement;
    
    for (int v : S) engagement[v] = 0.0L;
    
    vector<int> S_upper, S_lower;
    for (int v : S) {
        if (g.is_upper(v)) S_upper.push_back(v);
        else S_lower.push_back(v);
    }
    
    if ((int)S_upper.size() < P || (int)S_lower.size() < Q) return engagement;
    
    // 采样方式计算 engagement
    int max_samples = 1000;  // 减少采样数量以加快速度
    int sample_count = 0;
    mt19937 rng(12345);
    
    while (sample_count < max_samples) {
        vector<int> upper_sample = S_upper;
        shuffle(upper_sample.begin(), upper_sample.end(), rng);
        if ((int)upper_sample.size() > P) upper_sample.resize(P);
        if ((int)upper_sample.size() < P) break;
        
        unordered_set<int> common;
        bool first = true;
        for (int u : upper_sample) {
            unordered_set<int> neighbors;
            for (int v : S_lower) {
                if (g.has(u, v)) neighbors.insert(v);
            }
            if (first) {
                common = neighbors;
                first = false;
            } else {
                unordered_set<int> new_common;
                for (int v : common) {
                    if (neighbors.count(v)) new_common.insert(v);
                }
                common = new_common;
            }
        }
        
        if ((int)common.size() >= Q) {
            vector<int> lower_sample(common.begin(), common.end());
            shuffle(lower_sample.begin(), lower_sample.end(), rng);
            if ((int)lower_sample.size() > Q) lower_sample.resize(Q);
            
            for (int v : upper_sample) engagement[v] += 1.0L;
            for (int v : lower_sample) engagement[v] += 1.0L;
            sample_count++;
        }
        if (sample_count >= max_samples) break;
    }
    return engagement;
}

static void save_engagement_data(const unordered_map<int, long double>& engagement, const string& filename) {
    ofstream out(filename);
    out << "vertex_id,engagement\n";
    for (const auto& p : engagement) {
        out << p.first << "," << fixed << setprecision(6) << p.second << "\n";
    }
    out.close();
}


static long double estimate_density_variance(BiGraph& g, const vector<int>& S, int P, int Q,
                                              long double epsilon, int num_samples = 10) {
    if (S.empty()) return 0.0L;
    
    vector<long double> estimates;
    mt19937 rng(54321);
    
    for (int i = 0; i < num_samples; i++) {
        long double true_d = count_bicliques_in_set(g, S, P, Q) / (long double)S.size();
        long double noise_scale = 1.0L / epsilon;
        uniform_real_distribution<double> unif(0.0, 1.0);
        double u = unif(rng) - 0.5;
        long double noise = noise_scale * (u > 0 ? log(1 + 2*u) : -log(1 - 2*u));
        estimates.push_back(true_d + noise);
    }
    
    long double mean = 0.0L;
    for (auto e : estimates) mean += e;
    mean /= estimates.size();
    
    long double variance = 0.0L;
    for (auto e : estimates) variance += (e - mean) * (e - mean);
    variance /= estimates.size();
    
    return variance;
}

static void save_variance_data(const vector<tuple<int, long double, long double>>& variance_data,
                                const string& filename) {
    ofstream out(filename);
    out << "set_size,true_density,estimated_variance\n";
    for (const auto& t : variance_data) {
        out << get<0>(t) << "," << fixed << setprecision(6) 
            << get<1>(t) << "," << get<2>(t) << "\n";
    }
    out.close();
}

DDSResultV17 eldp_v17(BiGraph& g, int P, int Q, long double epsilon, unsigned long seed,
                      const string& dataset_name, int round_id) {
    DDSResultV17 result;
    
    cout << "  [ELDP v17] Data Collection Mode" << endl;
    cout << "  Dataset: " << dataset_name << ", Round: " << round_id << endl;
    cout << "  P=" << P << ", Q=" << Q << ", eps=" << fixed << setprecision(4) << epsilon << endl;
    
    namespace fs = std::filesystem;
    fs::path output_dir("/data/gengdaz/biclique_extension/experiment_results/v17_analysis");
    fs::create_directories(output_dir);
    
    stringstream ss;
    ss << dataset_name << "_P" << P << "_Q" << Q << "_eps" << fixed << setprecision(1) << epsilon 
       << "_r" << round_id;
    string base_name = ss.str();
    
    result.trace_file = (output_dir / (base_name + "_trace.csv")).string();
    result.engagement_file = (output_dir / (base_name + "_engagement.csv")).string();
    result.variance_file = (output_dir / (base_name + "_variance.csv")).string();
    
    cout << "  [v17] Running v11 core with multi-start..." << endl;
    
    unsigned long s1 = seed;
    unsigned long s2 = seed ^ 0x9E3779B97F4A7C15ULL;
    unsigned long s3 = seed ^ 0xC2B2AE3D27D4EB4FULL;
    
    DDSResult r1 = eldp_biclique_dds(g, P, Q, epsilon, s1);
    DDSResult r2 = eldp_biclique_dds(g, P, Q, epsilon, s2);
    DDSResult r3 = eldp_biclique_dds(g, P, Q, epsilon, s3);
    
    const DDSResult* best = &r1;
    if (r2.estimated_density > best->estimated_density) best = &r2;
    if (r3.estimated_density > best->estimated_density) best = &r3;
    
    result.vertex_set = best->vertex_set;
    result.estimated_density = best->estimated_density;
    result.real_density = best->real_density;
    
    cout << "  [v17] Best result: |S|=" << result.vertex_set.size() 
         << ", density=" << result.real_density << endl;
    
    cout << "  [v17] Computing engagement..." << endl;
    auto engagement = compute_engagement(g, result.vertex_set, P, Q);
    save_engagement_data(engagement, result.engagement_file);
    
    vector<long double> eng_values;
    for (const auto& p : engagement) eng_values.push_back(p.second);
    sort(eng_values.begin(), eng_values.end());
    
    if (!eng_values.empty()) {
        long double eng_mean = 0.0L;
        for (auto e : eng_values) eng_mean += e;
        eng_mean /= eng_values.size();
        long double eng_median = eng_values[eng_values.size() / 2];
        long double eng_p90 = eng_values[(size_t)(eng_values.size() * 0.9)];
        cout << "  [v17] Engagement: mean=" << fixed << setprecision(2) << eng_mean
             << ", median=" << eng_median << ", p90=" << eng_p90 << endl;
    }
    
    cout << "  [v17] Analyzing variance vs |S|..." << endl;
    vector<tuple<int, long double, long double>> variance_data;
    
    if (result.vertex_set.size() > 10) {
        vector<int> sizes = {10, 20, 50, 100, 200, 500, 1000};
        for (int sz : sizes) {
            if (sz > (int)result.vertex_set.size()) break;
            vector<int> subset(result.vertex_set.begin(), result.vertex_set.begin() + sz);
            long double true_d = count_bicliques_in_set(g, subset, P, Q) / (long double)sz;
            long double var = estimate_density_variance(g, subset, P, Q, epsilon, 10);
            variance_data.push_back(make_tuple(sz, true_d, var));
        }
        save_variance_data(variance_data, result.variance_file);
    }
    
    ofstream trace_out(result.trace_file);
    trace_out << "set_size,estimated_density,real_density,label\n";
    trace_out << result.vertex_set.size() << "," 
              << fixed << setprecision(6) << result.estimated_density << ","
              << result.real_density << ",final\n";
    trace_out.close();
    
    cout << "  [v17] Data saved to " << output_dir << endl;
    
    return result;
}
