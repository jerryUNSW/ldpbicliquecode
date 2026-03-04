/*
 * exact_dds.cpp
 * Exact (p,q)-biclique densest subgraph solver via greedy peeling.
 *
 * Usage: ./exact_dds <dataset> <P> <Q>
 * Output: writes to ../Exact_result/{dataset}_{P}_{Q}_exact.txt
 */

#include "biclique.h"
#include <fstream>
#include <sys/stat.h>

using namespace std;

// Globals DEFINED here (biclique.cpp externs them from main)
long double Eps = 0, Eps0 = 0, Eps1 = 0, Eps2 = 0, p = 0;
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0;
long double communication_cost = 0;
long double alpha0 = 0.05, alpha1 = 0.6, alpha2 = 0.35;
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

// Globals DEFINED in biclique.cpp — extern here
extern long double _cate, _wedge, _btf;
extern bool one_round, count_cate, sampling_noisy_graph, eva_comm;
extern bool edge_clipping, count_cc;
extern double p____;
extern int alpha;
extern stats::rand_engine_t engine;
extern bool two_noisy_graph_switch, multi_estimator_switch;
// These are extern in biclique.h, defined in main.cpp — define here
bool use_probability_filtering = false;
int post_processing_mode = 0;
long double winsorize_percentile = 95.0;
long double truncation_bound = 0.0;
extern long double gamma__;

/*
 * Count (p,q)-bicliques in the induced subgraph on vertex set S.
 */
static long double count_bicliques_in_subgraph(BiGraph& g, const vector<int>& S, int P, int Q) {
    vector<int> upper_verts, lower_verts;
    for (int v : S) {
        if (g.is_upper(v)) upper_verts.push_back(v);
        else lower_verts.push_back(v);
    }
    if ((int)upper_verts.size() < P || (int)lower_verts.size() < Q) return 0.0;
    
    unordered_map<int, int> upper_map, lower_map;
    for (int i = 0; i < (int)upper_verts.size(); i++) upper_map[upper_verts[i]] = i;
    for (int i = 0; i < (int)lower_verts.size(); i++) lower_map[lower_verts[i]] = i;
    
    BiGraph sub;
    sub.init(upper_verts.size(), lower_verts.size());
    for (int u : upper_verts) {
        for (int nb : g.neighbor[u]) {
            if (lower_map.find(nb) != lower_map.end())
                sub.addEdge(upper_map[u], lower_map[nb] + upper_verts.size());
        }
    }
    if (sub.num_edges == 0) return 0.0;
    
    biGraph bg = convertBiGraphTobiGraph(sub);
    BCListPlusPlus counter(&bg, P, Q);
    return counter.exactCount();
}

struct ExactResult {
    vector<int> vertex_set;
    long double density;
    long double biclique_count;
};

static ExactResult greedy_peeling_biclique_dds(BiGraph& g, int P, int Q) {
    int n = g.num_nodes();
    vector<int> current_set;
    for (int i = 0; i < n; i++) current_set.push_back(i);
    
    long double total_bicliques = count_bicliques_in_subgraph(g, current_set, P, Q);
    long double best_density = (current_set.size() > 0) ? total_bicliques / (long double)current_set.size() : 0;
    vector<int> best_set = current_set;
    long double best_count = total_bicliques;
    
    cout << "Initial: |S|=" << current_set.size()
         << ", bicliques=" << total_bicliques
         << ", density=" << best_density << endl;
    
    while ((int)current_set.size() > max(P, Q)) {
        int worst_vertex = -1;
        long double min_participation = 1e30;
        
        for (int v : current_set) {
            vector<int> set_without_v;
            for (int u : current_set) if (u != v) set_without_v.push_back(u);
            long double bicliques_without = count_bicliques_in_subgraph(g, set_without_v, P, Q);
            long double participation = total_bicliques - bicliques_without;
            if (participation < min_participation) {
                min_participation = participation;
                worst_vertex = v;
            }
        }
        if (worst_vertex == -1) break;
        
        current_set.erase(remove(current_set.begin(), current_set.end(), worst_vertex), current_set.end());
        total_bicliques -= min_participation;
        
        if (current_set.size() % 50 == 0)
            total_bicliques = count_bicliques_in_subgraph(g, current_set, P, Q);
        
        long double current_density = (current_set.size() > 0) ?
            total_bicliques / (long double)current_set.size() : 0;
        
        if (current_density > best_density) {
            best_density = current_density;
            best_set = current_set;
            best_count = total_bicliques;
        }
        if (current_set.size() % 100 == 0 || current_set.size() <= 20)
            cout << "  Peeling: |S|=" << current_set.size()
                 << ", density=" << current_density
                 << ", best=" << best_density << endl;
    }
    
    best_count = count_bicliques_in_subgraph(g, best_set, P, Q);
    best_density = best_count / (long double)best_set.size();
    return {best_set, best_density, best_count};
}

static void save_exact_result(const string& dataset, int P, int Q, const ExactResult& result) {
    size_t found = dataset.find_last_of("/\\");
    string ds_name = (found != string::npos) ? dataset.substr(found + 1) : dataset;
    
    string filename = "../Exact_result/" + ds_name + "_" + to_string(P) + "_" + to_string(Q) + "_exact.txt";
    ofstream fout(filename);
    fout << "dataset: " << ds_name << "\nP: " << P << "\nQ: " << Q
         << "\nset_size: " << result.vertex_set.size()
         << fixed << setprecision(10)
         << "\nbiclique_count: " << result.biclique_count
         << "\ndensity: " << result.density << "\nvertex_set:";
    for (int v : result.vertex_set) fout << " " << v;
    fout << endl;
    fout.close();
    cout << "Saved to " << filename << endl;
    
    sqlite3* db;
    string db_path = "../Exact_result/biclique_dds_exact.db";
    if (sqlite3_open(db_path.c_str(), &db) != SQLITE_OK) {
        cerr << "Error opening database: " << sqlite3_errmsg(db) << endl;
        return;
    }
    const char* create_sql =
        "CREATE TABLE IF NOT EXISTS exact_results ("
        "dataset TEXT, p INTEGER, q INTEGER, "
        "set_size INTEGER, density REAL, biclique_count REAL, "
        "PRIMARY KEY(dataset, p, q));";
    sqlite3_exec(db, create_sql, nullptr, nullptr, nullptr);
    
    const char* insert_sql =
        "INSERT OR REPLACE INTO exact_results (dataset, p, q, set_size, density, biclique_count) "
        "VALUES (?, ?, ?, ?, ?, ?);";
    sqlite3_stmt* stmt;
    if (sqlite3_prepare_v2(db, insert_sql, -1, &stmt, nullptr) == SQLITE_OK) {
        sqlite3_bind_text(stmt, 1, ds_name.c_str(), -1, SQLITE_STATIC);
        sqlite3_bind_int(stmt, 2, P);
        sqlite3_bind_int(stmt, 3, Q);
        sqlite3_bind_int(stmt, 4, (int)result.vertex_set.size());
        sqlite3_bind_double(stmt, 5, (double)result.density);
        sqlite3_bind_double(stmt, 6, (double)result.biclique_count);
        sqlite3_step(stmt);
        sqlite3_finalize(stmt);
    }
    sqlite3_close(db);
    cout << "Saved to " << db_path << endl;
}

static bool load_exact_result(const string& dataset, int P, int Q, long double& density) {
    size_t found = dataset.find_last_of("/\\");
    string ds_name = (found != string::npos) ? dataset.substr(found + 1) : dataset;
    
    sqlite3* db;
    string db_path = "../Exact_result/biclique_dds_exact.db";
    if (sqlite3_open(db_path.c_str(), &db) != SQLITE_OK) return false;
    
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
    return found_result;
}

static string exact_result_txt_path(const string& dataset, int P, int Q) {
    size_t found = dataset.find_last_of("/\\");
    string ds_name = (found != string::npos) ? dataset.substr(found + 1) : dataset;
    return "../Exact_result/" + ds_name + "_" + to_string(P) + "_" + to_string(Q) + "_exact.txt";
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: ./exact_dds <dataset> <P> <Q>" << endl;
        return 1;
    }
    string dataset = argv[1];
    P___ = atoi(argv[2]);
    K___ = atoi(argv[3]);
    
    cout << "=== Exact (p,q)-Biclique DDS Solver ===" << endl;
    cout << "Dataset: " << dataset << ", P=" << P___ << ", Q=" << K___ << endl;
    
    long double existing_density;
    if (load_exact_result(dataset, P___, K___, existing_density)) {
        string txt_path = exact_result_txt_path(dataset, P___, K___);
        ifstream fin(txt_path);
        if (fin.good()) {
            cout << "Already computed! Density = " << existing_density << endl;
            return 0;
        }
        cout << "Density exists in DB but exact set file missing. Recomputing to recover S*..." << endl;
    }
    
    BiGraph g(dataset);
    cout << "Graph: n1=" << g.num_v1 << ", n2=" << g.num_v2 << ", m=" << g.num_edges << endl;
    
    double t0 = omp_get_wtime();
    ExactResult result = greedy_peeling_biclique_dds(g, P___, K___);
    double t1 = omp_get_wtime();
    
    cout << "\n=== Result ===" << endl;
    cout << "Set size: " << result.vertex_set.size() << endl;
    cout << "Bicliques: " << fixed << setprecision(2) << result.biclique_count << endl;
    cout << "Density: " << fixed << setprecision(10) << result.density << endl;
    cout << "Time: " << fixed << setprecision(2) << (t1 - t0) << "s" << endl;
    
    save_exact_result(dataset, P___, K___, result);
    return 0;
}
