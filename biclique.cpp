#include "biclique.h"
#include "include/mt19937ar.h"
#include <unordered_set>
#include <fstream>

#include <sys/resource.h>
void printMemoryUsage() {

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "Memory usage: " << usage.ru_maxrss << " KB" << std::endl;
}

using namespace std;
long double _cate, _wedge, _btf;
vector<int> upper_sample, lower_sample;
double sample_ratio = 1.0;
unordered_map<int, bool> in_sampled_up, in_sampled_lo;
bool samling_one_round = false;
long double verified = 0, not_verified = 0;
extern long double Eps, Eps0, Eps1, Eps2, p, m3__, m2__, m1__, m0__, communication_cost;
extern long double alpha0, alpha1, alpha2;
extern vector<int> priv_deg;
extern vector<long double> naive_estis;
extern int priv_dmax_1, priv_dmax_2;
extern vector<vector<int>> up_options, lo_options;
extern int iteration;

// double noisy graph optimization?
bool two_noisy_graph_switch ; 
bool multi_estimator_switch ; 

extern int K___ ; 
extern int P___;

extern long double real; 
extern long double real_ld;  // High precision ground truth for large numbers 


extern long double avg_estimated_variance ;

bool one_round = false, 
    // multiR= false, wedge_based = false, 
    edge_clipping = true, count_cate = false, 
	sampling_noisy_graph = false, eva_comm = false;

bool count_cc = false ;  // this is computing the bipartite clustering coefficient. 
double p____ = 0.5;  // sampling ratio.
int alpha = 10; 
stats::rand_engine_t engine(std::time(0));  // used to be 1776 

bool vertex_pair_reduction = true; 

// bool averaging_f_estimates = true;
bool averaging_f_estimates = true;  

extern vector<vector<int>> up_options, lo_options; 


extern long double RR_time, server_side_time, naive_server_side, 
	local_count_time, deg_esti_time;


long double gamma__;

bool efficient_RR = true, skip_neg_deg = false ;
// why cannot we skip negative vertices

// convert your BiGraph instance g2 to the biGraph struct 
biGraph convertBiGraphTobiGraph(BiGraph& oldGraph) {
    biGraph newGraph;
    
    // Initialize basic properties
    newGraph.n1 = oldGraph.num_v1;
    newGraph.n2 = oldGraph.num_v2;
    newGraph.m = oldGraph.num_edges;


    // Allocate space
    newGraph.edges.resize(newGraph.m);
    newGraph.e1.resize(newGraph.m);
    newGraph.e2.resize(newGraph.m);
    newGraph.pU.resize(newGraph.n1 + 5);
    newGraph.pV.resize(newGraph.n2 + 5);
    
    int edge_index = 0;
    for (vid_t u = 0; u < oldGraph.num_v1; ++u) {
        // for each upper vertex.
        assert(oldGraph.is_upper(u));

        // Make a copy of neighbors to avoid destroying the original graph
        auto neighbors = oldGraph.neighbor[u];
        for (const auto& v__ : neighbors) {
            assert(oldGraph.is_lower(v__));
            vid_t v = v__ - oldGraph.num_v1;
            assert(v < oldGraph.num_v2); 

            newGraph.edges[edge_index].u = u;         
            newGraph.edges[edge_index].v = v;          
            edge_index++;
        }
        // Note: We don't clear neighbors anymore to preserve original graph
    }

    // cout<<"edge_index = "<<edge_index<<endl;
    // Note: edge_index might be less than newGraph.m due to duplicate edges in original graph
    if (edge_index != newGraph.m) {
        std::cout << "Warning: Edge count mismatch - expected: " << newGraph.m 
                  << ", actual: " << edge_index 
                  << " (likely due to duplicate edges in original graph)" << std::endl;
        newGraph.m = edge_index;  // Update to actual edge count
    }

    // compute the degree ordering in new graph
    newGraph.changeToDegreeOrder();

    return newGraph;
}

void randomized_response_single_bit(int u, int v, BiGraph& g, BiGraph& g2) {
	double keep_probability = g.has(u, v) ? 1 - p : p;
	// if(sampling_noisy_graph){
	// 	keep_probability *= p____;
	// }
	assert(u!=v);
	g2.edge_vector[min(u, v)][max(u, v)] = (genrand_real2() < keep_probability);
	// it is either 1 or 0
}

void private_estimate_of_degrees(BiGraph& g) {
    // private estimate degrees.
    priv_dmax_1 = 0;
    priv_dmax_2 = 0;
    priv_deg.resize(g.num_nodes());
    for (int i = 0; i < g.num_nodes(); i++) {
        priv_deg[i] = g.degree[i] + stats::rlaplace(0.0, 1 / (Eps0), engine);
        if (edge_clipping) {
            priv_deg[i] += alpha;
        }
        if (g.is_upper(i)) {
            priv_dmax_1 = priv_dmax_1 > priv_deg[i] ? priv_dmax_1 : priv_deg[i];
        } else {
            priv_dmax_2 = priv_dmax_2 > priv_deg[i] ? priv_dmax_2 : priv_deg[i];
        }
    }
}

double my_genrand_real2() { return genrand_real2(); }

void my_init_genrand(unsigned long seed) { init_genrand(seed); }

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed) {
    const int range_from = g2.num_v1;
    const int range_to = g2.num_nodes() - 1;
    std::mt19937 generator(seed);  // Use the provided seed instead of random_device
    std::uniform_int_distribution<int> distr(range_from, range_to);

    init_genrand(seed);
    int flip1 = 0;
    int visited_vertices = 0;
    int total_vertices = g2.num_v1;
    int ten_percent = total_vertices / 5;
    long double max_time_per_user = -1;

    int max_num_noisy_edges_per_vertex = -1;
    for (int i = 0; i < g2.num_v1; i++) {
        visited_vertices++;
        // if (visited_vertices % ten_percent == 0) {
        // 	int progress = visited_vertices * 100 / total_vertices;
        // 	cout << "Processed " << progress << "% of vertices" << endl;
        // }
        int num_noisy_edges_per_i = 0;
        double tx = omp_get_wtime();
        for (int j = g2.num_v1; j < g2.num_nodes(); j++) {
            if (std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) !=
                g.neighbor[i].end()) {
                if (genrand_real2() >= p) {  // 1  --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            } else {
                if (genrand_real2() < p) {  // 0 --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            }
        }
        if(num_noisy_edges_per_i > max_num_noisy_edges_per_vertex){
            max_num_noisy_edges_per_vertex =  num_noisy_edges_per_i; 
        }

        double ty = omp_get_wtime();
        max_time_per_user =max_time_per_user > (ty - tx) ? max_time_per_user : (ty - tx);
    }
    // RR_time += max_time_per_user;
    // the dominating cost is incurred in server side butterfly counting on the
    // dense noisy graph

    // cout << "noisy edges = " << flip1 << endl;

    if (eva_comm) {
        cout<<"computing the of Randomized responses (1)"<<endl;
        double byte_per_edge = 8*(log2(g.num_v1) + log2(g.num_v2));
        communication_cost += max_num_noisy_edges_per_vertex * byte_per_edge;
        // communication_cost += flip1 * sizeof(int);
    }

    long double expected_E =g.num_edges * (1 - p) + (g.num_v1 * g.num_v2 - g.num_edges) * p;
    // cout << "expected E = " << expected_E << endl;

    g2.computePriority();
}

// applying RR to lower vertices 
void construct_noisy_graph_2(BiGraph& g, BiGraph& g2, unsigned long seed) {
    const int range_from = g2.num_v1;
    const int range_to = g2.num_nodes() - 1;
    std::mt19937 generator(seed);  // Use the provided seed instead of random_device
    std::uniform_int_distribution<int> distr(range_from, range_to);

    init_genrand(seed);
    int flip1 = 0;
    int visited_vertices = 0;
    int total_vertices = g2.num_v1;
    int ten_percent = total_vertices / 5;
    long double max_time_per_user = -1;

    // applying RR on U(G)
    // for (int i = 0; i < g2.num_v1; i++) {

    // now applying RR to L(G)
    int max_num_noisy_edges_per_vertex = -1;
    for (int i = g2.num_v1; i < g2.num_nodes(); i++) {
        visited_vertices++;
        double tx = omp_get_wtime();

        // for (int j = g2.num_v1; j < g2.num_nodes(); j++) {
        int num_noisy_edges_per_i = 0;
        for (int j = 0; j < g2.num_v1; j++) {
            if (std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) !=
                g.neighbor[i].end()) {
                if (genrand_real2() >= p) {  // 1  --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            } else {
                if (genrand_real2() < p) {  // 0 --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            }
        }

        if(num_noisy_edges_per_i > max_num_noisy_edges_per_vertex){
            max_num_noisy_edges_per_vertex =  num_noisy_edges_per_i; 
        }
        double ty = omp_get_wtime();
        max_time_per_user =
            max_time_per_user > (ty - tx) ? max_time_per_user : (ty - tx);
    }
    // RR_time += max_time_per_user;
    // the dominating cost is incurred in server side butterfly counting on the
    // dense noisy graph

    // cout << "noisy edges = " << flip1 << endl;

    // if (eva_comm) {
    //     cout<<"computing the of Randomized responses (1)"<<endl;
    //     communication_cost += flip1 * sizeof(int);
    // }
    if (eva_comm) {
        cout<<"computing the of Randomized responses (1)"<<endl;
        double byte_per_edge = 8*(log2(g.num_v1) + log2(g.num_v2));
        communication_cost += max_num_noisy_edges_per_vertex * byte_per_edge;
        // communication_cost += flip1 * sizeof(int);
    }
    long double expected_E =g.num_edges * (1 - p) + (g.num_v1 * g.num_v2 - g.num_edges) * p;
    // cout << "expected E = " << expected_E << endl;

    g2.computePriority();
}


void compute_m3_m2(long double& m4, long double& m3, long double& m2,
                   long double& m1, long double& m0, BiGraph& g2) {
    long double caterpillar = 0;
    long double chopsticks = 0;
    long double wedges = 0;

    long double g2_num_edges = g2.num_edges;

    long double g2_num_v1 = g2.num_v1;
    long double g2_num_v2 = g2.num_v2;

#pragma omp parallel for reduction(+ : caterpillar, chopsticks)
    for (int i = 0; i < g2.num_v1; i++) {
        for (auto j : g2.neighbor[i]) {
            // for each edge (i,j)
            caterpillar += (g2.degree[i] - 1) * (g2.degree[j] - 1);
            chopsticks += g2_num_edges - g2.degree[i] - g2.degree[j] + 1;
        }
    }

    chopsticks /= 2;  // each is counted twice

#pragma omp parallel for reduction(+ : wedges)
    for (int i = 0; i < g2.num_nodes(); i++) {
        long double deg_i = g2.degree[i];
        if (deg_i == 0) continue;
        if (g2.is_upper(i)) {
            wedges += (g2_num_v1 - 1) * deg_i * (deg_i - 1) /
                      2;  // is this correct? I am curious.
        } else {
            wedges += (g2_num_v2 - 1) * deg_i * (deg_i - 1) / 2;
        }
    }

    m3 = caterpillar - 4 * m4;

    long double m21 = wedges - 4 * m4 - 2 * m3;

    long double m22 = chopsticks - 2 * m4 - m3;

    m2 = m21 + m22;

    m1 = g2_num_edges * (g2_num_v1 - 1) * (g2_num_v2 - 1) - 4 * m4 - 3 * m3 -
         2 * m2;

    m0 = (g2_num_v1 * (g2_num_v1 - 1) / 2) * (g2_num_v2 * (g2_num_v2 - 1) / 2) -
         m4 - m3 - m2 - m1;

    cout << "m4 = " << m4 << endl;
    cout << "m3 = " << m3 << endl;
    cout << "m2 = " << m2 << endl;
    cout << "m1 = " << m1 << endl;
    cout << "m0 = " << m0 << endl;
}

long double one_round_btf(BiGraph& g, unsigned long seed) {
    long double t0 = omp_get_wtime();

    std::mt19937 rng(std::random_device{}());

    BiGraph g2(g);

    // user side:
    construct_noisy_graph(g, g2, seed);

    long double t1 = omp_get_wtime();

    RR_time += t1 - t0;

    // server side:
    long double m4;

    if (sampling_noisy_graph) {
        cout << "sampling ratio = " << p____ << endl;
        long double sum__ = 0;
        int num_itr = 1;
        for (int xxx = 0; xxx < num_itr; xxx++) {
            // get a sampled subgraph:
            BiGraph g__(g);  // sampled subgraph
            init_genrand(rng());
            for (int i = 0; i < g2.num_v1; i++) {
                for (auto j : g2.neighbor[i]) {
                    if (genrand_real2() < p____) {
                        g__.addEdge(i, j);
                    }
                }
            }
            g__.computePriority();
            cout << "num edges in sampled noisy graph = " << g__.num_edges
                 << endl;
            m4 = BFC_EVP(g__);
            m4 /= (p____ * p____ * p____ * p____);
            sum__ += m4;
        }
        m4 = sum__ / num_itr;

        cout << "estimated m4' = " << m4 << endl;

    } else {
        cout<<"computing m4"<<endl;
        m4 = BFC_EVP(g2);
    }

    cout << "m4 is ready" << endl;

    long double t1x = omp_get_wtime();

    // should we estimate cate and other things independently?

    long double m3 = 0, m2 = 0, m1 = 0, m0 = 0, estimate, mu = exp(Eps);

    compute_m3_m2(m4, m3, m2, m1, m0, g2);

    cout << "computing m3, m2, m1, m0" << endl;

    // this actually means getting bipartite clustering coefficient
    /*
    if(count_cc){
        long double cate_estimate = -4 * mu * mu * mu * m4 + mu * mu * (mu * mu + 3) * m3 -
                   2 * mu * (mu * mu + 1) * m2 + (3 * mu * mu + 1) * m1 -
                   4 * mu * m0;
        cate_estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        // 
        long double btf_estimate = mu * mu * mu * mu * m4 - mu * mu * mu * m3 + mu * mu * m2 -
                   mu * m1 + m0;
        btf_estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        estimate = 4 * btf_estimate / cate_estimate; 

        naive_estis[iteration] = m4 * 4/(m4 * 4 + m3);


    }else 
    */
    if (count_cate) {
        cout << "\tcaterpillar estimation: " << endl;
        estimate = -4 * mu * mu * mu * m4 
                   + mu * mu * (mu * mu + 3) * m3 -
                   2 * mu * (mu * mu + 1) * m2 
                   + (3 * mu * mu + 1) * m1 
                   -4 * mu * m0;
        estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        naive_estis[iteration] = m4 * 4 + m3;
    } 
    else {
        cout << "\tbtf estimation: " << endl;
        estimate = mu * mu * mu * mu * m4 - mu * mu * mu * m3 + mu * mu * m2 -
                   mu * m1 + m0;
        estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        naive_estis[iteration] = m4;
        // the time needed to get m4 on G2 alone.
        naive_server_side += t1x - t1;
    }
    long double t2 = omp_get_wtime();
    // record time elapsed.
    // RR_time += t1-t0; // compute the total time from randomized response.
    server_side_time += t2 - t1;

    return estimate;
}

// the challenge lies in how to compute deg(u, w) = {v \in N(u), v < w } for
// each u < w combination.
long double BFC_EVP(BiGraph& g) {
    long double BTF = 0;
#pragma omp parallel for reduction(+ : BTF)
    for (int u = 0; u < g.num_nodes(); u++) {
        if (g.degree[u] <= 1) continue;

        // cout<<"u  = "<<u<<endl;
        unordered_map<vid_t, int> count_wedge(0);
        for (auto v : g.neighbor[u]) {
            // cout<<"\t v = "<<v<<endl;
            // u -> v
            for (auto w : g.neighbor[v]) {
                // u->v->w
                // cout<<"u, v, w= "<<u<<" "<<"v"<<" "<<w<<endl;
                // this step is not invoked for g3.
                // what if we just use id.
                if (g.com_p(w, v) & g.com_p(w, u)) {  // this is a lot faster.
                    count_wedge[w] = count_wedge[w] + 1;
                }
            }
        }
        // long double btf_u = 0;
        for (auto ele : count_wedge) {
            // cout<<"wedge count = " <<ele.second <<endl;
            if (ele.second >= 2) {
                BTF += (ele.second - 1) * ele.second / 2;
            }
        }
    }

    // cout<<"BTF = "<<BTF<<endl;
    return BTF;
}

long double get_wedges(BiGraph& g) {
    long double wedges = 0;
    for (int i = 0; i < g.num_nodes(); i++) {
        long double deg_i = g.degree[i];
        if (deg_i == 0) continue;
        if (g.is_upper(i)) {
            wedges += deg_i * (deg_i - 1) / 2;
        } else {
            wedges += deg_i * (deg_i - 1) / 2;
        }
    }
    return wedges;
}

long double get_laplace(long double parameter) {
    return stats::rlaplace(0.0, parameter, engine);
}

long double get_cate(BiGraph& g) {
    long double caterpillar = 0;
    for (int i = 0; i < g.num_v1; i++) {
        for (auto j : g.neighbor[i]) {
            caterpillar += (g.degree[i] - 1) * (g.degree[j] - 1);
        }
    }
    return caterpillar;
}

// new approach: 
// we might also need to adopt some budget optimization strategy here. 
// this is a two phase algorithm. 
// first construct the whole noisy graph, and then use it
// is it possible to combine this with the two-round algorithm? 


// this function handles the case when p = 2 
long double wedge_based_two_round_2_K_biclique(BiGraph& g, unsigned long seed) {
    // Phase 0. deg_esti_time records the maximum degree perturbation time.
    // double t0 = omp_get_wtime();
    // cout << "private_estimate_of_degrees(g); " << endl;
    Eps0 = Eps * 0.05;
    // private_estimate_of_degrees(g);
	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}
    // if (eva_comm) {
    //     cout<<"communication cots of uploading degree estimates"<<endl;
    //     communication_cost += g.num_nodes() * sizeof(int);
    // }

    // upload noisy degrees
    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
    cout<<"constructing g2\n";
	construct_noisy_graph(g, g2, seed);  // upload noisy edges
    // unfortunately, this step cannot be run in parallel


    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

    // Phase 2. local counting
    double t2 = omp_get_wtime();
    cout << "local counting" << endl;

    /*
    if (eva_comm) {
        // for each vertex, download all vertex degrees.
        // communication_cost += g.num_nodes() * g.num_nodes() * sizeof(int);
        // for each vertex, download the whole noisy graph
        double byte_per_edge = 8*(log2(g.num_v1) + log2(g.num_v2));
        communication_cost += g2.num_edges * byte_per_edge;
        if(two_noisy_graph_switch){
            communication_cost += g3.num_edges * byte_per_edge;
        }

        // upload common neighbor estimates:
        size_t pairwise_count = g.num_v1 * (g.num_v1 - 1) / 2;
        communication_cost += sizeof(long double) * pairwise_count;
        if (multi_estimator_switch) {
            communication_cost += sizeof(long double) * pairwise_count;
        }

        return 0;
    }
    */

	Eps2 = Eps - Eps1 - Eps0;
    
	gamma__ = (1-p) / (1-2*p);

	long double res___ = 0; 

	// what if we only consider upper vertices ?  --> better efficiency and effect  
	int K = K___;  // we are considering (2, K)-biclique right now

    cout<<"K___ = "<<K___ <<endl;

	int start__, end__;

	start__ = g.num_v1 < g.num_v2 ? 0 : g.num_v1; 
	end__ = g.num_v1 < g.num_v2 ? g.num_v1 : g.num_nodes(); 

    // Probability filtering logic
    if (use_probability_filtering) {
        cout << "Using probability-based filtering..." << endl;
        
        // First pass: compute all f_u_w and esti_var_f values
        struct PairData {
            int u, w;
            long double f_u_w;
            long double esti_var_f;
        };
        
        vector<PairData> all_pairs;
        
        for(int u = start__; u < end__; u++) {
            for(int w = start__; w < end__; w++) {
                if(u <= w) continue;
                
                long double f_u_w, f_w_u;    
                if(two_noisy_graph_switch){
                    f_u_w = locally_compute_f_given_q_and_x_two_graphs(u, w, g, g2, g3);
                    f_w_u = locally_compute_f_given_q_and_x_two_graphs(w, u, g, g2, g3);
                } else {
                    f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
                    if(multi_estimator_switch){
                        f_w_u = locally_compute_f_given_q_and_x(w, u, g, g2);
                    }
                }
                
                long double esti_var_f;
                if(!multi_estimator_switch){
                    esti_var_f = 2 * pow(gamma__,2) / pow(Eps2,2) + p * (1 - p) * deg_estis[u] / pow(1-2*p,2); 
                } else {
                    f_u_w = (f_u_w + f_w_u)/2;
                    long double variance_f_u = 2 * pow(gamma__,2) / pow(Eps2,2);
                    long double variance_f_w = 2 * pow(gamma__,2) / pow(Eps2,2);
                    long double main_fu = p * (1 - p) * deg_estis[u] / pow(1-2*p,2);
                    long double main_fw = p * (1 - p) * deg_estis[w] / pow(1-2*p,2);
                    if(two_noisy_graph_switch){
                        main_fu/=2;
                        main_fw/=2;
                    }
                    variance_f_u += main_fu;
                    variance_f_w += main_fw;
                    esti_var_f = (variance_f_u + variance_f_w) / 4;
                }
                
                all_pairs.push_back({u, w, f_u_w, esti_var_f});
            }
        }
        
        // Filter pairs based on probability P(f_true >= K) >= 90%
        vector<PairData> pairs_to_process;
        long double prob_threshold = 0.90;
        
        for(const auto& pair : all_pairs) {
            // Clamp variance to avoid numerical issues
            long double f_variance = max(pair.esti_var_f, (long double)1e-6);
            
            // Compute P(f_true >= K) using normal approximation
            long double z_score = (K - pair.f_u_w) / sqrt(f_variance);
            long double prob_f_ge_K = 0.5 * erfc(z_score / sqrt(2.0));
            
            if (prob_f_ge_K >= prob_threshold) {
                pairs_to_process.push_back(pair);
            }
        }
        
        cout << "Selected " << pairs_to_process.size() << " pairs out of " << all_pairs.size() 
             << " based on probability filtering (P(f>=" << K << ") >= " << prob_threshold << ")" << endl;
        
        // Second pass: compute local_res for selected pairs
        for(const auto& pair : pairs_to_process) {
            long double local_res = compute_local_res(K, pair.f_u_w, pair.esti_var_f);
            res___ += local_res;
        }
        
    } else {
        // Original logic without filtering
	#pragma omp parallel
	{
	#pragma omp for schedule(static)
		for(int u =start__ ; u <end__ ; u++) {
			for(int w =start__ ; w <end__ ; w++) {
                if(u<=w) continue;

                long double f_u_w, f_w_u;
                if(two_noisy_graph_switch){
                    // when this switch is on, by default we expect multiple estimators.
                    f_u_w = locally_compute_f_given_q_and_x_two_graphs(u, w, g, g2, g3);
                    f_w_u = locally_compute_f_given_q_and_x_two_graphs(w, u, g, g2, g3);

                    // basically getting the same thing using two noisy graphs.
                }else{
                    f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
                    if(multi_estimator_switch){
                        f_w_u = locally_compute_f_given_q_and_x(w, u, g, g2);
                    }
                }

                long double diff1 =0, diff2 = 0;
                
                long double local_res = 0;

                // define some variables
                long double esti_var_f, variance_f_u, variance_f_w, main_fu, main_fw;

                if(!multi_estimator_switch){
                    // single source estimator: 
                    esti_var_f = 2 * pow(gamma__,2)  / pow(Eps2,2) + p * (1 - p) * deg_estis[u] / pow(1-2*p,2); 
                }else{
                    // multi source estimator:  
                    f_u_w = (f_u_w + f_w_u)/2;
                    // esitmate the variance of f_u_w
                    esti_var_f = 0;
                    // these are variance from laplace, not affected by two_noisy_graph_switch
                    variance_f_u = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    variance_f_w = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    // main_fu and main_fw are the variance from local counts
                    main_fu =  p * (1 - p) * deg_estis[u] / pow(1-2*p,2);
                    main_fw =  p * (1 - p) * deg_estis[w] / pow(1-2*p,2);
                    if(two_noisy_graph_switch){
                        // maybe we should increase epsilon 2 to reduce the impact of laplace.
                        main_fu/=2;
                        main_fw/=2;
                    }
                    variance_f_u += main_fu;
                    variance_f_w += main_fw;
                    esti_var_f = (variance_f_u + variance_f_w) / 4;
                }


                // (2, K)-biclique need these moments of the unbiased estimate of f^2:
                // when f' and var(f') are known, we can estimate (P,K)-biclique for any Q.
                local_res = compute_local_res(K, f_u_w, esti_var_f);

				#pragma omp critical
				res___ += local_res; 
			}

		}
		}
	}
    return res___;
}

// Batch version that estimates (2,Q)-biclique counts for multiple Q values efficiently
// This function computes f' and var(f') once, then estimates biclique counts for Q = [4, 5, 6, 7, 8, 9, 10]
std::vector<long double> wedge_based_two_round_2_K_biclique_batch(BiGraph& g, unsigned long seed) {
    // Phase 0. deg_esti_time records the maximum degree perturbation time.
    Eps0 = Eps * 0.05;
	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
    for(int i=0;i<g.num_nodes();i++){
        deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
    }

    // upload noisy degrees
    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
    cout<<"constructing g2\n";
	construct_noisy_graph(g, g2, seed);  // upload noisy edges

    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

    // Phase 2. local counting
    double t2 = omp_get_wtime();
    cout << "local counting (batch mode)" << endl;

    Eps2 = Eps - Eps1 - Eps0;
	gamma__ = (1-p) / (1-2*p);

    // Initialize results for Q = [4, 5, 6, 7, 8, 9, 10]
    std::vector<long double> results(7, 0.0);  // results[0] = Q=4, results[1] = Q=5, etc.

	// what if we only consider upper vertices ?  --> better efficiency and effect  
    int start__, end__; 
    start__ = g.num_v1 < g.num_v2 ? 0 : g.num_v1; 
    end__ = g.num_v1 < g.num_v2 ? g.num_v1 : g.num_nodes(); 

    cout<<"start__ = "<<start__ <<endl;
    cout<<"end__ = "<<end__ <<endl;

    if(start__==0){
        cout<<"n1 n2 tells me to use U(G)"<<endl;
    }else{
        cout<<"n1 n2 tells me to use L(G)"<<endl;
    }

	#pragma omp parallel
	{
        // Local results for each thread
        std::vector<long double> local_results(7, 0.0);
        
	#pragma omp for schedule(static)
		for(int u =start__ ; u <end__ ; u++) {
			for(int w =start__ ; w <end__ ; w++) {
                if(u<=w) continue;

                long double f_u_w, f_w_u;
                if(two_noisy_graph_switch){
                    f_u_w = locally_compute_f_given_q_and_x_two_graphs(u, w, g, g2, g3);
                    f_w_u = locally_compute_f_given_q_and_x_two_graphs(w, u, g, g2, g3);
                }else{
                    f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
                    if(multi_estimator_switch){
                        f_w_u = locally_compute_f_given_q_and_x(w, u, g, g2);
                    }
                }
                
                long double local_res = 0;

                // define some variables
                long double esti_var_f, variance_f_u, variance_f_w, main_fu, main_fw;

                if(!multi_estimator_switch){
                    // single source estimator: 
                    esti_var_f = 2 * pow(gamma__,2)  / pow(Eps2,2) + p * (1 - p) * deg_estis[u] / pow(1-2*p,2); 
                }else{
                    // multi source estimator:  
                    f_u_w = (f_u_w + f_w_u)/2;
                    // esitmate the variance of f_u_w
                    esti_var_f = 0;
                    // these are variance from laplace, not affected by two_noisy_graph_switch
                    variance_f_u = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    variance_f_w = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    // main_fu and main_fw are the variance from local counts
                    main_fu =  p * (1 - p) * deg_estis[u] / pow(1-2*p,2);
                    main_fw =  p * (1 - p) * deg_estis[w] / pow(1-2*p,2);
                    if(two_noisy_graph_switch){
                        // maybe we should increase epsilon 2 to reduce the impact of laplace.
                        main_fu/=2;
                        main_fw/=2;
                    }
                    variance_f_u += main_fu;
                    variance_f_w += main_fw;
                    esti_var_f = (variance_f_u + variance_f_w) / 4;
                }

                // Compute results for all Q values [4, 5, 6, 7, 8, 9, 10]
                for(int q_idx = 0; q_idx < 7; q_idx++) {
                    int Q = q_idx + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
                    long double local_res_q = compute_local_res(Q, f_u_w, esti_var_f);
                    local_results[q_idx] += local_res_q;
                }
            }
        }
        
        // Reduce results from all threads
				#pragma omp critical
        {
            for(int q_idx = 0; q_idx < 7; q_idx++) {
                results[q_idx] += local_results[q_idx];
            }
        }
    }
    
    cout << "Batch results for Q=[4,5,6,7,8,9,10]: ";
    for(int q_idx = 0; q_idx < 7; q_idx++) {
        cout << "Q" << (q_idx + 4) << "=" << results[q_idx] << " ";
    }
    cout << endl;
    
    return results;
}

// Batch version of naive algorithm that estimates (P,Q)-biclique counts for multiple Q values efficiently
// This function builds the noisy graph once, then estimates biclique counts for Q = [4, 5, 6, 7, 8, 9, 10]
std::vector<long double> naive_biclique_batch(BiGraph& g, unsigned long seed, int P___) {
    cout << "Running naive batch algorithm..." << endl;
    
    // Build noisy graph once
    p = 1.0 / (exp(Eps) + 1.0);
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);
    
    // cout << "noisy edges = " << g2.num_edges << endl;
    // cout << "expected E = " << (g.num_edges * (1 - p) + (g.num_v1 * g.num_v2 - g.num_edges) * p) << endl;
    
    // Convert to the format needed by naive algorithm (do this once)
    biGraph convertedGraph = convertBiGraphTobiGraph(g2);
    cout << "Converted graph: n1=" << convertedGraph.n1 << ", n2=" << convertedGraph.n2 << ", m=" << convertedGraph.m << std::endl;
    
    // Initialize results for Q = [4, 5, 6, 7, 8, 9, 10]
    std::vector<long double> results(7, 0.0);
    
    // For each Q value, run exact counting on the same converted graph
    for (int q_idx = 0; q_idx < 7; q_idx++) {
        int Q = q_idx + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
        
        // Run exact counting on the same converted graph
        BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, P___, Q);
        results[q_idx] = counter->exactCount();
        
        cout << "Q" << Q << " naive estimate = " << results[q_idx] << endl;
        
        delete counter;
    }
    
    cout << "Naive batch results for Q=[4,5,6,7,8,9,10]: ";
    for (int q_idx = 0; q_idx < 7; q_idx++) {
        cout << "Q" << (q_idx + 4) << "=" << results[q_idx] << " ";
    }
    cout << endl;
    
    return results;
}

// this function is here to handle when p = 3, sample size = 2*10^6
long double wedge_based_two_round_3_K_biclique(BiGraph& g, unsigned long seed) {
    double t1 = omp_get_wtime();

    Eps0 = Eps * 0.05;

	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  


    // two noisy graph technique
    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

	Eps2 = Eps - Eps1 - Eps0;
    
	long double res___ = 0; 


	int K = K___;  // we are considering (2, K)-biclique right now

    cout<<"p = "<<3 <<endl;
    cout<<"q = "<<K___ <<endl;


    // Calculate the size of the smaller partition
    int smaller_partition_size = std::min(g.num_v1, g.num_v2);

    // Calculate the total number of possible triples in the smaller partition
    long long total_triples = static_cast<long long>(smaller_partition_size) * 
                            (smaller_partition_size - 1) * 
                            (smaller_partition_size - 2) / 6; 


    // Determine how many triples to sample based on the fraction
    // double sample_fraction = 1e-4;

    // Set T to 2 million for better effectiveness
    long long T = 2000000LL;

    bool use_exhaustive = (total_triples <= T);
    double sample_fraction;

    if (use_exhaustive) {
        sample_fraction = 1.0;  // Process all triplets
        cout<<"Using EXHAUSTIVE enumeration (total = " << total_triples << ", T = " << T << ")" << endl;
    } else {
        sample_fraction = (double)T / total_triples;  // Sample 2M triplets
        cout<<"Using REJECTION SAMPLING (T = " << T << " of " << total_triples << " total = " << (sample_fraction * 100) << "%)" << endl;
    }
    
    long long num_triples_to_sample = static_cast<long long>(total_triples * sample_fraction);

    bool is_upper_smaller = (g.num_v1 < g.num_v2 );
    gamma__ = (1-p) / (1-2*p);

    long double res = 0; 

    // Determine the start and end bounds based on partition sizes
    // by default, right now we only consider upper layer.
    int start__ = 0; 
    int end__ = g.num_v1; 

    if (use_exhaustive) {
        // EXHAUSTIVE ENUMERATION: Process all triplets
        cout<<"total triplets: "<<total_triples <<endl;
        cout << "Processing all triplets" <<endl;
        
        #pragma omp parallel
        {
        #pragma omp for schedule(static)
        for (int v1 = start__; v1 < end__ - 2; ++v1) {
            for (int v2 = v1 + 1; v2 < end__ - 1; ++v2) {
                for (int v3 = v2 + 1; v3 < end__; ++v3) {


                // from the neighbors of v1:
                long double f1 = 0, f2= 0, f3 = 0, f12 = 0, f13=0;
                long double local_res = 0, esti_var_f_uvw = 0, fuvw = 0 ;

                if(multi_estimator_switch){
                    // multi-source estimator
                    for(auto nb: g.neighbor[v1]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);

                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                        }
                        f1 += A1 * A2; 
                        f12 += A1; 
                        f13 += A2; 
                    }
                    // two_noisy_graph_switch does not change GS and Lap noise
                    f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                    f12 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f13 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    
                    long double f21 = 0, f23=0;
                    for(auto nb: g.neighbor[v2]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);
                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                        }
                        f2 += A1 * A2; 
                        f21 += A1; 
                        f23 += A2;          
                    }
                    // two_noisy_graph_switch does not change GS and Lap noise
                    f2 += stats::rlaplace(0.0,  (gamma__*gamma__/Eps2), engine); 
                    f21 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f23 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    long double f31 = 0, f32=0;
                    for(auto nb: g.neighbor[v3]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                        }
                        f3 += A1 * A2; 
                        f31 += A1; 
                        f32 += A2;    

                    }
                    f3 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                    f31 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f32 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

                    // averaging
                    fuvw = (f1 + f2 + f3 )/3 ; 

                    long double var_phi = p * (1-p) / pow(1-2*p, 2); 
                    if (two_noisy_graph_switch) {
                        var_phi/=2;
                    }

                    long double esti_var_f1, esti_var_f2, esti_var_f3; 

                    // this should always be true
                    bool improvement = true;
                    if(improvement){
                        // averaging f12 and f21 is useful too!
                        esti_var_f1 = var_phi * (f12 + f21 + f13 + f31)/2 ; 
                        esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                        esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                        esti_var_f2 = var_phi * (f12 + f21 + f23 + f32)/2 ; 
                        esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                        esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                        esti_var_f3 = var_phi * (f13 + f31 + f23 + f32)/2 ; 
                        esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                        esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                    }else{ 
                        // this will be true when switch = 3
                        esti_var_f1 = var_phi * (f12 + f13) ; 
                        esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                        esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                        esti_var_f2 = var_phi * (f12 + f23) ; 
                        esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                        esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                        esti_var_f3 = var_phi * (f13 + f23) ; 
                        esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                        esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                    }

                    
                    // include C2(vi, vj)
                    f12 = (f12 + f21)/2;
                    f13 = (f13 + f31)/2;
                    f23 = (f23 + f32)/2;

                    // compute terms from var(f1), var(f2), and var(f3)
                    esti_var_f_uvw = (esti_var_f1 + esti_var_f2 + esti_var_f3);

                    // consider the co-variance of f1, f2, and f3.
                    esti_var_f_uvw += 2 * var_phi * (f12 + f13 + f23); 
                    esti_var_f_uvw /= 9;
                }else{
                    // single-source estimator: f1
                    // bool use_min_degree_vertex = false; // Switch to enable/disable min degree selection
                    
                    // if (use_min_degree_vertex) {
                    //     // Choose vertex with minimum estimated degree
                    //     // Swap vertices so that the one with minimum degree becomes v1
                    //     if (deg_estis[v2] < deg_estis[v1]) {
                    //         std::swap(v1, v2);
                    //     }
                    //     if (deg_estis[v3] < deg_estis[v1]) {
                    //         std::swap(v1, v3);
                    //     }
                    // }
                    
                    // Now v1 is either the original v1 or the vertex with minimum estimated degree
                    for(auto nb: g.neighbor[v1]){
                        long double A1 = g2.has(nb, v2) ? 1 : 0 ; 
                        A1 = (A1-p) / (1-2*p); 

                        long double A2 = g2.has(nb, v3) ? 1 : 0 ; 
                        A2 = (A2-p) / (1-2*p); 

                        f1 += A1 * A2;
                    }
                    f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 

                    // no averaging
                    fuvw = f1;
                    
                    // estimate the variance of f(u, v, w)

                    long double var_phi = p * (1-p) / pow(1-2*p, 2); 

                    long double esti_var_f1 = var_phi * (f12 + f13) ; 
                    esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                    esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                    esti_var_f_uvw = esti_var_f1;
                }

                local_res = compute_local_res(K, fuvw, esti_var_f_uvw);

                #pragma omp critical
                res += local_res; 
            }
        }
    }
    }
    } else {
        // REJECTION SAMPLING: Sample exactly T triplets
        cout<<"total triplets: "<<total_triples <<endl;
        cout << "Target sample size = " << T <<endl;
        
        // Random number generation for rejection sampling
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> vertex_dis(0, end__ - 1);
        
        // Hash set to store selected triplets (for uniqueness)
        std::unordered_set<long long> selected_triplets;
        
        // Function to encode triplet as unique ID
        auto encode_triplet = [end__](int v1, int v2, int v3) -> long long {
            // Ensure v1 < v2 < v3 for consistent encoding
            if (v1 > v2) std::swap(v1, v2);
            if (v2 > v3) std::swap(v2, v3);
            if (v1 > v2) std::swap(v1, v2);
            
            // Encode as: v1 * n^2 + v2 * n + v3 (where n = end__)
            long long n = end__;
            return v1 * n * n + v2 * n + v3;
        };

        // Rejection sampling to get exactly T unique triplets
        while (selected_triplets.size() < T) {
            // Generate random triplet
            int v1 = vertex_dis(gen);
            int v2 = vertex_dis(gen);
            int v3 = vertex_dis(gen);
            
            // Ensure v1 < v2 < v3 (valid triplet)
            if (v1 >= v2 || v2 >= v3 || v1 >= v3) continue;
            
            long long triplet_id = encode_triplet(v1, v2, v3);
            
            // Check if already selected
            if (selected_triplets.find(triplet_id) != selected_triplets.end()) continue;
            
            // Add to selected set
            selected_triplets.insert(triplet_id);
            
            // Process this triplet
            // from the neighbors of v1:
            long double f1 = 0, f2= 0, f3 = 0, f12 = 0, f13=0;
            long double local_res = 0, esti_var_f_uvw = 0, fuvw = 0 ;

            if(multi_estimator_switch){
                // multi-source estimator
                for(auto nb: g.neighbor[v1]){
                    long double A1 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                    long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);

                    if (two_noisy_graph_switch) {
                        A1 = (A1 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                        A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                    }
                    f1 += A1 * A2; 
                    f12 += A1; 
                    f13 += A2; 
                }
                // two_noisy_graph_switch does not change GS and Lap noise
                f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                f12 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                f13 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                
                long double f21 = 0, f23=0;
                for(auto nb: g.neighbor[v2]){
                    long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                    long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);
                    if (two_noisy_graph_switch) {
                        A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                        A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                    }
                    f2 += A1 * A2; 
                    f21 += A1; 
                    f23 += A2;          
                }
                // two_noisy_graph_switch does not change GS and Lap noise
                f2 += stats::rlaplace(0.0,  (gamma__*gamma__/Eps2), engine); 
                f21 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                f23 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                long double f31 = 0, f32=0;
                for(auto nb: g.neighbor[v3]){
                    long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                    long double A2 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                    if (two_noisy_graph_switch) {
                        A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                        A2 = (A2 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                    }
                    f3 += A1 * A2; 
                    f31 += A1; 
                    f32 += A2;    

                }
                f3 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                f31 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                f32 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

                // averaging
                fuvw = (f1 + f2 + f3 )/3 ; 

                long double var_phi = p * (1-p) / pow(1-2*p, 2); 
                if (two_noisy_graph_switch) {
                    var_phi/=2;
                }

                long double esti_var_f1, esti_var_f2, esti_var_f3; 

                // this should always be true
                bool improvement = true;
                if(improvement){
                    // averaging f12 and f21 is useful too!
                    esti_var_f1 = var_phi * (f12 + f21 + f13 + f31)/2 ; 
                    esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                    esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                    esti_var_f2 = var_phi * (f12 + f21 + f23 + f32)/2 ; 
                    esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                    esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                    esti_var_f3 = var_phi * (f13 + f31 + f23 + f32)/2 ; 
                    esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                    esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                }else{ 
                    // this will be true when switch = 3
                    esti_var_f1 = var_phi * (f12 + f13) ; 
                    esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                    esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                    esti_var_f2 = var_phi * (f12 + f23) ; 
                    esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                    esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                    esti_var_f3 = var_phi * (f13 + f23) ; 
                    esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                    esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                }

                
                // include C2(vi, vj)
                f12 = (f12 + f21)/2;
                f13 = (f13 + f31)/2;
                f23 = (f23 + f32)/2;

                // compute terms from var(f1), var(f2), and var(f3)
                esti_var_f_uvw = (esti_var_f1 + esti_var_f2 + esti_var_f3);

                // consider the co-variance of f1, f2, and f3.
                esti_var_f_uvw += 2 * var_phi * (f12 + f13 + f23); 
                esti_var_f_uvw /= 9;
            }else{
                // single-source estimator: f1
                // bool use_min_degree_vertex = true; // Switch to enable/disable min degree selection
                
                // if (use_min_degree_vertex) {
                //     // Choose vertex with minimum estimated degree
                //     // Swap vertices so that the one with minimum degree becomes v1
                //     if (deg_estis[v2] < deg_estis[v1]) {
                //         std::swap(v1, v2);
                //     }
                //     if (deg_estis[v3] < deg_estis[v1]) {
                //         std::swap(v1, v3);
                //     }
                // }
                
                // Now v1 is either the original v1 or the vertex with minimum estimated degree
                for(auto nb: g.neighbor[v1]){
                    long double A1 = g2.has(nb, v2) ? 1 : 0 ; 
                    A1 = (A1-p) / (1-2*p); 

                    long double A2 = g2.has(nb, v3) ? 1 : 0 ; 
                    A2 = (A2-p) / (1-2*p); 

                    f1 += A1 * A2;
                }
                f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 

                // no averaging
                fuvw = f1;
                
                // estimate the variance of f(u, v, w)

                long double var_phi = p * (1-p) / pow(1-2*p, 2); 

                long double esti_var_f1 = var_phi * (f12 + f13) ; 
                esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                esti_var_f_uvw = esti_var_f1;
            }

            local_res = compute_local_res(K, fuvw, esti_var_f_uvw);
            res += local_res;
        }
    }

    // assert(count == num_triples_to_sample);
    double t2 = omp_get_wtime();
    cout<<"time  = "<<t2 - t1 <<endl;
    return res / sample_fraction;

}

// No-sampling batch version that processes all triplets for Q = 4 to 10
std::vector<long double> wedge_based_two_round_3_K_biclique_no_sampling_batch(BiGraph& g, unsigned long seed) {
    double t1 = omp_get_wtime();

    Eps0 = Eps * 0.05;

    vector<long double> deg_estis; 
    deg_estis.resize(g.num_nodes());
    for(int i=0;i<g.num_nodes();i++){
        deg_estis[i] = g.degree[i] + stats::rlaplace(0.0, 1/(Eps0), engine);
    }

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Set the flip probability for randomized response
    p = 1.0 / (exp(Eps1) + 1.0);

    // two noisy graph technique
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);  // upload noisy edges

    // two noisy graph technique
    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

    Eps2 = Eps - Eps1 - Eps0;
    
    // Initialize results for Q = [4, 5, 6, 7, 8, 9, 10]
    std::vector<long double> results(7, 0.0);

    cout<<"p = "<<3 <<endl;
    cout<<"q = [4,5,6,7,8,9,10]" <<endl;

    // Calculate the size of the smaller partition
    int smaller_partition_size = std::min(g.num_v1, g.num_v2);

    // Calculate the total number of possible triples in the smaller partition
    long long total_triples = static_cast<long long>(smaller_partition_size) * 
                            (smaller_partition_size - 1) * 
                            (smaller_partition_size - 2) / 6; 

    cout<<"total triplets: "<<total_triples <<endl;
    cout << "Processing ALL triplets (no sampling)" <<endl;

    bool is_upper_smaller = (g.num_v1 < g.num_v2 );
    gamma__ = (1-p) / (1-2*p);

    // Process all triplets using exhaustive enumeration
    #pragma omp parallel
    {
        // Local results for each thread
        std::vector<long double> local_results(7, 0.0);
        
        #pragma omp for schedule(static)
        for (int v1 = 0; v1 < smaller_partition_size - 2; ++v1) {
            for (int v2 = v1 + 1; v2 < smaller_partition_size - 1; ++v2) {
                for (int v3 = v2 + 1; v3 < smaller_partition_size; ++v3) {
                    
                    // Process this triplet with proper multi_estimator_switch logic
                    long double f1 = 0, f2 = 0, f3 = 0, f12 = 0, f13 = 0, f23 = 0;
                    long double f21 = 0, f31 = 0, f32 = 0;
                    long double fuvw = 0;
                    long double esti_var_f_uvw = 0;
                    
                    if(multi_estimator_switch){
                    
                        // Multi-source estimator
                        // Compute f1 (same as batch function)
                    for (auto nb : g.neighbor[v1]) {
                        long double A1 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);

                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                        }
                        f1 += A1 * A2; 
                        f12 += A1; 
                        f13 += A2; 
                    }
                    f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                    f12 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f13 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    
                    // Compute f2 (same as batch function)
                    for(auto nb: g.neighbor[v2]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);
                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                        }
                        f2 += A1 * A2; 
                        f21 += A1; 
                        f23 += A2;          
                    }
                    f2 += stats::rlaplace(0.0,  (gamma__*gamma__/Eps2), engine); 
                    f21 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f23 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    
                    // Compute f3 (same as batch function)
                    for(auto nb: g.neighbor[v3]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                        }
                        f3 += A1 * A2; 
                        f31 += A1; 
                        f32 += A2;    
                    }
                    f3 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                    f31 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f32 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

                        // averaging
                        fuvw = (f1 + f2 + f3 )/3 ; 

                        long double var_phi = p * (1-p) / pow(1-2*p, 2); 
                        if (two_noisy_graph_switch) {
                            var_phi/=2;
                        }

                        long double esti_var_f1, esti_var_f2, esti_var_f3; 

                        bool improvement = true;
                        if(improvement){
                            // averaging f12 and f21 is useful too!
                            esti_var_f1 = var_phi * (f12 + f21 + f13 + f31)/2 ; 
                            esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                            esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                            esti_var_f2 = var_phi * (f12 + f21 + f23 + f32)/2 ; 
                            esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                            esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                            esti_var_f3 = var_phi * (f13 + f31 + f23 + f32)/2 ; 
                            esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                            esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                        }

                        esti_var_f_uvw = (esti_var_f1 + esti_var_f2 + esti_var_f3) / 3;
                    } else {
                        // Single-source estimator: f1 only
                        for (auto nb : g.neighbor[v1]) {
                            long double A1 = g2.has(nb, v2) ? 1 : 0 ; 
                            A1 = (A1-p) / (1-2*p); 

                            long double A2 = g2.has(nb, v3) ? 1 : 0 ; 
                            A2 = (A2-p) / (1-2*p); 

                            f1 += A1 * A2;
                        }
                        f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 

                        // no averaging
                        fuvw = f1;
                        
                        // estimate the variance of f(u, v, w)
                        long double var_phi = p * (1-p) / pow(1-2*p, 2); 

                        long double esti_var_f1 = var_phi * (f12 + f13) ; 
                        esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                        esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                        esti_var_f_uvw = esti_var_f1;
                    }
                    
                    // Compute results for all Q values [4, 5, 6, 7, 8, 9, 10]
                    for (int q_idx = 0; q_idx < 7; q_idx++) {
                        int K = q_idx + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
                        long double local_res_q = compute_local_res(K, fuvw, esti_var_f_uvw);
                        local_results[q_idx] += local_res_q;
                    }
                }
            }
        }
        
        // Reduce results from all threads
        #pragma omp critical
        {
            for (int q_idx = 0; q_idx < 7; q_idx++) {
                results[q_idx] += local_results[q_idx];
            }
        }
    }
        
    double t2 = omp_get_wtime();
    cout<<"time  = "<<t2 - t1 <<endl;
    
    cout << "No-sampling batch results for Q=[4,5,6,7,8,9,10]: ";
    for(int q_idx = 0; q_idx < 7; q_idx++) {
        cout << "Q" << (q_idx + 4) << "=" << results[q_idx] << " ";
    }
    cout << endl;
    
    return results;
}

// Improved version with rejection sampling for P=3
long double wedge_based_two_round_3_K_biclique_rejection_sampling(BiGraph& g, unsigned long seed) {
    double t1 = omp_get_wtime();

    Eps0 = Eps * 0.05;

    vector<long double> deg_estis; 
    deg_estis.resize(g.num_nodes());
    for(int i=0;i<g.num_nodes();i++){
        deg_estis[i] = g.degree[i];
    }

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Set the flip probability for randomized response
    p = 1.0 / (exp(Eps1) + 1.0);

    // two noisy graph technique
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);  // upload noisy edges

    // two noisy graph technique
    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

    Eps2 = Eps - Eps1 - Eps0;
    
    long double res___ = 0; 

    int K = K___;  // we are considering (3, K)-biclique right now

    cout<<"p = "<<3 <<endl;
    cout<<"q = "<<K___ <<endl;

    // Calculate the size of the smaller partition
    int smaller_partition_size = std::min(g.num_v1, g.num_v2);

    // Calculate the total number of possible triples in the smaller partition
    long long total_triples = static_cast<long long>(smaller_partition_size) * 
                            (smaller_partition_size - 1) * 
                            (smaller_partition_size - 2) / 6; 

    // Target number of triplets to sample (2 million, or all if less)
    long long T = std::min(2000000LL, total_triples);  // Use 2M or all triplets if less
    
    // Random number generation for rejection sampling
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> vertex_dis(0, smaller_partition_size - 1);
    
    // Hash set to store selected triplets (for uniqueness)
    std::unordered_set<long long> selected_triplets;
    
    // Function to encode triplet as unique ID
    auto encode_triplet = [](int v1, int v2, int v3) -> long long {
        // Ensure v1 < v2 < v3 for consistent encoding
        if (v1 > v2) std::swap(v1, v2);
        if (v2 > v3) std::swap(v2, v3);
        if (v1 > v2) std::swap(v1, v2);
        
        // Encode as: v1 * n^2 + v2 * n + v3
        return static_cast<long long>(v1) * 1000000 + v2 * 1000 + v3;
    };

    // Rejection sampling to get exactly T unique triplets
    cout<<"total triplets: "<<total_triples <<endl;
    if (T == total_triples) {
        cout << "Target sample size = " << T << " (all triplets)" <<endl;
    } else {
        cout << "Target sample size = " << T << " (2M max)" <<endl;
    }
    
    while (selected_triplets.size() < T) {
        // Generate random triplet
        int v1 = vertex_dis(gen);
        int v2 = vertex_dis(gen);
        int v3 = vertex_dis(gen);
        
        // Ensure v1 < v2 < v3 (valid triplet)
        if (v1 >= v2 || v2 >= v3 || v1 >= v3) continue;
        
        long long triplet_id = encode_triplet(v1, v2, v3);
        
        // Only add if not already selected (rejection sampling)
        if (selected_triplets.find(triplet_id) == selected_triplets.end()) {
            selected_triplets.insert(triplet_id);
        }
    }
    
    cout << "Successfully sampled " << selected_triplets.size() << " unique triplets" << endl;

    bool is_upper_smaller = (g.num_v1 < g.num_v2 );
    gamma__ = (1-p) / (1-2*p);

    long double res = 0; 

    // Process the selected triplets
    #pragma omp parallel
    {
        long double local_res = 0;
        
        // Convert set to vector for parallel processing
        std::vector<long long> triplet_ids(selected_triplets.begin(), selected_triplets.end());
        
        #pragma omp for schedule(static)
        for (size_t i = 0; i < triplet_ids.size(); ++i) {
            long long triplet_id = triplet_ids[i];
            
            // Decode triplet ID back to vertices
            int v3 = triplet_id % 1000;
            triplet_id /= 1000;
            int v2 = triplet_id % 1000;
            triplet_id /= 1000;
            int v1 = triplet_id;

            // Process this triplet with proper multi_estimator_switch logic
            long double f1 = 0, f2 = 0, f3 = 0, f12 = 0, f13 = 0, f23 = 0;
            long double f21 = 0, f31 = 0, f32 = 0;
            long double fuvw = 0;
            long double esti_var_f_uvw = 0;
            
            if(multi_estimator_switch){
                // Multi-source estimator
                // Compute f1 (same as batch function)
            for (auto nb : g.neighbor[v1]) {
                long double A1 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);

                if (two_noisy_graph_switch) {
                    A1 = (A1 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                    A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                }
                f1 += A1 * A2; 
                f12 += A1; 
                f13 += A2; 
            }
            f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
            f12 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            f13 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            
            // Compute f2 (same as batch function)
            for(auto nb: g.neighbor[v2]){
                long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);
                if (two_noisy_graph_switch) {
                    A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                    A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                }
                f2 += A1 * A2; 
                f21 += A1; 
                f23 += A2;          
            }
            f2 += stats::rlaplace(0.0,  (gamma__*gamma__/Eps2), engine); 
            f21 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            f23 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            
            // Compute f3 (same as batch function)
            for(auto nb: g.neighbor[v3]){
                long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                long double A2 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                if (two_noisy_graph_switch) {
                    A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                    A2 = (A2 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                }
                f3 += A1 * A2; 
                f31 += A1; 
                f32 += A2;    
            }
            f3 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
            f31 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            f32 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

                // averaging
                fuvw = (f1 + f2 + f3 )/3 ; 

                long double var_phi = p * (1-p) / pow(1-2*p, 2); 
                if (two_noisy_graph_switch) {
                    var_phi/=2;
                }

                long double esti_var_f1, esti_var_f2, esti_var_f3; 

                bool improvement = true;
                if(improvement){
                    // averaging f12 and f21 is useful too!
                    esti_var_f1 = var_phi * (f12 + f21 + f13 + f31)/2 ; 
                    esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                    esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                    esti_var_f2 = var_phi * (f12 + f21 + f23 + f32)/2 ; 
                    esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                    esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                    esti_var_f3 = var_phi * (f13 + f31 + f23 + f32)/2 ; 
                    esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                    esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                }

                esti_var_f_uvw = (esti_var_f1 + esti_var_f2 + esti_var_f3) / 3;
            } else {
                // Single-source estimator: f1 only
                for (auto nb : g.neighbor[v1]) {
                    long double A1 = g2.has(nb, v2) ? 1 : 0 ; 
                    A1 = (A1-p) / (1-2*p); 

                    long double A2 = g2.has(nb, v3) ? 1 : 0 ; 
                    A2 = (A2-p) / (1-2*p); 

                    f1 += A1 * A2;
                }
                f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 

                // no averaging
                fuvw = f1;
                
                // estimate the variance of f(u, v, w)
                long double var_phi = p * (1-p) / pow(1-2*p, 2); 

                long double esti_var_f1 = var_phi * (f12 + f13) ; 
                esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                esti_var_f_uvw = esti_var_f1;
            }

            local_res += compute_local_res(K, fuvw, esti_var_f_uvw);
        }
        
        #pragma omp critical
        res += local_res;
    }

    double t2 = omp_get_wtime();
    cout<<"time  = "<<t2 - t1 <<endl;
    
    // Scale by the sampling ratio (T / total_triples)
    long double sample_fraction = static_cast<long double>(T) / total_triples;
    return res / sample_fraction;
}


// Batch version of wedge_based_two_round_3_K_biclique_rejection_sampling
std::vector<long double> wedge_based_two_round_3_K_biclique_rejection_sampling_batch(BiGraph& g, unsigned long seed) {
    double t1 = omp_get_wtime();

    Eps0 = Eps * 0.05;

    vector<long double> deg_estis; 
    deg_estis.resize(g.num_nodes());
    for(int i=0;i<g.num_nodes();i++){
        deg_estis[i] = g.degree[i];
    }

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Set the flip probability for randomized response
    p = 1.0 / (exp(Eps1) + 1.0);

    // two noisy graph technique
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);  // upload noisy edges

    // two noisy graph technique
    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

    Eps2 = Eps - Eps1 - Eps0;
    
    // Initialize results for Q = [4, 5, 6, 7, 8, 9, 10]
    std::vector<long double> results(7, 0.0);

    cout<<"p = "<<3 <<endl;
    cout<<"q = [4,5,6,7,8,9,10]" <<endl;

    // Calculate the size of the smaller partition
    int smaller_partition_size = std::min(g.num_v1, g.num_v2);

    // Calculate the total number of possible triples in the smaller partition
    long long total_triples = static_cast<long long>(smaller_partition_size) * 
                            (smaller_partition_size - 1) * 
                            (smaller_partition_size - 2) / 6; 

    // Target number of triplets to sample (2 million, or all if less)
    long long T = std::min(2000000LL, total_triples);  // Use 2M or all triplets if less
    
    // Random number generation for rejection sampling
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> vertex_dis(0, smaller_partition_size - 1);
    
    // Hash set to store selected triplets (for uniqueness)
    std::unordered_set<long long> selected_triplets;
    
    // Function to encode triplet as unique ID
    auto encode_triplet = [](int v1, int v2, int v3) -> long long {
        // Ensure v1 < v2 < v3 for consistent encoding
        if (v1 > v2) std::swap(v1, v2);
        if (v2 > v3) std::swap(v2, v3);
        if (v1 > v2) std::swap(v1, v2);
        
        // Encode as: v1 * n^2 + v2 * n + v3
        return static_cast<long long>(v1) * 1000000 + v2 * 1000 + v3;
    };

    // Rejection sampling to get exactly T unique triplets
    cout<<"total triplets: "<<total_triples <<endl;
    if (T == total_triples) {
        cout << "Target sample size = " << T << " (all triplets)" <<endl;
    } else {
        cout << "Target sample size = " << T << " (2M max)" <<endl;
    }
    
    while (selected_triplets.size() < T) {
        // Generate random triplet
        int v1 = vertex_dis(gen);
        int v2 = vertex_dis(gen);
        int v3 = vertex_dis(gen);
        
        // Ensure v1 < v2 < v3 (valid triplet)
        if (v1 >= v2 || v2 >= v3 || v1 >= v3) continue;
        
        long long triplet_id = encode_triplet(v1, v2, v3);
        
        // Only add if not already selected (rejection sampling)
        if (selected_triplets.find(triplet_id) == selected_triplets.end()) {
            selected_triplets.insert(triplet_id);
        }
    }
    
    cout << "Successfully sampled " << selected_triplets.size() << " unique triplets" << endl;

    bool is_upper_smaller = (g.num_v1 < g.num_v2 );
    gamma__ = (1-p) / (1-2*p);

    // Process the selected triplets
    #pragma omp parallel
    {
        // Local results for each thread
        std::vector<long double> local_results(7, 0.0);
        
        // Convert set to vector for parallel processing
        std::vector<long long> triplet_ids(selected_triplets.begin(), selected_triplets.end());
        
        #pragma omp for schedule(static)
        for (size_t i = 0; i < triplet_ids.size(); ++i) {
            long long triplet_id = triplet_ids[i];
            
            // Decode triplet ID back to vertices
            int v3 = triplet_id % 1000;
            int v2 = (triplet_id / 1000) % 1000;
            int v1 = triplet_id / 1000000;
            
            // Process this triplet with proper multi_estimator_switch logic
            long double f1 = 0, f2 = 0, f3 = 0, f12 = 0, f13 = 0, f23 = 0;
            long double f21 = 0, f31 = 0, f32 = 0;
            long double fuvw = 0;
            long double esti_var_f_uvw = 0;
            
            if(multi_estimator_switch){
                // Multi-source estimator
                // Compute f1 (same as original)
            for (auto nb : g.neighbor[v1]) {
                long double A1 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);

                if (two_noisy_graph_switch) {
                    A1 = (A1 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                    A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                }
                f1 += A1 * A2; 
                f12 += A1; 
                f13 += A2; 
            }
            f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
            f12 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            f13 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            
            // Compute f2 (same as original)
            for(auto nb: g.neighbor[v2]){
                long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);
                if (two_noisy_graph_switch) {
                    A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                    A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                }
                f2 += A1 * A2; 
                f21 += A1; 
                f23 += A2;          
            }
            f2 += stats::rlaplace(0.0,  (gamma__*gamma__/Eps2), engine); 
            f21 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            f23 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            
            // Compute f3 (same as original)
            for(auto nb: g.neighbor[v3]){
                long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                long double A2 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                if (two_noisy_graph_switch) {
                    A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                    A2 = (A2 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                }
                f3 += A1 * A2; 
                f31 += A1; 
                f32 += A2;    
            }
            f3 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
            f31 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
            f32 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

                // averaging
                fuvw = (f1 + f2 + f3 )/3 ; 

                long double var_phi = p * (1-p) / pow(1-2*p, 2); 
                if (two_noisy_graph_switch) {
                    var_phi/=2;
                }

                long double esti_var_f1, esti_var_f2, esti_var_f3; 

                bool improvement = true;
                if(improvement){
                    // averaging f12 and f21 is useful too!
                    esti_var_f1 = var_phi * (f12 + f21 + f13 + f31)/2 ; 
                    esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                    esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                    esti_var_f2 = var_phi * (f12 + f21 + f23 + f32)/2 ; 
                    esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                    esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                    esti_var_f3 = var_phi * (f13 + f31 + f23 + f32)/2 ; 
                    esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                    esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                }

                esti_var_f_uvw = (esti_var_f1 + esti_var_f2 + esti_var_f3) / 3;
            } else {
                // Single-source estimator: f1 only
                for (auto nb : g.neighbor[v1]) {
                    long double A1 = g2.has(nb, v2) ? 1 : 0 ; 
                    A1 = (A1-p) / (1-2*p); 

                    long double A2 = g2.has(nb, v3) ? 1 : 0 ; 
                    A2 = (A2-p) / (1-2*p); 

                    f1 += A1 * A2;
                }
                f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 

                // no averaging
                fuvw = f1;
                
                // estimate the variance of f(u, v, w)
                long double var_phi = p * (1-p) / pow(1-2*p, 2); 

                long double esti_var_f1 = var_phi * (f12 + f13) ; 
                esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                esti_var_f_uvw = esti_var_f1;
            }
            
            // Compute results for all Q values [4, 5, 6, 7, 8, 9, 10]
            for (int q_idx = 0; q_idx < 7; q_idx++) {
                int K = q_idx + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
                long double local_res_q = compute_local_res(K, fuvw, esti_var_f_uvw);
                local_results[q_idx] += local_res_q;
            }
        }
        
        // Reduce results from all threads
        #pragma omp critical
        {
            for (int q_idx = 0; q_idx < 7; q_idx++) {
                results[q_idx] += local_results[q_idx];
            }
        }
    }
        
    // Scale by the sampling ratio (same as original function)
    long double sample_fraction = static_cast<long double>(T) / total_triples;
    for (int q_idx = 0; q_idx < 7; q_idx++) {
        results[q_idx] = results[q_idx] / sample_fraction;
    }
    
    double t2 = omp_get_wtime();
    cout<<"time  = "<<t2 - t1 <<endl;
    
    cout << "Batch results for Q=[4,5,6,7,8,9,10]: ";
    for(int q_idx = 0; q_idx < 7; q_idx++) {
        cout << "Q" << (q_idx + 4) << "=" << results[q_idx] << " ";
    }
    cout << endl;
    
    return results;
}


// this function handles p > 3. In experiment, we focus on p = 4.
long double wedge_based_two_round_general_biclique(BiGraph& g, 
    unsigned long seed, int P___, int K___, long long T) {
    double t1 = omp_get_wtime();

    Eps0 = Eps * 0.05;

	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    gamma__ = (1-p) / (1-2*p);
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  

	Eps2 = Eps - Eps1 - Eps0;
    
	long double res___ = 0; 

    int N = g.num_v1 ; 
    
    // For general P, always use sampling since p-tuples grow exponentially
    // T is now passed as a parameter
    // long long T = 10000LL;
    cout << "Using SAMPLING for P=" << P___ << " (T = " << T << " p-tuples)" << endl;
    
    // Calculate total possible p-tuples for reference
    long long total_p_tuples = 1;
    for (int i = 0; i < P___; i++) {
        total_p_tuples = total_p_tuples * (N - i) / (i + 1);
    }
    cout << "Total p-tuples: " << total_p_tuples << ", sampling " << T << endl;
    
    // Sample T p-tuples using rejection sampling
    vector<vector<int>> subsets;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> vertex_dis(0, N - 1);
    
    // Hash set to store selected p-tuples (for uniqueness)
    std::unordered_set<long long> selected_p_tuples;
    
    // Function to encode p-tuple as unique ID
    auto encode_p_tuple = [N](const vector<int>& p_tuple) -> long long {
        vector<int> sorted_tuple = p_tuple;
        std::sort(sorted_tuple.begin(), sorted_tuple.end());
        
        // Encode as: v1 * n^(p-1) + v2 * n^(p-2) + ... + vp
        long long encoding = 0;
        for (int i = 0; i < sorted_tuple.size(); i++) {
            encoding = encoding * N + sorted_tuple[i];
        }
        return encoding;
    };
    
    // Rejection sampling to get exactly T unique p-tuples
    while (selected_p_tuples.size() < T) {
        // Generate random p-tuple
        vector<int> p_tuple;
        std::unordered_set<int> used_vertices;
        
        // Generate P___ distinct vertices
        while (p_tuple.size() < P___) {
            int vertex = vertex_dis(gen);
            if (used_vertices.find(vertex) == used_vertices.end()) {
                p_tuple.push_back(vertex);
                used_vertices.insert(vertex);
            }
        }
        
        long long p_tuple_id = encode_p_tuple(p_tuple);
        
        // Check if already selected
        if (selected_p_tuples.find(p_tuple_id) != selected_p_tuples.end()) continue;
        
        // Add to selected set
        selected_p_tuples.insert(p_tuple_id);
        subsets.push_back(p_tuple);
    }
    
    double sample_fraction = (double)T / total_p_tuples;
    cout << "Sample fraction: " << (sample_fraction * 100) << "%" << endl;



    #pragma omp parallel
	{
    #pragma omp for schedule(static)
    for (int i__ = 0;  i__ < subsets.size(); ++i__) {
        const auto& subset = subsets[i__];


        // picking the min degree vertex as source
        int v1 = subset[0];
        long double min_deg = deg_estis[v1];
        for (int i = 1; i < P___; ++i) {
            if (deg_estis[subset[i]] < min_deg) {
                min_deg = deg_estis[subset[i]];
                v1 = subset[i];
            }
        }
        
        // std::cout << "Min deg vertex (v1): " << v1 << ", deg_est: " << min_deg << "\n";

        long double local_res;


        // if we estimate f_avg, then var(f_avg) will be hard to estimate due to a lot of covariances. 
        // what we can do instead, is to ask each vertex to estimate f' and var(f'). 
        // use compute_local_res to estimate local biclique count for each vertex. 
        // at last aggregate. 
        
        if(multi_estimator_switch) {
            // Multi-source estimator: use all vertices in subset as sources
            vector<long double> f_estimates(P___);
            vector<long double> f1_estimates(P___);
            
            // Compute true common neighbor count for the current subset
            // Use the first vertex in subset as reference
            int reference_vertex = subset[0];
            int true_common_neighbors_subset = 0;
            for(auto nb: g.neighbor[reference_vertex]) {
                bool is_common_neighbor = true;
                for(int i = 1; i < P___; i++) {  // Start from 1 since we already know nb is neighbor of subset[0]
                    if(!g.has(nb, subset[i])){
                        is_common_neighbor = false;
                        break;
                    }
                }
                if(is_common_neighbor){
                    true_common_neighbors_subset++;
                }
            }
            
            // Print subset information to file (single-threaded)
            static std::ofstream debug_file("/tmp/biclique_debug.txt");
            debug_file << "--- 4-tuple: " << subset[0] << "," << subset[1] << "," << subset[2] << "," << subset[3] 
                       << " true=" << true_common_neighbors_subset << " ---" << std::endl;
            debug_file.flush();
            
            for(int source_idx = 0; source_idx < P___; source_idx++) {
                int v1 = subset[source_idx];
            
                // construct the X set based on v1 (everything minus v1)
                vector<int> X;
                for (int i = 0; i < P___; ++i) {
                    if (subset[i] != v1){
                        X.emplace_back(subset[i]);
                    }
                }
                
                // construct estimator using v1
                // step 1: estimate the number of common neighbors among all vertices in subset
                
                long double f1 = 0;
                for(auto nb: g.neighbor[v1]){
                    long double fv1 = 1.0; 

                    // computing a product
                    for(auto xx : X){
                        long double A_ = g2.has(nb, xx) ? 1 : 0 ; 
                        A_ = (A_-p) / (1-2*p); 
                        fv1 = fv1 * A_;
                    }
                    f1 += fv1; 
                }
                f1 += stats::rlaplace(0.0, (pow(gamma__,X.size())/Eps2), engine);
                
                // Debug output: show f1 estimate from each vertex
                debug_file << "  vertex " << v1 << " f1=" << f1 << std::endl;
                debug_file.flush();
                
                // step 2: estimate the variance of f1.
                long double esti_var_f = 0; 
                long double theta = p * (1-p) / pow(1-2*p, 2);

                long double esti_var_f_noisy = 0; 
                int X_size = X.size();
                for (int mask = 0; mask < (1 << X_size) - 1; ++mask) {
                    // Include all subsets, from 0 to (1 << X_size) - 1
                    vector<int> Y;
                    for (int i = 0; i < X_size; ++i) {
                        if (mask & (1 << i)) {
                            Y.push_back(X[i]);
                        }
                    }
                    // processing subset Y:
                    // estimate the number of common neighbors among Y \cup v1.
                    long double f1Y = 0; 
                    if(Y.size()==0){
                        f1Y = 1; 
                    }else{
                        for(auto nb: g.neighbor[v1]){
                            long double fv1 = 1.0; 
                            for(auto xx : Y){
                                long double A_ = g2.has(nb, xx) ? 1 : 0 ; 
                                A_ = (A_-p) / (1-2*p);
                                fv1 = fv1 * A_; 
                            }
                            f1Y += fv1; 
                        }
                    }
                    esti_var_f_noisy += pow(theta, P___-1-Y.size()) * f1Y;
                }
                esti_var_f_noisy += 2 *pow(gamma__,2*P___-2) / pow(Eps2, 2);
                esti_var_f = esti_var_f_noisy;

                // Each vertex computes its own local result
                long double vertex_local_res = compute_local_res(K___, f1, esti_var_f);
                f_estimates[source_idx] = vertex_local_res;
                f1_estimates[source_idx] = f1;
            }
            
            // Multi-source: average all local results
            long double aggregated_local_res = 0;
            for(int i = 0; i < P___; i++) {
                aggregated_local_res += f_estimates[i];
            }
            aggregated_local_res /= P___;
            
            local_res = aggregated_local_res;
            
            // Calculate and show average of f1 estimates
            long double avg_f1 = 0;
            for(int i = 0; i < P___; i++) {
                avg_f1 += f1_estimates[i];
            }
            avg_f1 /= P___;
            debug_file << "  avg_f1=" << avg_f1 << std::endl;
            
            // Empty line to separate p-tuples
            debug_file << std::endl;
            debug_file.flush();
        } else {
            // Single-source estimator: use only the minimum degree vertex
            // construct the X set based on v1 (everything minus v1)
            vector<int> X;
            for (int i = 0; i < P___; ++i) {
                if (subset[i] != v1){
                    X.emplace_back(subset[i]);
                }
            }
            
            // construct estimator using v1
            long double f1 = 0;
            for(auto nb: g.neighbor[v1]){
                long double fv1 = 1.0; 

                // computing a product
                for(auto xx : X){
                    long double A_ = g2.has(nb, xx) ? 1 : 0 ; 
                    A_ = (A_-p) / (1-2*p); 
                    fv1 = fv1 * A_;
                }
                f1 += fv1; 
            }
            f1 += stats::rlaplace(0.0, (pow(gamma__,X.size())/Eps2), engine);
            
            // step 2: estimate the variance of f1.
            long double esti_var_f = 0; 
            long double theta = p * (1-p) / pow(1-2*p, 2);

            long double esti_var_f_noisy = 0; 
            int X_size = X.size();
            for (int mask = 0; mask < (1 << X_size) - 1; ++mask) {
                vector<int> Y;
                for (int i = 0; i < X_size; ++i) {
                    if (mask & (1 << i)) {
                        Y.push_back(X[i]);
                    }
                }
                long double f1Y = 0; 
                if(Y.size()==0){
                    f1Y = 1; 
                }else{
                    for(auto nb: g.neighbor[v1]){
                        long double fv1 = 1.0; 
                        for(auto xx : Y){
                            long double A_ = g2.has(nb, xx) ? 1 : 0 ; 
                            A_ = (A_-p) / (1-2*p);
                            fv1 = fv1 * A_; 
                        }
                        f1Y += fv1; 
                    }
                }
                esti_var_f_noisy += pow(theta, P___-1-Y.size()) * f1Y;
            }
            esti_var_f_noisy += 2 *pow(gamma__,2*P___-2) / pow(Eps2, 2);
            esti_var_f = esti_var_f_noisy;

            local_res = compute_local_res(K___, f1, esti_var_f);
        }

        #pragma omp critical
        res___ += local_res;
        // res___ += aggregated_local_res;
    }
    }

    double t2 = omp_get_wtime();
    cout << "time = " << t2 - t1 << endl;
    
    return res___ / sample_fraction;
}

long double compute_local_res(int K, long double f_u_w, long double esti_var_f) {
    long double f = f_u_w;
    long double s2 = esti_var_f;  // σ²
    long double f2 = f * f;
    long double f3 = f2 * f;
    long double f4 = f3 * f;
    long double f5 = f4 * f;
    long double f6 = f5 * f;
    long double f7 = f6 * f;
    long double f8 = f7 * f;
    long double f9 = f8 * f;
    long double f10 = f9 * f;
    
    long double s4 = s2 * s2;
    long double s6 = s4 * s2;
    long double s8 = s4 * s4;
    
    if (K == 1) {
        return f;
    }

    if (K == 2) {
        return f2 / 2.0 - 0.5 * f - s2 / 2.0;
    }

        if (K == 3) {
        return f3 / 6.0 - 0.5 * f2 - f * s2 / 2.0 + 0.333333333333333 * f + 0.5 * s2;
    }

        if (K == 4) {
        return f4 / 24.0 - 0.25 * f3 - f2 * s2 / 4.0 + 0.458333333333333 * f2 
               + 0.75 * f * s2 - 0.25 * f + s4 / 8.0 - 0.458333333333333 * s2;
    }

        if (K == 5) {
        return f5 / 120.0 - 0.0833333333333333 * f4 - f3 * s2 / 12.0 + 0.291666666666667 * f3 
               + 0.5 * f2 * s2 - 0.416666666666667 * f2 + f * s4 / 8.0 - 0.875 * f * s2 
               + 0.2 * f - 0.25 * s4 + 0.416666666666667 * s2;
    }

        if (K == 6) {
        return f6 / 720.0 - 0.0208333333333333 * f5 - f4 * s2 / 48.0 + 0.118055555555556 * f4 
               + 0.208333333333333 * f3 * s2 - 0.3125 * f3 + f2 * s4 / 16.0 
               - 0.708333333333333 * f2 * s2 + 0.380555555555556 * f2 
               - 0.3125 * f * s4 + 0.9375 * f * s2 - 0.166666666666667 * f 
               - s6 / 48.0 + 0.354166666666667 * s4 - 0.380555555555556 * s2;
    }
    
        if (K == 7) {
        return f7 / 5040.0 - 0.00416666666666667 * f6 - f5 * s2 / 240.0 + 0.0347222222222222 * f5 
               + 0.0625 * f4 * s2 - 0.145833333333333 * f4 + f3 * s4 / 48.0 
               - 0.347222222222222 * f3 * s2 + 0.322222222222222 * f3 
               - 0.1875 * f2 * s4 + 0.875 * f2 * s2 - 0.35 * f2 
               - f * s6 / 48.0 + 0.520833333333333 * f * s4 - 0.966666666666667 * f * s2 
               + 0.142857142857143 * f + 0.0625 * s6 - 0.4375 * s4 + 0.35 * s2;
    }
    
        if (K == 8) {
        return f8 / 40320.0 - 0.000694444444444444 * f7 - f6 * s2 / 1440.0 
               + 0.00798611111111111 * f6 + 0.0145833333333333 * f5 * s2 
               - 0.0486111111111111 * f5 + f4 * s4 / 192.0 - 0.119791666666667 * f4 * s2 
               + 0.167881944444444 * f4 - 0.0729166666666667 * f3 * s4 
               + 0.486111111111111 * f3 * s2 - 0.325694444444444 * f3 
               - f2 * s6 / 96.0 + 0.359375 * f2 * s4 - 1.00729166666667 * f2 * s2 
               + 0.324107142857143 * f2 + 0.072916666666667 * f * s6 
               - 0.729166666666667 * f * s4 + 0.977083333333333 * f * s2 - 0.125 * f 
               + s8 / 384.0 - 0.119791666666667 * s6 + 0.503645833333333 * s4 
               - 0.324107142857143 * s2;
    }
    
        if (K == 9) {
        return f9 / 362880.0 - 9.92063492063492e-5 * f8 - f7 * s2 / 10080.0 
               + 0.00150462962962963 * f7 + 0.00277777777777778 * f6 * s2 - 0.0125 * f6 
               + f5 * s4 / 960.0 - 0.0315972222222222 * f5 * s2 + 0.0618634259259259 * f5 
               - 0.0208333333333333 * f4 * s4 + 0.1875 * f4 * s2 - 0.185416666666667 * f4 
               - f3 * s6 / 288.0 + 0.157986111111111 * f3 * s4 - 0.618634259259259 * f3 * s2 
               + 0.325518077601411 * f3 + 0.0416666666666667 * f2 * s6 
               - 0.5625 * f2 * s4 + 1.1125 * f2 * s2 - 0.301984126984127 * f2 
               + f * s8 / 384.0 - 0.157986111111111 * f * s6 + 0.927951388888889 * f * s4 
               - 0.976554232804233 * f * s2 + 0.111111111111111 * f 
               - 0.0104166666666667 * s8 + 0.1875 * s6 - 0.55625 * s4 
               + 0.301984126984127 * s2;
    }

    if (K == 10) {
        long double s10 = s8 * s2;
        return f10 / 3628800.0 - 1.24007936507936e-5 * f9 - f8 * s2 / 80640.0 
               + 0.000239748677248677 * f8 + 0.000446428571428571 * f7 * s2 
               - 0.00260416666666667 * f7 + f6 * s4 / 5760.0 - 0.00671296296296296 * f6 * s2 
               + 0.0174363425925926 * f6 - 0.0046875 * f5 * s4 + 0.0546875 * f5 * s2 
               - 0.07421875 * f5 - f4 * s6 / 1152.0 + 0.0503472222222222 * f4 * s4 
               - 0.261545138888889 * f4 * s2 + 0.199426807760141 * f4 
               + 0.015625 * f3 * s6 - 0.2734375 * f3 * s4 + 0.7421875 * f3 * s2 
               - 0.323164682539683 * f3 + f2 * s8 / 768.0 - 0.100694444444444 * f2 * s6 
               + 0.784635416666667 * f2 * s4 - 1.19656084656085 * f2 * s2 
               + 0.282896825396825 * f2 - 0.01171875 * f * s8 + 0.2734375 * f * s6 
               - 1.11328125 * f * s4 + 0.969494047619048 * f * s2 - 0.1 * f 
               - s10 / 3840.0 + 0.0251736111111111 * s8 - 0.261545138888889 * s6 
               + 0.598280423280423 * s4 - 0.282896825396825 * s2;
    }

    // Unsupported K
    return 0;
}

double locally_compute_f_given_q_and_x(int q, int x, BiGraph& g, BiGraph& g2) {

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

        // only consider when priority nb < priority q

        // cout<<"nb of q: "<<nb<<endl;
		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
	}
    // cout<<endl;

	// if x is not a neighbor of q, then 
	// if(g2.is_bipartite){
		// in the case of bipartite graphs
    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }
	// }else{
	// 	// in the case of general graphs 
	// 	if(x_is_a_nb_of_q){
	// 		assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x-1);
	// 	}else{
	// 		assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x); 
	// 	}
	// }
	// locally the degree of q is known. 
	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; 

	return res;
}

double locally_compute_f_given_q_and_x_ad_hoc(int q, int x, BiGraph& g, BiGraph& g2) {

	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;
    // cout<<"x = "<<x <<endl;q
    // cout<<"visiting nb of q:"<< q<<endl;
	for(auto nb: g.neighbor[q]){

		if(nb == x){
			x_is_a_nb_of_q = true;
		}

        // need to check (nb, x) \in g2.
        if (!g2.has_computed(nb, x)){
            randomized_response_single_bit(nb, x, g, g2);
        }

        // if(g2.edge_vector[min(nb,x)][max(nb,x)]){
		// 	Nx_cap_Nq_minus_x++;
		// }else{
		// 	Nq_minus_Nx_minus_x++;
		// }
        unsigned int smaller = (nb < x) ? nb : x;
        unsigned int larger = (nb < x) ? x : nb;

        if (g2.edge_vector[smaller][larger]) {
            Nx_cap_Nq_minus_x++;
        } else {
            Nq_minus_Nx_minus_x++;
        }
	}


    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; // why without noise it is even bigger? 

	return res;
}


double locally_compute_f_given_q_and_x_two_graphs(int q, int x, BiGraph& g, BiGraph& g2, BiGraph& g3) {

	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;

	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;

    double Nx_cap_Nq_minus_x_dup =0, Nq_minus_Nx_minus_x_dup =0;

	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
		if(g3.has(nb, x)){
			Nx_cap_Nq_minus_x_dup++;
		}else{
			Nq_minus_Nx_minus_x_dup++;
		}
	}

    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

    double res2 = Nx_cap_Nq_minus_x_dup * gamma__ + Nq_minus_Nx_minus_x_dup * (-p)/(1-2*p); 

    res = (res + res2)/2;

    // the GS is the same.
	res += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	return res;
}


// count the common neighbor of q and x 
// where prority is less than q 
double locally_compute_f_given_q_and_x_vp(int q, int x, BiGraph& g, BiGraph& g2, int& res__) {

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

        // only consider when priority nb < priority q
        if(g.prio[q] <= g.prio[nb] ) 
            continue;

        assert(g.prio[q] > g.prio[nb]);

        // cout<<"nb of q: "<<nb<<endl;
		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
	}
    // cout<<endl;

	// if x is not a neighbor of q, then 
	// if(g2.is_bipartite){
		// in the case of bipartite graphs
    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	// locally the degree of q is known. 
	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

    res__ = Nx_cap_Nq_minus_x + Nq_minus_Nx_minus_x;

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; 

	return res;
}

double locally_compute_f_given_q_and_x_vp_2(int q, int x, BiGraph& g, BiGraph& g2, int& res__) {

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

        // only consider when priority nb < priority q
        if(g.prio[x] <= g.prio[nb] ) 
            continue;

        assert(g.prio[x] > g.prio[nb]);

        // cout<<"nb of q: "<<nb<<endl;
		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
	}
    // cout<<endl;

	// if x is not a neighbor of q, then 
	// if(g2.is_bipartite){
		// in the case of bipartite graphs
    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	// locally the degree of q is known. 
	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

    res__ = Nx_cap_Nq_minus_x + Nq_minus_Nx_minus_x;

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; 

	return res;
}

// this works

long double binomial(int n, int k) {
    if (n < 0 || k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    k = std::min(k, n - k); // Optimize by using smaller k
    long double res = 1.0;
    for (int i = 0; i < k; ++i) {
        res *= static_cast<long double>(n - i) / (i + 1);
    }
    return res;
}



// one-round biclique counting: 
// _switch = btf:0, cate:1, biclique:2, quasi-biclique: 3.
// Function to test f distribution for Gaussian assumption
void test_f_distribution_p2(string dataset, int num_samples) {
    cout << "=== Testing f_avg distribution for Gaussian assumption ===" << endl;
    cout << "Method: 10 fixed pairs × " << num_samples << " noise iterations" << endl;
    
    // Load dataset
    BiGraph g;
    string filename = "/data/yizhangh/bidata/" + dataset;
    g.loadGraph(filename);
    cout << "Loaded dataset: " << dataset << " with " << g.num_nodes() << " nodes" << endl;
    
    // Set up parameters like the existing function
    Eps = 1.0;
    Eps0 = Eps * 0.05;
    Eps1 = Eps * 0.6;
    Eps2 = Eps * 0.35;
    p = 1.0 / (exp(Eps1) + 1.0);
    gamma__ = (1-p) / (1-2*p);
    
    // Initialize degree estimates
    vector<long double> deg_estis;
    deg_estis.resize(g.num_nodes());
    for(int i=0; i<g.num_nodes(); i++){
        deg_estis[i] = g.degree[i] + stats::rlaplace(0.0, 1/(Eps0), engine);
    }
    
    // Find vertex pairs with the most common neighbors
    vector<pair<int,int>> fixed_pairs;
    vector<pair<pair<int,int>, int>> all_pairs_with_cn;  // Store (pair, cn_count)
    
    cout << "Searching for vertex pairs with the most common neighbors..." << endl;
    
    // First, find all pairs and their common neighbor counts
    for (int u = 0; u < g.num_v1; u++) {
        for (int w = u + 1; w < g.num_v1; w++) {
            // Count common neighbors for this pair
            int common_neighbors = 0;
            for (int v = 0; v < g.num_v2; v++) {
                if (g.has(u, v) && g.has(w, v)) {
                    common_neighbors++;
                }
            }
            
            if (common_neighbors > 0) {
                all_pairs_with_cn.push_back({{u, w}, common_neighbors});
            }
        }
    }
    
    // Sort by common neighbor count (descending)
    sort(all_pairs_with_cn.begin(), all_pairs_with_cn.end(), 
         [](const pair<pair<int,int>, int>& a, const pair<pair<int,int>, int>& b) {
             return a.second > b.second;
         });
    
    cout << "Found " << all_pairs_with_cn.size() << " pairs with common neighbors" << endl;
    cout << "Top 20 pairs by common neighbor count:" << endl;
    for (int i = 0; i < min(20, (int)all_pairs_with_cn.size()); i++) {
        cout << "  Rank " << i+1 << ": vertices (" << all_pairs_with_cn[i].first.first 
             << "," << all_pairs_with_cn[i].first.second << ") with " 
             << all_pairs_with_cn[i].second << " common neighbors" << endl;
    }
    
    // Select top 10 pairs
    for (int i = 0; i < min(10, (int)all_pairs_with_cn.size()); i++) {
        fixed_pairs.push_back(all_pairs_with_cn[i].first);
        cout << "Selected Pair " << i << ": vertices (" << all_pairs_with_cn[i].first.first 
             << "," << all_pairs_with_cn[i].first.second << ") with " 
             << all_pairs_with_cn[i].second << " common neighbors" << endl;
    }
    
    cout << "Selected 10 fixed pairs for p=2 analysis" << endl;
    
    // For each pair, run num_samples iterations with different noise
    for (int pair_idx = 0; pair_idx < 10; pair_idx++) {
        int u = fixed_pairs[pair_idx].first;
        int w = fixed_pairs[pair_idx].second;
        
        // Compute real number of common neighbors for this pair
        int real_common_neighbors = 0;
        for (int v = 0; v < g.num_v2; v++) {
            if (g.has(u, v) && g.has(w, v)) {
                real_common_neighbors++;
            }
        }
        
        vector<double> f_u_vals, f_w_vals;
        
        for (int iter = 0; iter < num_samples; iter++) {
            // Create noisy graph with different seed each time
            BiGraph g2(g);
            construct_noisy_graph(g, g2, iter);
            
            // Compute f_u and f_w using the same approach as the working function
            long double f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
            long double f_w_u = locally_compute_f_given_q_and_x(w, u, g, g2);
            
            f_u_vals.push_back(f_u_w);
            f_w_vals.push_back(f_w_u);
        }
        
        // Output for this pair
        cout << "Pair_" << pair_idx << "_real_common_neighbors: " << real_common_neighbors << endl;
        cout << "Pair_" << pair_idx << "_degree_u: " << g.degree[u] << endl;
        cout << "Pair_" << pair_idx << "_degree_w: " << g.degree[w] << endl;
        cout << "Pair_" << pair_idx << "_p2_f_u: ";
        for (auto v : f_u_vals) cout << v << " ";
        cout << endl;
        
        cout << "Pair_" << pair_idx << "_p2_f_w: ";
        for (auto v : f_w_vals) cout << v << " ";
        cout << endl;
        
    }
    
    cout << "\n=== Data collection complete ===" << endl;
    cout << "Generated: 10 pairs × " << num_samples << " = " << (10 * num_samples) << " samples for p=2" << endl;
    cout << "This isolates the noise distribution from graph structure" << endl;
}

long double one_round_biclique(BiGraph& g, unsigned long seed, 
    int p__, int q__){

	BiGraph g2(g); 
    
	construct_noisy_graph(g, g2, seed);

    cout<<"p__ = "<<p__ <<endl;
    cout<<"q__ = "<<q__ <<endl;


    std::vector<int> U, L; 
    for (int i = 0; i < g.num_v1; ++i) U.push_back(i); 
    
    for (int i = g.num_v1; i < g.num_nodes(); ++i) L.push_back(i);
    

    // Generate all combinations of p vertices from U
    generate_combinations(U, p__, up_options);
    generate_combinations(L, q__, lo_options);
	
	cout<<"counting biclique on noisy graph"<<endl;

    
    std::vector<long double> m__(p__ * q__ + 1, 0);
    // Convert adjacency lists to unordered_set for O(1) edge lookup
    std::vector<std::unordered_set<int>> adj(g2.neighbor.size());
    for (size_t i = 0; i < g2.neighbor.size(); i++) {
        
        adj[i] = std::unordered_set<int>( g2.neighbor[i].begin(), g2.neighbor[i].end());

    }


    std::cout << "Old way: Counting mi numbers\n";
    #pragma omp parallel
    {
        std::vector<long double> local_m(p__ * q__ + 1, 0);  // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t up_idx = 0; up_idx < up_options.size(); up_idx++) {
            for (size_t lo_idx = 0; lo_idx < lo_options.size(); lo_idx++) {
                const auto& xxx = up_options[up_idx];
                const auto& yyy = lo_options[lo_idx];

                int num_edges = 0;
                for (int u : xxx) {
                    for (int v : yyy) {
                        if (adj[u].count(v)) num_edges++; // O(1) edge lookup
                    }
                }
                local_m[num_edges]++;
            }
        }
        // Reduce results
        #pragma omp critical
        for (size_t i = 0; i <= p__ * q__; i++) { // Start from 0
            #pragma omp atomic
            m__[i] += local_m[i];
        }
    }
    

    
    /*
    vector<long double> m__(p__ * q__ + 1, 0);
	cout<<"counting mi numbers"<<endl;
	#pragma omp parallel for collapse(2)
	for(auto xxx:up_options){ 
		for(auto yyy:lo_options){

            // for each motif, check how many incident edges are there
			int num_edges = 0;
			for(auto i:xxx){
				for(auto j:yyy){
                    // they all need to be checked using g2. 
					if(std::find(g2.neighbor[i].begin(), g2.neighbor[i].end(), j) != g2.neighbor[i].end() ){
						num_edges++;
					}
				}
			}
			#pragma omp atomic
			m__[num_edges]++;
		}
	}
    */
    

	for(int i=0;i<m__.size();i++){
		cout<<"edge = "<<i<<" num = "<<m__[i]<<endl;
	}



	long double res = 0, mu = exp(Eps); 
	
	// if(_switch==2){
    // hard to obtain the other motif counts
    for(int i=0;i<m__.size();i++){
        res +=  power(-mu, i) * m__[i]; 
    }
    res /= power(1-mu, p__*q__);
    naive_estis[iteration] = m__[m__.size()-1];
	// }

    /*
	if(_switch==3){

		res += -6*mu* m__[0];

		res += (5*power(mu,2) + 1) * m__[1];

		res += -2*mu*(2*power(mu,2) + 1) * m__[2];

		res += 3*power(mu,2)*(power(mu,2) + 1) * m__[3];

		res += -2*power(mu,3)*(power(mu,2) + 2) * m__[4];

		res += power(mu,4)*(power(mu,2) + 5) * m__[5];

		res += -6*power(mu,5) * m__[6];

		res /= power(mu-1, 6);

		naive_estis[iteration] = m__[m__.size()-2];
	}
    */
	// double t2 = omp_get_wtime();

	// record time elapsed.
	// RR_time += t1-t0;
	// server_side_time += t2-t1;
	//
    cout<<"naive esti = "<<naive_estis[iteration] <<endl;
	return res;
}


long double one_round_biclique_2_K(BiGraph& g, int K, unsigned long seed) {
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);

    int p__ = 2;
    int q__ = K;
    std::cout << "p__ = " << p__ << "\n";
    std::cout << "q__ = " << q__ << "\n";

    // Get upper and lower vertex indices
    std::vector<int> U, L;
    for (int i = 0; i < g.num_v1; ++i) U.push_back(i);
    for (int i = g.num_v1; i < g.num_nodes(); ++i) L.push_back(i);
    int n1 = U.size();
    int n2 = L.size();

    std::cout << "Counting biclique on noisy graph\n";

    // Convert adjacency lists to unordered_set for O(1) lookup
    std::vector<std::unordered_set<int>> adj(g2.neighbor.size());
    for (size_t i = 0; i < g2.neighbor.size(); ++i) {
        adj[i] = std::unordered_set<int>(g2.neighbor[i].begin(), g2.neighbor[i].end());
    }

    // Array to store motif counts (B_i for i = 0 to 2*K)
    std::vector<long double> m__(2 * K + 1, 0.0);

    std::cout << "Counting Bi numbers (optimized)\n";
    #pragma omp parallel
    {
        std::vector<long double> local_m(2 * K + 1, 0.0); // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t u1_idx = 0; u1_idx < n1; ++u1_idx) {
            for (size_t u2_idx = u1_idx + 1; u2_idx < n1; ++u2_idx) {
                int u1 = U[u1_idx];
                int u2 = U[u2_idx];

                // Compute s2 = |N(u1) ∩ N(u2)|
                int s2 = 0;
                for (int v : g2.neighbor[u1]) {
                    if (adj[u2].count(v)) ++s2;
                }

                // Compute degrees
                int deg_u1 = g2.neighbor[u1].size();
                int deg_u2 = g2.neighbor[u2].size();

                // Compute s1 = deg(u1) + deg(u2) - 2 * s2
                int s1 = deg_u1 + deg_u2 - 2 * s2;

                // Compute s0 = n2 - (deg(u1) + deg(u2) - s2)
                int s0 = n2 - (deg_u1 + deg_u2 - s2);

                // Compute contributions to m__[i] for i = 0 to 2*K
                for (int c = 0; c <= K; ++c) {
                    for (int b = 0; b <= K - c; ++b) {
                        int a = K - b - c;
                        if (a >= 0) {
                            int i = 2 * c + b;
                            long double contrib = binomial(s0, a) * binomial(s1, b) * binomial(s2, c);
                            local_m[i] += contrib;
                        }
                    }
                }
            }
        }

        // Reduce results
        #pragma omp critical
        for (size_t i = 0; i <= 2 * K; ++i) {
            #pragma omp atomic
            m__[i] += local_m[i];
        }
    }

    cout<<"motif count distribution"<<endl;
    for (size_t i = 0; i <= 2 * K; ++i) {
        cout<<"i = "<<i<<" mi = "<< m__[i] <<endl;
    }

    // Output motif counts and verify sum
    long double sum___ = std::accumulate(m__.begin(), m__.end(), 0.0);
    long double target = binomial(n1, 2) * binomial(n2, K);
    std::cout << "Sum of motif counts: " << sum___ << "\n";
    std::cout << "Expected total: " << target << "\n";
    if (std::abs(sum___ - target) > 1e-6 * target) {
        // when this happens, there is the issue of overflow
        std::cout << "Sum of motif counts: " << sum___ << "\n";
        std::cout << "Expected total: " << target << "\n";
        std::cerr << "Warning: Sum of motif counts does not match expected total.\n";
    }

    // Compute unbiased estimate using Theorem 2
    long double res = 0, mu = std::exp(Eps);
    for (size_t i = 0; i <= 2 * K; ++i) {
        res += std::pow(-mu, i) * m__[i];
    }
    res /= std::pow(1 - mu, 2 * K);
    naive_estis[iteration] = m__[2 * K];

    std::cout << "naive esti = " << naive_estis[iteration] << "\n";
    return res;
}

// Batch version of one_round_biclique_2_K that handles Q = [4, 5, 6, 7, 8, 9, 10]
// This function builds the noisy graph once, then estimates biclique counts for all Q values
std::vector<long double> one_round_biclique_2_K_batch(BiGraph& g, unsigned long seed) {
    cout << "Running one_round_biclique_2_K batch algorithm..." << endl;
    
    // Build noisy graph once
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);

    int p__ = 2;
    std::cout << "p__ = " << p__ << "\n";

    // Get upper and lower vertex indices
    std::vector<int> U, L;
    for (int i = 0; i < g.num_v1; ++i) U.push_back(i);
    for (int i = g.num_v1; i < g.num_nodes(); ++i) L.push_back(i);
    int n1 = U.size();
    int n2 = L.size();

    std::cout << "Counting biclique on noisy graph\n";

    // Convert adjacency lists to unordered_set for O(1) lookup
    std::vector<std::unordered_set<int>> adj(g2.neighbor.size());
    for (size_t i = 0; i < g2.neighbor.size(); ++i) {
        adj[i] = std::unordered_set<int>(g2.neighbor[i].begin(), g2.neighbor[i].end());
    }

    // Initialize results for Q = [4, 5, 6, 7, 8, 9, 10]
    std::vector<long double> results(7, 0.0);
    
    // For each Q value, compute the estimate
    for (int q_idx = 0; q_idx < 7; q_idx++) {
        int K = q_idx + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
        std::cout << "Processing Q = " << K << "\n";

    // Array to store motif counts (B_i for i = 0 to 2*K)
    std::vector<long double> m__(2 * K + 1, 0.0);

        std::cout << "Counting Bi numbers (optimized) for Q=" << K << "\n";
    #pragma omp parallel
    {
        std::vector<long double> local_m(2 * K + 1, 0.0); // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t u1_idx = 0; u1_idx < n1; ++u1_idx) {
            for (size_t u2_idx = u1_idx + 1; u2_idx < n1; ++u2_idx) {
                int u1 = U[u1_idx];
                int u2 = U[u2_idx];

                // Compute s2 = |N(u1) ∩ N(u2)|
                int s2 = 0;
                for (int v : g2.neighbor[u1]) {
                    if (adj[u2].count(v)) ++s2;
                }

                // Compute degrees
                int deg_u1 = g2.neighbor[u1].size();
                int deg_u2 = g2.neighbor[u2].size();

                // Compute s1 = deg(u1) + deg(u2) - 2 * s2
                int s1 = deg_u1 + deg_u2 - 2 * s2;

                // Compute s0 = n2 - (deg(u1) + deg(u2) - s2)
                int s0 = n2 - (deg_u1 + deg_u2 - s2);

                // Compute contributions to m__[i] for i = 0 to 2*K
                for (int c = 0; c <= K; ++c) {
                    for (int b = 0; b <= K - c; ++b) {
                        int a = K - b - c;
                        if (a >= 0) {
                            int i = 2 * c + b;
                            long double contrib = binomial(s0, a) * binomial(s1, b) * binomial(s2, c);
                            local_m[i] += contrib;
                        }
                    }
                }
            }
        }

        // Reduce results
        #pragma omp critical
        for (size_t i = 0; i <= 2 * K; ++i) {
            #pragma omp atomic
            m__[i] += local_m[i];
        }
    }

    // Compute unbiased estimate using Theorem 2
    long double res = 0, mu = std::exp(Eps);
    for (size_t i = 0; i <= 2 * K; ++i) {
        res += std::pow(-mu, i) * m__[i];
    }
    res /= std::pow(1 - mu, 2 * K);
        
        results[q_idx] = res;
        std::cout << "Q" << K << " one_round estimate = " << res << "\n";
    }
    
    cout << "One_round batch results for Q=[4,5,6,7,8,9,10]: ";
    for (int q_idx = 0; q_idx < 7; q_idx++) {
        cout << "Q" << (q_idx + 4) << "=" << results[q_idx] << " ";
    }
    cout << endl;
    
    return results;
}



// Naive biclique count
long double naive_biclique(BiGraph& g, unsigned long seed, 
    int p__, int q__){

	BiGraph g2(g); 
    
    long double res = 0; 

	construct_noisy_graph(g, g2, seed);

    // if(p__ == 2 && q__ ==2){
    //     res = BFC_EVP(g2);
    //     cout<<"btf res = "<<res<<endl;
    //     return res;
    // }

    // only do this when larger biclique
    biGraph convertedGraph = convertBiGraphTobiGraph(g2);

    cout << "Converted graph: n1=" << convertedGraph.n1 << ", n2=" << convertedGraph.n2 << ", m=" << convertedGraph.m << std::endl;
    
    // Create and use BCListPlusPlus
    BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, p__, q__);

    res = counter->exactCount();

    cout<<"res = "<<res<<endl;

	return res;
}


// Naive biclique counting with vertex sampling for speed
long double naive_biclique_with_vertex_sampling(BiGraph& g, unsigned long seed, 
    int p__, int q__, double sampling_ratio, int num_samples){

    BiGraph g2(g);
    
    long double res = 0; 

	construct_noisy_graph(g, g2, seed);

    cout << "Original noisy graph: n1=" << g2.num_v1 << ", n2=" << g2.num_v2 << ", edges=" << g2.num_edges << endl;

    // Apply vertex sampling if sampling_ratio < 1.0
    if (sampling_ratio < 1.0) {
        cout << "Applying vertex sampling with ratio: " << sampling_ratio << " (" << num_samples << " samples)" << endl;
        
        long double total_estimate = 0.0;
        
        for (int sample = 0; sample < num_samples; sample++) {
            // Create vertex-sampled graph
            BiGraph g2_sampled;
            g2_sampled.init(g2.num_v1, g2.num_v2);
            
            // Initialize random number generator with different seed for each sample
            init_genrand(seed + 12345 + sample * 1000);
            
            // Sample vertices from both partitions
            vector<bool> sampled_v1(g2.num_v1, false);
            vector<bool> sampled_v2(g2.num_v2, false);
            
            int sampled_v1_count = 0;
            int sampled_v2_count = 0;
            
            // Sample vertices from partition 1 (upper)
            for (int u = 0; u < g2.num_v1; u++) {
                if (genrand_real2() < sampling_ratio) {
                    sampled_v1[u] = true;
                    sampled_v1_count++;
                }
            }
            
            // Sample vertices from partition 2 (lower)
            for (int v = 0; v < g2.num_v2; v++) {
                if (genrand_real2() < sampling_ratio) {
                    sampled_v2[v] = true;
                    sampled_v2_count++;
                }
            }
            
            // Add edges between sampled vertices
            int sampled_edges = 0;
            for (int u = 0; u < g2.num_v1; u++) {
                if (sampled_v1[u]) {
                    for (auto v : g2.neighbor[u]) {
                        int v_lower = v - g2.num_v1;  // Convert to lower partition index
                        if (sampled_v2[v_lower]) {
                            g2_sampled.addEdge(u, v);
                            sampled_edges++;
                        }
                    }
                }
            }
            
            if (sample == 0) {  // Only print details for first sample
                cout << "Sample " << (sample+1) << " - Sampled vertices: " << sampled_v1_count << " (upper), " << sampled_v2_count << " (lower)" << endl;
                cout << "Sample " << (sample+1) << " - Sampled graph edges: " << sampled_edges << endl;
                cout << "Sample " << (sample+1) << " - Vertex sampling ratio achieved: " << (double)sampled_v1_count / g2.num_v1 << " (upper), " << (double)sampled_v2_count / g2.num_v2 << " (lower)" << endl;
            }
            
            // Use sampled graph for counting
            biGraph convertedGraph = convertBiGraphTobiGraph(g2_sampled);
            
            // Count bicliques on sparse graph
            BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, p__, q__);
            long double sparse_count = counter->exactCount();
            delete counter;
            
            // Unbiased estimate: scale by 1/sampling_ratio^(p+q) (for vertex sampling)
            long double scaling_factor = pow(sampling_ratio, p__ + q__);
            long double sample_estimate = sparse_count / scaling_factor;
            
            total_estimate += sample_estimate;
            
            if (sample == 0) {  // Only print details for first sample
                cout << "Sample " << (sample+1) << " - Sparse graph biclique count: " << sparse_count << endl;
                cout << "Sample " << (sample+1) << " - Scaling factor: " << scaling_factor << endl;
                cout << "Sample " << (sample+1) << " - Unbiased estimate: " << sample_estimate << endl;
            }
        }
        
        // Average the estimates
        res = total_estimate / num_samples;
        cout << "Average of " << num_samples << " samples: " << res << endl;
        
    } else {
        // No sampling - use original method
        biGraph convertedGraph = convertBiGraphTobiGraph(g2);
        cout << "Converted graph: n1=" << convertedGraph.n1 << ", n2=" << convertedGraph.n2 << ", m=" << convertedGraph.m << std::endl;
        
        BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, p__, q__);
        res = counter->exactCount();
        delete counter;
        
        cout << "Exact count: " << res << endl;
    }

	return res;
}
// todo: implement the vertex-priority-based wedge butterfly counting 

// TODO: (1) use estimated priority instead
// (2) use estimated value of number of priority-obeying neighbors
// (3) if we do not care about degree, just use vertex id, what will happen?

void fetch_or_compute_biclique_count(int P___, int K___, 
    string dataset, BiGraph& g){

    sqlite3* db;
    if (sqlite3_open("../biclq_counts.db", &db) != SQLITE_OK) {
        std::cerr << "Error opening database: " << sqlite3_errmsg(db) << std::endl;
        exit(1);
    }
    // Dataset, p, and q values to filter
    size_t found = dataset.find_last_of("/\\");  // Find the last slash
    std::string dataset_to_find = dataset.substr(found + 1);  // Extract part after the last slash

    // std::string dataset_to_find = "unicode";  // Example dataset
    int p_to_find = P___;                   // Example p value
    int q_to_find = K___;                   // Example q value

    // Query to retrieve one row based on dataset, p, and q
    const char* sql = "SELECT dataset, count FROM pqbiclique_counts WHERE dataset = ? AND p = ? AND q = ?;";
    sqlite3_stmt* stmt;

    // Prepare statement
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        std::cerr << "Error preparing statement: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        exit(1);
    }

    // Bind parameters
    sqlite3_bind_text(stmt, 1, dataset_to_find.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int(stmt, 2, p_to_find);
    sqlite3_bind_int(stmt, 3, q_to_find);

    // Execute the query
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        std::string dataset = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        
        // Get the column type to handle both INTEGER and REAL values
        int column_type = sqlite3_column_type(stmt, 1);
        long double count_value;
        
        if (column_type == SQLITE_INTEGER) {
            count_value = static_cast<long double>(sqlite3_column_int64(stmt, 1));
        } else if (column_type == SQLITE_FLOAT) {
            count_value = static_cast<long double>(sqlite3_column_double(stmt, 1));
        } else {
            std::cerr << "Error: Unexpected column type for count" << std::endl;
            sqlite3_close(db);
            exit(1);
        }
        
        // Store the ground truth value (use long double consistently)
        real_ld = count_value;
        real = count_value;
        
        std::cout << "Dataset: " << dataset << ", p = " << P___ << ", q = " << K___ 
                  << ", biclique count = " << std::scientific << std::setprecision(6) << count_value << std::endl;
    } else {
        std::cout << "No matching data found." << std::endl;

        // when this happens, we need to compute count ad hoc.
        biGraph convertedGraph = convertBiGraphTobiGraph(g);
        std::cout << "Converted graph: n1=" << convertedGraph.n1 
                  << ", n2=" << convertedGraph.n2 
                  << ", m=" << convertedGraph.m << std::endl;
        BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, P___, K___);
        long double computed_count = counter->exactCount();
        cout<<"cliq count = "<<computed_count<<endl;

        // Store the computed value in high precision
        real_ld = computed_count;
        
        // Also store as unsigned long long for backward compatibility (with truncation warning)
        if (computed_count > static_cast<long double>(ULLONG_MAX)) {
            std::cout << "WARNING: Computed count exceeds unsigned long long limit for dataset: " << dataset_to_find 
                      << ", p = " << P___ << ", q = " << K___ << std::endl;
            std::cout << "Original value: " << std::scientific << std::setprecision(6) << computed_count << std::endl;
            std::cout << "Will use long double precision for relative error calculation" << std::endl;
            real = computed_count;  // Store actual value in long double (no truncation)
        } else {
            real = computed_count;  // Use long double consistently
        }

        // Insert the new value into the database
        const char* insert_sql = "INSERT INTO pqbiclique_counts (dataset, p, q, count) VALUES (?, ?, ?, ?);";
        sqlite3_stmt* insert_stmt;

        // Prepare the INSERT statement
        if (sqlite3_prepare_v2(db, insert_sql, -1, &insert_stmt, nullptr) != SQLITE_OK) {
            std::cerr << "Error preparing insert statement: " << sqlite3_errmsg(db) << std::endl;
            sqlite3_close(db);
            exit(1);
        }

        // Bind parameters for the insert statement
        sqlite3_bind_text(insert_stmt, 1, dataset_to_find.c_str(), -1, SQLITE_STATIC);
        sqlite3_bind_int(insert_stmt, 2, P___);
        sqlite3_bind_int(insert_stmt, 3, K___);
        sqlite3_bind_double(insert_stmt, 4, computed_count);  // Store as REAL (double precision)

        // Execute the insert statement
        if (sqlite3_step(insert_stmt) != SQLITE_DONE) {
            std::cerr << "Error inserting data: " << sqlite3_errmsg(db) << std::endl;
        } else {
            std::cout << "Inserted new biclique count into database." << std::endl;
        }
        // Cleanup the insert statement
        sqlite3_finalize(insert_stmt);

    }
    // Cleanup
    sqlite3_finalize(stmt);
    sqlite3_close(db);
}
// biclique related code
bool BCListPlusPlus::costEstimateRaw() {
    std::default_random_engine e(10007);
    std::uniform_int_distribution<unsigned> chooseU(0, g->n1 - 1);
    std::uniform_int_distribution<unsigned> chooseV(0, g->n2 - 1);
    std::vector<uint32_t> sum(std::max(g->n1, g->n2));
    std::vector<uint32_t> tmp(std::max(g->n1, g->n2));
    int l = 0;
    uint32_t rd = 10;
    uint32_t Du = 0;
    uint32_t u = chooseU(e);
    for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
        uint32_t v = g->e1[i];
        for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
            uint32_t w = g->e2[j];

            if(w != u) {
                if(sum[w] == 0) tmp[l++] = w;
                sum[w]++;
            }
        }
    }
    for(int i = 0; i < l; i++) {
        if(sum[tmp[i]] >= (uint32_t)q) {
            Du++;
        }
        sum[tmp[i]] = 0;
    }
    l = 0;
    
    // t = rd;
    uint32_t Dv = 0;
    // while(t--) {
    uint32_t v = chooseV(e);
    // for(uint32_t v = 0; v < g->n2; v += rd) {
    for(uint32_t i = g->pV[v]; i < g->pV[v + 1]; i++) {
        uint32_t u = g->e2[i];

        for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
            uint32_t w = g->e1[j];

            if(w != v) {
                if(sum[w] == 0) tmp[l++] = w;
                sum[w]++;
            }
        }
    }
    

    for(int i = 0; i < l; i++) {
        if(sum[tmp[i]] >= (uint32_t)p) {
            Dv++;
        }
        sum[tmp[i]] = 0;
    }
    l = 0;
    
    

    double costU = pow(1.0 * Du / (g->n1 / 10), p - 2) * Du * g->maxDu;
    double costV = pow(1.0 * Dv / (g->n2 / 10), q - 2) * Dv * g->maxDv;
    // printf("cost:%.2f %.2f, du %u, dv %u\n", costU, costV, Du, Dv);

    return costU > costV;
}

bool BCListPlusPlus::costEstimate() {

    std::vector<uint32_t> sum(std::max(g->n1, g->n2));
    std::vector<uint32_t> tmp(std::max(g->n1, g->n2));
    int l = 0;

    uint32_t rd = 100;
    // uint32_t t = rd;
    uint32_t Du = 0, maxDu, rdu = 0;
    double sumDu = 0.0;
    // while(t--) {
        // uint32_t u = chooseU(e);
    for(uint32_t u = 0; u < g->n1; u += rd) {
        rdu++;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];

                if(w > u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(int i = 0; i < l; i++) {
            if(sum[tmp[i]] >= (uint32_t)q) {
                Du++;
            }
            sum[tmp[i]] = 0;
        }
        l = 0;

        maxDu = std::max(maxDu, Du);
        sumDu += Du;
    }
    
    uint32_t Dv = 0, maxDv, rdv = 0;
    double sumDv = 0.0;
    // while(t--) {
        // uint32_t v = chooseV(e);
    for(uint32_t v = 0; v < g->n2; v += rd) {
        rdv++;
        for(uint32_t i = g->pV[v]; i < g->pV[v + 1]; i++) {
            uint32_t u = g->e2[i];
    // if(u >= g->n1) {
    //     printf("n2 %u, %u, %u, i %u\n", g->n2, v, u, i);fflush(stdout);
    // }
            for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
                uint32_t w = g->e1[j];

                if(w > v) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(int i = 0; i < l; i++) {
            if(sum[tmp[i]] >= (uint32_t)p) {
                Dv++;
            }
            sum[tmp[i]] = 0;
        }
        l = 0;

        maxDv = std::max(maxDv, Dv);
        sumDv += Dv;
    }
    
    sumDu = sumDu / rdu * g->n1;
    sumDv = sumDv / rdv * g->n2;

    double avgDu = std::max(2.0, sumDu / g->n1);
    double avgDv = std::max(2.0, sumDv / g->n2);

    double costU = pow(avgDu, p - 2) * sumDu;
    double costV = pow(avgDv, q - 2) * sumDv;
    // printf("cost:%.2f %.2f, du %u, dv %u\n", costU, costV, Du, Dv);

    return costU > costV;
}

long double BCListPlusPlus::exactCount() {

    collect2HopNeighbors();
    // printf("collect 2\n"); fflush(stdout);

    S.resize(g->n2);

    uint32_t maxDegree = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        maxDegree = std::max(maxDegree, (uint32_t)H.lists[u].size());
    }


    tmpNodes.resize(p + 1);
    tmpNodes[0].resize(g->n1);
    for(int i = 1; i <= p; i++) {
        tmpNodes[i].resize(maxDegree);
    }
    H.nodes.resize(g->n1);

    H.d.resize(p + 1);
    for(int i = 0; i <= p; i++) {
        H.d[i].resize(g->n1);
    }
    for(uint32_t u = 0; u < g->n1; u++) {
        H.d[0][u] = H.lists[u].size();
    }

    ans = 0;

    layerBasedListing(0, g->n1, g->n2);

    printf("ans %.2Lf\n", ans);

    // unsigned long long int res = ans;

    return ans ; 
    // fflush(stdout);
}

void BCListPlusPlus::layerBasedListing(int l, int pH, int pS) {


    if(l == p) {
        ans += C[pS][q];
        return;
    }
    H.nodes.copy(tmpNodes[l].data(), pH);
// for(int i = 0; i < pH; i++) {
//     printf("%u ", tmpNodes[l][i]);
// }printf("\n");fflush(stdout);

auto t1 = std::chrono::steady_clock::now();
    for(int j = 0; j < pH; j++) {

        uint32_t u = tmpNodes[l][j];
        // uint32_t u = H.nodes[l][j];
        if(H.lists[u].size() < uint32_t(p - l - 1)) {
            continue;
        }

        int pSS = 0;
        if(l == 0) {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                S.changeTo(g->e1[i], pSS++);
            }
        }
        else {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                if(S.idx(g->e1[i]) < (uint32_t)pS) {
                    S.changeTo(g->e1[i], pSS++);
                }
            }
        }

        if(pSS < q) continue;
        
        int pHH = 0;

        for(int i = 0; i < H.d[l][u]; i++) {
            auto v = H.lists[u][i];
            if(H.nodes.idx(v) < (uint32_t)pH) {
                H.nodes.changeTo(v, pHH++);
            }
            // if(H.nodes[l].idx(v) < (uint32_t)pH) {
            //     H.nodes[l + 1].changeTo(v, pHH++);
            // }
        }

        if(l + 1 < p)
        for(int i = 0; i < pHH; i++) {
            // uint32_t u = H.nodes[l + 1][i];
            uint32_t u = H.nodes[i];
            int d = H.d[l][u];
            for(int k = 0; k < d; k++) {
                auto v = H.lists[u][k];
                if(H.nodes.idx(v) >= pHH) {
                    std::swap(H.lists[u][k], H.lists[u][--d]);
                    --k;
                }
                // if(H.nodes[l + 1].idx(v) >= pHH) {
                //     std::swap(H.lists[u][k], H.lists[u][--d]);
                //     --k;
                // }
            }
            H.d[l + 1][u] = d;
        }

        layerBasedListing(l + 1, pHH, pSS);
    }
}

void BCListPlusPlus::collect2HopNeighbors() {
    H.lists.resize(g->n1);

    std::vector<uint32_t> sum(g->n1);
    std::vector<uint32_t> tmp(g->n1);
    uint32_t l = 0;

double twotwo = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(uint32_t i = 0; i < l; i++) {
            uint32_t w = tmp[i];
            if(sum[w] >= q) {
twotwo += C[sum[w]][q];
                H.lists[u].push_back(w);
            }
            sum[w] = 0;
        }
        l = 0;
    }
// printf("2-2 clique %.0f\n", twotwo);
}

// Batch version of naive_biclique_with_vertex_sampling for P=3
std::vector<long double> naive_biclique_with_vertex_sampling_batch(BiGraph& g, unsigned long seed, int P___, double sampling_ratio, int num_samples, int round) {
    cout << "Running naive vertex sampling batch algorithm for P=3 (round " << round << ")..." << endl;
    
    BiGraph g2(g); 
    // Use the provided seed directly for noisy graph construction
    unsigned long noisy_graph_seed = seed;
    construct_noisy_graph(g, g2, noisy_graph_seed);

    cout << "Original noisy graph: n1=" << g2.num_v1 << ", n2=" << g2.num_v2 << ", edges=" << g2.num_edges << endl;

    // Initialize results for Q = [4, 5, 6, 7, 8, 9, 10]
    std::vector<long double> results(7, 0.0);

    // Apply vertex sampling if sampling_ratio < 1.0
    if (sampling_ratio < 1.0) {
        cout << "Applying vertex sampling with ratio: " << sampling_ratio << " (" << num_samples << " samples)" << endl;
        
        // For each sample, create one sampled subgraph and estimate all Q values on it
        for (int sample = 0; sample < num_samples; sample++) {
            // Create vertex-sampled graph (once per sample)
            BiGraph g2_sampled;
            g2_sampled.init(g2.num_v1, g2.num_v2);
            
            // Initialize random number generator with different seed for each sample
            init_genrand(seed + 12345 + sample * 1000);
            
            // Sample vertices from both partitions
            vector<bool> sampled_v1(g2.num_v1, false);
            vector<bool> sampled_v2(g2.num_v2, false);
            
            int sampled_v1_count = 0;
            int sampled_v2_count = 0;
            
            // Sample vertices from partition 1 (upper)
            for (int u = 0; u < g2.num_v1; u++) {
                if (genrand_real2() < sampling_ratio) {
                    sampled_v1[u] = true;
                    sampled_v1_count++;
                }
            }
            
            // Sample vertices from partition 2 (lower)
            for (int v = 0; v < g2.num_v2; v++) {
                if (genrand_real2() < sampling_ratio) {
                    sampled_v2[v] = true;
                    sampled_v2_count++;
                }
            }
            
            // Add edges between sampled vertices
            int sampled_edges = 0;
            for (int u = 0; u < g2.num_v1; u++) {
                if (sampled_v1[u]) {
                    for (auto v : g2.neighbor[u]) {
                        int v_lower = v - g2.num_v1;  // Convert to lower partition index
                        if (sampled_v2[v_lower]) {
                            g2_sampled.addEdge(u, v);
                            sampled_edges++;
                        }
                    }
                }
            }
            
            if (sample == 0) {  // Only print details for first sample
                cout << "Sample " << (sample+1) << " - Sampled vertices: " << sampled_v1_count << " (upper), " << sampled_v2_count << " (lower)" << endl;
                cout << "Sample " << (sample+1) << " - Sampled graph edges: " << sampled_edges << endl;
            }
            
            // Convert sampled graph once per sample
            biGraph convertedGraph = convertBiGraphTobiGraph(g2_sampled);
            
            // For this sample, estimate all Q values on the same sampled subgraph
            for (int q_idx = 0; q_idx < 7; q_idx++) {
                int Q = q_idx + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
                
                // Count bicliques on the same sampled graph for this Q
                BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, P___, Q);
                long double sparse_count = counter->exactCount();
                delete counter;
                
                // Unbiased estimate: scale by 1/sampling_ratio^(p+q) (for vertex sampling)
                long double scaling_factor = pow(sampling_ratio, P___ + Q);
                long double sample_estimate = sparse_count / scaling_factor;
                
                // Accumulate this sample's estimate
                results[q_idx] += sample_estimate;
            }
        }
        
        // Average the estimates across all samples
        for (int q_idx = 0; q_idx < 7; q_idx++) {
            results[q_idx] = results[q_idx] / num_samples;
        }
        
    } else {
        // No sampling - use original method for all Q values
        biGraph convertedGraph = convertBiGraphTobiGraph(g2);
        cout << "Converted graph: n1=" << convertedGraph.n1 << ", n2=" << convertedGraph.n2 << ", m=" << convertedGraph.m << std::endl;
        
        for (int q_idx = 0; q_idx < 7; q_idx++) {
            int Q = q_idx + 4;  // Q = 4, 5, 6, 7, 8, 9, 10
            cout << "Processing Q = " << Q << " (exact counting)" << endl;
            
            BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, P___, Q);
            results[q_idx] = counter->exactCount();
            delete counter;
        }
    }
    
    cout << "Batch results for Q=[4,5,6,7,8,9,10]: ";
    for(int q_idx = 0; q_idx < 7; q_idx++) {
        cout << "Q" << (q_idx + 4) << "=" << results[q_idx] << " ";
    }
    cout << endl;
    
    return results;
}
