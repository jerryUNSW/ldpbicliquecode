#include "biclique.h"
using namespace std;
using namespace std::chrono;

extern long double _cate, _wedge, _btf;
vector<long double> estis, naive_estis;
extern bool one_round, count_cate; 

long double p, Eps, Eps0, Eps1, Eps2;
// privacy budget allocation ratios
long double alpha0, alpha1, alpha2;
// private degrees:
vector<int> priv_deg;
int priv_dmax_1, priv_dmax_2, iteration, num_rounds; 
bool input_epsilon =true; 
extern int alpha;
extern bool count_cc; 
extern stats::rand_engine_t engine;  // Declare the engine as extern
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0, real_stars = 0,
			RR_time, server_side_time, naive_server_side, local_count_time = 0, deg_esti_time = 0,
			communication_cost = 0;


extern bool eva_comm ;


long double real = 0.0;  // Use long double consistently for ground truth
long double real_ld = 0.0;  // High precision ground truth for large numbers (kept for compatibility)

int P___, K___ ; 
// biclique related 
vector<vector<int>> up_options, lo_options; 

long double avg_estimated_variance = 0 ;

extern bool two_noisy_graph_switch; 

extern bool multi_estimator_switch ; 

bool use_probability_filtering = false;

// Noisy edge sampling ratio for naive algorithm
double noisy_edge_sampling_ratio = 0.1;  // Default: 10% edge sampling for speed

// 12th parameter: T (number of p-tuples to sample for general P)
long long T = 100000LL;  // Default: 100,000 p-tuples

// I  think we can have another way of counting caterpillars. 
// for each vertex, we consider the wedges. 
int main(int argc, char *argv[]) {

    eva_comm = false;

    if(eva_comm){
        cout<<"evaluating communication costs"<<endl;
    }

    // input parser:
    long double param = stold(argv[1]);
    string dataset = argv[2];
    num_rounds = atoi(argv[3]);
    int algorithm_switch = atoi(argv[4]);

    P___ = atoi(argv[5]);
    K___ = atoi(argv[6]);

    // Privacy budget allocation ratios (with default values)
    if (argc >= 10) {
        alpha0 = stold(argv[7]);
        alpha1 = stold(argv[8]);
        alpha2 = stold(argv[9]);
        
        // Validate that they sum to 1.0
        long double sum = alpha0 + alpha1 + alpha2;
        if (abs(sum - 1.0) > 1e-6) {
            cerr << "Error: alpha0 + alpha1 + alpha2 must equal 1.0, got " << sum << endl;
            return 1;
        }
    } else {
        // Default allocation: (0.05, 0.6, 0.35)
        alpha0 = 0.05;
        alpha1 = 0.6;
        alpha2 = 0.35;
    }

    // 11th parameter: probability filtering (0=disabled, 1=enabled)
    // if (argc >= 11) {
    //     use_probability_filtering = (atoi(argv[10]) == 1);
    // } else {
    //     use_probability_filtering = false; // Default to disabled
    // }
    use_probability_filtering = false; // Default to disabled

    // 11th parameter: noisy edge sampling ratio (0.0-1.0, default 0.1 = 10% sampling for speed)
    if (argc >= 11) {
        noisy_edge_sampling_ratio = stod(argv[10]);
        if (noisy_edge_sampling_ratio <= 0.0 || noisy_edge_sampling_ratio > 1.0) {
            cerr << "Error: noisy edge sampling ratio must be in (0.0, 1.0], got " << noisy_edge_sampling_ratio << endl;
            return 1;
        }
    } else {
        noisy_edge_sampling_ratio = 1.0; // Default: no sampling
    }

    // 12th parameter: T (number of p-tuples to sample for general P)
    if (argc >= 12) {
        T = stoll(argv[11]);
        if (T <= 0) {
            cerr << "Error: T must be positive, got " << T << endl;
            return 1;
        }
    } else {
        T = 100000LL; // Default: 100,000 p-tuples
    }


    cout<<"P___ = "<<P___ <<endl;
    cout<<"K___ = "<<K___ <<endl;
    cout<<"Budget allocation: alpha0="<<alpha0<<", alpha1="<<alpha1<<", alpha2="<<alpha2<<endl;
    cout<<"Probability filtering: "<<(use_probability_filtering ? "ENABLED" : "DISABLED")<<endl;
    cout<<"Noisy edge sampling ratio: "<<noisy_edge_sampling_ratio<<endl;
    cout<<"T (p-tuples to sample): "<<T<<endl;

    // initialize time
    RR_time = 0, server_side_time = 0, naive_server_side = 0;

    std::mt19937 rng(std::random_device{}());  // for seeding

    BiGraph g(dataset);

    long double fill_rate = g.num_edges * 1.0 / ((double)g.num_v1 * (double)g.num_v2);

    cout << "fill_rate = " << fill_rate << endl;

    Eps = param;
    cout << "Eps = " << Eps << endl;

    vector<long double> m__;

    // unsigned long long int btf = BFC_EVP(g);
    // cout<<"BTF =  "<<btf<<endl;

    // grab exact biclique coutns from the sqlite database
    fetch_or_compute_biclique_count(P___, K___, dataset, g); 

    estis.resize(num_rounds);
    naive_estis.resize(num_rounds);
    vector<long double> rel_err;

    double t0 = omp_get_wtime();

    for (iteration = 0; iteration < num_rounds; iteration++) {

        // the naive algorithm for pq bcilqiue counting.
        if (algorithm_switch == 0) {
            cout << "Naive algorithm for biclique counting" << endl;
            cout << "EPS = " << Eps << endl;
            if (noisy_edge_sampling_ratio < 1.0) {
                cout << "Using VERTEX SAMPLING approach: " << (noisy_edge_sampling_ratio * 100) << "% vertex sampling with 10 samples per iteration" << endl;
            } else {
                cout << "Using EXACT COUNTING approach: no sampling" << endl;
            }     

            p = 1.0 / (exp(Eps) + 1.0);
            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;


            // printMemoryUsage();
            if (noisy_edge_sampling_ratio < 1.0) {
                cout<<"naive_biclique_with_vertex_sampling"<<endl;
                estis[iteration] = naive_biclique_with_vertex_sampling(g, seed, P___, K___, noisy_edge_sampling_ratio, 10);
            } else {
                cout<<"naive_biclique"<<endl;
                estis[iteration] = naive_biclique(g, seed, P___, K___);
            }
            // printMemoryUsage();

            cout << "estimate = " << estis[iteration] << endl;

            long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

            cout << "relative error = " << relative_error << endl;
            rel_err.push_back(relative_error);
            cout << endl;
        }
        if (algorithm_switch == 1) {
            // one round algorithms (same communication costs as naive)
            if(P___ == 2 && K___== 2){     
                cout << "Oneround-BTF (existing DBE algorithm)" << endl;
                cout << "EPS = " << Eps << endl;
                // flip probability
                p = 1.0 / (exp(Eps) + 1.0);
                unsigned int seed = rng();
                cout << "seed = " << seed << endl;

                bool avg_btf_switch = false;
                if(avg_btf_switch){
                    // why can we do this? 
                    // in theory, we are able to build two noisy graphs. 
                    // based on each, we can have a BTF estimate. 
                    // however, in experiment evaluations, we dont consider this case 
                    // because the double noisy graph technique
                    // is our contribution to the ADV algorithm
                    long double esti1 = one_round_btf(g, seed);
                    unsigned int seed2 = rng();
                    long double esti2 = one_round_btf(g, seed2);
                    estis[iteration] = (esti1 + esti2)/2;
                }else{
                    // this is a specific algorithm, so it is faster
                    estis[iteration] = one_round_btf(g, seed);

                    // let's see if ths is faster 
                    // estis[iteration] = one_round_biclique_2_K(g, K___, seed); 
                }
                long double rel = abs(estis[iteration] - real) * 1.0 / real;
                cout << "estimate = " << estis[iteration] << endl;
                cout << "relative error = " << rel << endl;
                rel_err.push_back(rel);
                cout << endl;
            }
            if(P___ == 2 && K___ > 2){
                // one-round biclique algorithm for general P, Q values
                cout << "Oneround" << endl;
                cout<<"P___ == " << P___ <<endl;
                cout<<"K___ == " << K___ <<endl;
                cout << "epsilon = " << Eps << endl;
                unsigned int seed = rng();
                cout << "random seed = " << seed << endl;

                p = 1.0 / (exp(Eps) + 1.0);

                // estis[iteration] = one_round_biclique_2_3(g, seed);


                estis[iteration] = one_round_biclique_2_K(g, K___, seed); 

                // estis[iteration] = one_round_biclique(g, seed, P___, K___);

                std::cout << "estimate = " << std::fixed << std::setprecision(10) << estis[iteration] << std::endl;

                cout<<"real = "<< real<<endl;

                long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

                cout << "relative error = " << relative_error << endl;
                rel_err.push_back(relative_error);
                cout << endl;

            }

            if(P___ > 2 ){
                // one-round biclique algorithm for general P, Q values
                cout << "Oneround, P___ > 2" << endl;
                cout << "EPS = " << Eps << endl;
                unsigned int seed = rng();
                cout << "random seed = " << seed << endl;

                // estis[iteration] = wedge_based_btf_avg(g, seed);
                // estis[iteration] = VP_wedge_based_two_round_btf(g, seed);

                p = 1.0 / (exp(Eps) + 1.0);

                estis[iteration] = one_round_biclique(g, seed, P___, K___);
                // cout << "estimate = " << estis[iteration] << endl;
                
                std::cout << "estimate = " << std::fixed << std::setprecision(10) << estis[iteration] << std::endl;

                cout<<"real = "<< real<<endl;

                long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

                cout << "relative error = " << relative_error << endl;
                rel_err.push_back(relative_error);
                cout << endl;
            }
        } 
        if (algorithm_switch >= 2) {

            // manipulate swith
            if(algorithm_switch == 2){
                cout << "\nADV" << endl;
                multi_estimator_switch = false; 
                two_noisy_graph_switch = false ;
            }
            if(algorithm_switch == 3){
                cout << "\nADV+" << endl;
                multi_estimator_switch = true; 
                two_noisy_graph_switch = false ;
            }
            if(algorithm_switch == 4){
                cout << "\nADV++" << endl;
                multi_estimator_switch = true; 
                two_noisy_graph_switch = true ;
            }
            
            // in this algorithm, we only use one vertex as the source of estimation.
            cout << "EPS = " << Eps << endl;
            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;
            if(P___ == 2){
                // this will be the basis of the working example
                estis[iteration] = wedge_based_two_round_2_K_biclique(g, seed);

                // layer_based 
                // estis[iteration] = layer_based_wedge_based_two_round_2_K_biclique(g, seed);
            }
            else if(P___ == 3){
                // this is ready.
                // estis[iteration] = wedge_based_two_round_3_K_biclique_rejection_sampling(g, seed);
                estis[iteration] = wedge_based_two_round_3_K_biclique(g, seed);
            }
            // need to implement two_noisy_graph_switch optimization for P in general
            else{
                cout<<"P = "<<P___ <<endl;
                cout<<"Q = "<<K___ <<endl;
                // we do not test the communication of this
                estis[iteration] = wedge_based_two_round_general_biclique(g, seed, P___, K___, T);            
                // cout<<"Need to implement baseline for general P values"<<endl;
            }
            
            // cout << "estimate = " << estis[iteration] << endl;
            std::cout << "estimate = " << std::fixed << std::setprecision(10) << estis[iteration] << std::endl;
            cout<<"real = "<<real <<endl;
            long double relative_error = abs(estis[iteration] - real) * 1.0 / real;
            cout << "relative error = " << relative_error << endl;
            rel_err.push_back(relative_error);
            cout << endl;
        }
    }

    double t1 = omp_get_wtime();
    double seconds = t1 - t0;
    printf("time:%f\n", seconds);


    if(eva_comm){
        // byte --> metabyte 
        long double communication_cost_MB = communication_cost / (1024.0 * 1024.0);
        printf("com cost (MB) = %Lf\n", communication_cost_MB);
        return 0;
    }
    
    printf("# Mean = %Lf\n", calculateMean(estis));
    cout << "real count = " << real << endl;

    cout << "adv rel err = " << calculateMean(rel_err) << endl;
    return 0;
}