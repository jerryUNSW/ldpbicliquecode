#include "biclique.h"
using namespace std;
using namespace std::chrono;

extern long double _cate, _wedge, _btf;
vector<long double> estis, naive_estis;
extern bool one_round, count_cate; 

long double p, Eps, Eps0, Eps1, Eps2;
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

unsigned long long int  real ;

int P___, K___ ; 
// biclique related 
vector<vector<int>> up_options, lo_options; 

long double avg_estimated_variance = 0 ;

extern bool two_noisy_graph_switch; 

// I  think we can have another way of counting caterpillars. 
// for each vertex, we consider the wedges. 
int main(int argc, char *argv[]) {
    // input parser:
    long double param = stold(argv[1]);
    string dataset = argv[2];
    num_rounds = atoi(argv[3]);
    int algorithm_switch = atoi(argv[4]);

    P___ = atoi(argv[5]);
    K___ = atoi(argv[6]);

    cout<<"P___ = "<<P___ <<endl;
    cout<<"K___ = "<<K___ <<endl;

    // initialize time
    RR_time = 0, server_side_time = 0, naive_server_side = 0;

    std::mt19937 rng(std::random_device{}());  // for seeding

    BiGraph g(dataset);

    long double fill_rate = g.num_edges * 1.0 / ((double)g.num_v1 * (double)g.num_v2);

    cout << "fill_rate = " << fill_rate << endl;

    Eps = param;
    cout << "Eps = " << Eps << endl;

    vector<long double> m__;

    unsigned long long int btf = BFC_EVP(g);
    cout<<"BTF =  "<<btf<<endl;

    // biGraph convertedGraph = convertBiGraphTobiGraph(g);

    // std::cout << "Converted graph: n1=" << convertedGraph.n1 
    //           << ", n2=" << convertedGraph.n2 
    //           << ", m=" << convertedGraph.m << std::endl;
    // // Create and use BCListPlusPlus
    // BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, P___, K___);
    // cout<<"cliq count = "<<counter->exactCount()<<endl;
        
    // grab exact biclique coutns from the sqlite database
    check_exact_result_in_DB(P___, K___, dataset, g); 

    estis.resize(num_rounds);
    naive_estis.resize(num_rounds);
    vector<long double> rel_err;

    double t0 = omp_get_wtime();

    for (iteration = 0; iteration < num_rounds; iteration++) {

        // we need to implement the naive algorithm for pq bcilqiue counting.
        // option 1: RR -> noisy graph G' -> return 
        // option 2: 

        if (algorithm_switch == 0) {
            cout << "Naive algorithm for biclique counting" << endl;
            cout << "EPS = " << Eps << endl;     

            p = 1.0 / (exp(Eps) + 1.0);
            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;


            // printMemoryUsage();
            estis[iteration] = naive_biclique(g, seed, P___, K___);
            // printMemoryUsage();


            cout << "estimate = " << estis[iteration] << endl;

            long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

            cout << "relative error = " << relative_error << endl;
            rel_err.push_back(relative_error);
            cout << endl;
        }
        if (algorithm_switch == 1) {

            if(P___ == 2 && K___==2){     
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
                    // because this is our contribution to the ADV algorithm
                    long double esti1 = one_round_btf(g, seed);
                    unsigned int seed2 = rng();
                    long double esti2 = one_round_btf(g, seed2);
                    estis[iteration] = (esti1 + esti2)/2;
                }else{
                    estis[iteration] = one_round_btf(g, seed);
                }
                long double rel = abs(estis[iteration] - real) * 1.0 / real;
                cout << "estimate = " << estis[iteration] << endl;
                cout << "relative error = " << rel << endl;
                rel_err.push_back(rel);
                cout << endl;
            }
            else{
                // one-round biclique algorithm for general P, Q values
                cout << "Oneround, handling genral P, and Q values" << endl;
                cout << "EPS = " << Eps << endl;
                unsigned int seed = rng();
                cout << "random seed = " << seed << endl;

                // estis[iteration] = wedge_based_btf_avg(g, seed);
                // estis[iteration] = VP_wedge_based_two_round_btf(g, seed);

                // real = 1775 ; 
                // real = 5400;

                p = 1.0 / (exp(Eps) + 1.0);

                estis[iteration] = one_round_biclique(g, 
                    seed, P___, K___);
                // cout << "estimate = " << estis[iteration] << endl;
                
                std::cout << "estimate = " << std::fixed << std::setprecision(10) << estis[iteration] << std::endl;

                cout<<"real = "<< real<<endl;

                long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

                cout << "relative error = " << relative_error << endl;
                rel_err.push_back(relative_error);
                cout << endl;
            }
        } 
        /*
        if (algorithm_switch == 2) {
            cout << "\nOld Multi-round algorithm for BTF" << endl;
            cout << "EPS = " << Eps << endl;

            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;

            estis[iteration] = two_round_btf(g, seed);

            cout << "estimate = " << estis[iteration] << endl;

            long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

            cout << "relative error = " << relative_error << endl;
            rel_err.push_back(relative_error);
            cout << endl;
        }
        */
        if (algorithm_switch >= 3) {
            
            // the common neighbor estimation inspired algorithm 
            cout << "\nWedge-based biclique counting algorithm" << endl;

            cout << "EPS = " << Eps << endl;

            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;


            if(P___ == 2){
                // this has the average version.
                if(algorithm_switch==4){
                    two_noisy_graph_switch = true ;
                }
                else{
                    two_noisy_graph_switch = false ;
                }
                estis[iteration] = wedge_based_two_round_2_K_biclique(g, seed);
                // estis[iteration] = wedge_based_btf_avg(g, seed); // this is the avg version.
            }
            // need to implement two_noisy_graph_switch optimization for P == 3 
            else if(P___ == 3){
                estis[iteration] = wedge_based_two_round_3_K_biclique(g, seed);
            }
            // need to implement two_noisy_graph_switch optimization for P in general
            else{
                cout<<"P = "<<P___ <<endl;
                cout<<"Q = "<<K___ <<endl;

                estis[iteration] = wedge_based_two_round_general_biclique(g, seed, P___, K___);
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
    
    printf("# Mean = %Lf\n", calculateMean(estis));
    cout << "real count = " << real << endl;

    cout << "adv rel err = " << calculateMean(rel_err) << endl;

    // cout<<"esti_var_f mean: = "<<calculateMean(naive_estis) << endl;
    // cout<<"observed var(f_test): = "<<computeVariance(naive_estis) << endl;

    // this is the variance of f1
    // cout<<"sample variance = "<<computeVariance(estis)<<endl;

    // cout<<"esti varinace =  "<<avg_estimated_variance / estis.size()<<endl;


    // vector<long double> naive_err;
    // for (auto xxx : naive_estis) {
    //     naive_err.push_back(abs(xxx - real) * 1.0 / real);
    // }
    // if (algorithm_switch == 4) {
    //     cout << "naive rel err = " << calculateMean(naive_err) << endl;
    // }




    // efficiency evaluations:
    /*
    if (algorithm_switch == 1) {
        printf("RR time = %Lf\n", RR_time/num_rounds);
        printf("adv server time = %Lf\n", server_side_time);
        printf("naive server time = %Lf\n", naive_server_side);
        printf("rr cost = %Lf\n", communication_cost);
    } else {
        printf("RR time = %Lf\n", RR_time/num_rounds);
        printf("deg time = %Lf\n", deg_esti_time/num_rounds);
        printf("local time = %Lf\n", local_count_time/num_rounds);
        printf("com cost = %Lf\n", communication_cost);
    }
    */
    return 0;
}