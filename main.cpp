#include "ldp-btf.h"
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
extern stats::rand_engine_t engine;  // Declare the engine as extern
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0, real_stars = 0,
			RR_time, server_side_time, naive_server_side,
			real = 0, local_count_time = 0, deg_esti_time = 0,
			communication_cost = 0;

int main(int argc, char *argv[]) {
    // input parser:
    long double param = stold(argv[1]);
    string dataset = argv[2];
    num_rounds = atoi(argv[3]);
    int algorithm_switch = atoi(argv[4]);

    if (algorithm_switch == 1) {
        one_round = true;
    } else {
        one_round = false;
    }
    // check the type of motif
    count_cate = false;

    // initialize time
    RR_time = 0, server_side_time = 0, naive_server_side = 0;

    std::mt19937 rng(std::random_device{}());  // for seeding

    BiGraph g(dataset);

    long double fill_rate = g.num_edges * 1.0 / (g.num_v1 * g.num_v2);

    cout << "fill_rate = " << fill_rate << endl;

    if (!input_epsilon) {
        Eps = compute_epsilon(param, g.num_edges, g.num_v1, g.num_v2);
        cout << "Eps = " << Eps << endl;
    } else {
        Eps = param;
        cout << "Eps = " << Eps << endl;
    }

    double t0 = omp_get_wtime();

    vector<long double> m__;

    _btf = BFC_EVP(g);
    // _cate = get_cate(g);
    cout << "btf = " << _btf << endl;
    real = _btf;

    estis.resize(num_rounds);
    naive_estis.resize(num_rounds);
    vector<long double> rel_err;

    for (iteration = 0; iteration < num_rounds; iteration++) {
        if (one_round) {
            cout << "\nOne round algorithm" << endl;
            cout << "EPS = " << Eps << endl;
            // flip probability
            p = 1.0 / (exp(Eps) + 1.0);
            unsigned int seed = rng();
            cout << "seed = " << seed << endl;

            estis[iteration] = one_round_btf(g, seed);

            long double rel = abs(estis[iteration] - real) * 1.0 / real;
            cout << "estimate = " << estis[iteration] << endl;
            cout << "relative error = " << rel << endl;
            rel_err.push_back(rel);
            cout << endl;
        } else {
            cout << "\nMultiple round algorithms" << endl;
            cout << "EPS = " << Eps << endl;

            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;

            estis[iteration] = two_round_btf(g, seed);

            cout << "estimate = " << estis[iteration] << endl;

            double relative_error =
                count_cate ? abs(estis[iteration] - _cate) * 1.0 / _cate
                           : abs(estis[iteration] - _btf) * 1.0 / _btf;
            cout << "relative error = " << relative_error << endl;
            rel_err.push_back(relative_error);
            cout << endl;
        }
    }
    printf("# Mean = %Lf\n", calculateMean(estis));
    cout << "real count = " << real << endl;
    cout << "adv rel err = " << calculateMean(rel_err) << endl;
    vector<long double> naive_err;
    for (auto xxx : naive_estis) {
        naive_err.push_back(abs(xxx - real) * 1.0 / real);
    }
    if (one_round) {
        cout << "naive rel err = " << calculateMean(naive_err) << endl;
    }

    // efficiency evaluations:

    if (one_round) {
        // printf("RR time = %Lf\n", RR_time/num_rounds);
        // printf("adv server time = %Lf\n", server_side_time);
        // printf("naive server time = %Lf\n", naive_server_side);
        printf("rr cost = %Lf\n", communication_cost);
    } else {
        // printf("RR time = %Lf\n", RR_time/num_rounds);
        // printf("deg time = %Lf\n", deg_esti_time/num_rounds);
        // printf("local time = %Lf\n", local_count_time/num_rounds);
        printf("com cost = %Lf\n", communication_cost);
    }
    double t1 = omp_get_wtime();
    double seconds = t1 - t0;
    printf("time:%f\n", seconds);
    return 0;
}