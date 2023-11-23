#include "ldp-btf.h"
#include "mt19937ar.h"

using namespace std;

long double _cate, _wedge, _btf;
vector<int> upper_sample, lower_sample;
double sample_ratio = 1.0; 
unordered_map<int, bool> in_sampled_up, in_sampled_lo;  
bool samling_one_round = false;
long double verified = 0 , not_verified = 0;

extern long double Eps, Eps0, Eps1, Eps2, p, m3__, m2__, m1__, m0__;
// priv deg related. 
extern vector<int> priv_deg; 
extern vector<long double> naive_estis;
extern int priv_dmax_1, priv_dmax_2; 
extern vector<vector<int>> up_options, lo_options; 

extern long double communication_cost; 

extern int iteration;

stats::rand_engine_t engine(std::time(0)); // used to be 1776

bool one_round = false, edge_clipping = true;

bool count_cate = false;

bool sampling_noisy_graph = false ; 

double p____ = 0.5; // this is the sampling ratio.

int alpha = 10; // maybe this will make two-round worse on BX.

bool eva_comm = false;

extern long double RR_time, server_side_time, naive_server_side;

extern long double local_count_time, deg_esti_time; 

void private_estimate_of_degrees(BiGraph& g){
	// Eps0 = 0.05;

	// private estimate degrees. 
	priv_dmax_1 = 0 ;       
	priv_dmax_2 = 0 ; 
	priv_deg.resize(g.num_nodes());

	for(int i=0;i<g.num_nodes();i++){

		priv_deg[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 

		if(edge_clipping){
			priv_deg[i]+=alpha;
		}
		// long double tmp = add_geometric_noise(g.degree[i], 1, Eps0);  
		// we found that geometric 
		// cout<<"deg = "<<g.degree[i]<<",\t"<<priv_deg[i]<<",\t"; cout<<tmp<<endl; 
		if(g.is_upper(i)){
			priv_dmax_1 = priv_dmax_1 > priv_deg[i] ? priv_dmax_1 : priv_deg[i] ; 
		}else{
			priv_dmax_2 = priv_dmax_2 > priv_deg[i] ? priv_dmax_2 : priv_deg[i] ; 
		}
		// if(choose_upper && g.is_upper(i))x1+=priv_deg[i];
		// if(!choose_upper && g.is_lower(i))x1+=priv_deg[i];
	}

	// cout<<"esitmated  = "<<x1<<endl;
	// g.num_edges = x1;
}

double my_genrand_real2(){return genrand_real2(); }

void my_init_genrand(unsigned long seed){init_genrand(seed);}

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed){

    const int range_from  = g2.num_v1;
    const int range_to    = g2.num_nodes()-1;
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(range_from, range_to);

	init_genrand(seed);
	int flip1 = 0; 
	int visited_vertices = 0;
	int total_vertices = g2.num_v1;
	int ten_percent = total_vertices / 5;
	long double max_time_per_user = -1;

	for(int i=0;i<g2.num_v1;i++){
		visited_vertices++; 
		// if (visited_vertices % ten_percent == 0) {
		// 	int progress = visited_vertices * 100 / total_vertices;
		// 	cout << "Processed " << progress << "% of vertices" << endl;
		// }
		double tx = omp_get_wtime();
		for(int j=g2.num_v1;j<g2.num_nodes();j++){
			if(std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) != g.neighbor[i].end() ){
				if(genrand_real2() >= p ){ // 1  --> 1 
					g2.addEdge(i,j);
					flip1++;
				}	
			}else{
				if(genrand_real2() < p){   // 0 --> 1
					g2.addEdge(i,j);
					flip1++;	
				}
			}
		}
		double ty = omp_get_wtime();
		max_time_per_user = max_time_per_user > (ty-tx) ? max_time_per_user : (ty-tx);
	}
	/*
	#pragma omp parallel
    {
		BiGraph g2_thread_copy(g);

        #pragma omp for schedule(dynamic) // You can choose a different scheduling strategy
		for(int i=0;i<g2_thread_copy.num_v1;i++){
			// double tx = omp_get_wtime();
			for(int j=g2_thread_copy.num_v1;j<g2_thread_copy.num_nodes();j++){

				g2_thread_copy.addEdge(i,j);

				// this checking requires involvement with g. 
				// if(std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) != g.neighbor[i].end() ){
				// 	if(genrand_real2() >= p ){ // 1  --> 1 
				// 		g2_thread_copy.addEdge(i,j);
				// 	}	
				// }else{
				// 	if(genrand_real2() < p){   // 0 --> 1
				// 		g2_thread_copy.addEdge(i,j);
				// 	}
				// }
			}
			// double ty = omp_get_wtime();
			// max_time_per_user = max_time_per_user > (ty-tx) ? max_time_per_user : (ty-tx);
		}
		#pragma omp critical
		{
			for(int i=0;i<g.num_v1;i++){
				for(auto j:g2_thread_copy.neighbor[i]){
					g2.addEdge(i,j);
				}
			}
		}
    }
	flip1 = g2.num_edges ; 
	*/

	// RR_time += max_time_per_user; 
	// the dominating cost is incurred in server side butterfly counting on the dense noisy graph

	cout<<"noisy edges = "<<flip1<<endl;

	communication_cost+= flip1*sizeof(int); 

	long double expected_E = g.num_edges*(1-p) +(g.num_v1*g.num_v2-g.num_edges)*p; 

	// long double expected_E = (g.num_edges * (1 - p)) + ((static_cast<long double>(g.num_v1) * static_cast<long double>(g.num_v2) - g.num_edges) * p);

	cout<<"expected E = "<< expected_E<<endl;

	g2.computePriority();
}

long double two_round_btf(BiGraph& g, unsigned long seed){
	
	// Phase 0. deg_esti_time records the maximum degree perturbation time.
	double t0 = omp_get_wtime();
	cout<<"private_estimate_of_degrees(g); "<<endl;
	Eps0 = Eps*0.1; 
	private_estimate_of_degrees(g); 

	// upload noisy degrees
	if(eva_comm) communication_cost+=g.num_nodes()*sizeof(int); 

	// Phase 1. RR 
	double t1 = omp_get_wtime();
	cout<<"construct_noisy_graph(g); "<<endl;
	Eps1 = Eps*0.5;
	p = 1.0 / (exp(Eps1) + 1.0);

	BiGraph g2(g);
	if(sampling_noisy_graph){

		cout<<"sampling ratio = "<<p____ <<endl;
		const int range_from  = g2.num_v1;
		const int range_to    = g2.num_nodes()-1;
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<int>  distr(range_from, range_to);
		init_genrand(seed);
		int flip1 = 0; 
		int visited_vertices = 0;
		int total_vertices = g2.num_v1;
		int ten_percent = total_vertices / 5;
		// long double max_time_per_user = -1;
		for(int i=0;i<g2.num_v1;i++){
			// visited_vertices++; 
			// if (visited_vertices % ten_percent == 0) {
			// 	int progress = visited_vertices * 100 / total_vertices;
			// 	cout << "Processed " << progress << "% of vertices" << endl;
			// }
			// double tx = omp_get_wtime();
			for(int j=g2.num_v1;j<g2.num_nodes();j++){
				if(std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) != g.neighbor[i].end() ){
					if(genrand_real2() >= p ){ // 1  --> 1 
						if(genrand_real2() >= p____) continue; // keep with probability p____
						g2.addEdge(i,j);
						flip1++;
					}	
				}else{
					if(genrand_real2() < p){   // 0 --> 1
						if(genrand_real2() >= p____) continue; // keep with probability p____
						g2.addEdge(i,j);
						flip1++;	
					}
				}
			}
			// double ty = omp_get_wtime();
			// max_time_per_user = max_time_per_user > (ty-tx) ? max_time_per_user : (ty-tx);
		}
		cout<<"noisy edges = "<<flip1<<endl;
		long double expected_E = g.num_edges*(1-p) +(g.num_v1*g.num_v2-g.num_edges)*p; 
		cout<<"expected E = "<< expected_E<<endl;
		g2.computePriority(); // is this necessary? 
	}else{
		construct_noisy_graph(g, g2, seed); // upload noisy edges
	}
	

	// Phase 2. local counting records the counting time 
	double t2 = omp_get_wtime();
	cout<<"local counting"<<endl;
	if(eva_comm){
		// for each vertex, it needs to download all vertex degrees. 
		communication_cost += g.num_nodes()* g.num_nodes() * sizeof(int); 

		// for each vertex, it needs to the whol noisy graph
		communication_cost += g.num_nodes() * g2.num_edges * sizeof(int); 

		// for each vertex, it needs to upload local count wi
		communication_cost += g.num_nodes() * sizeof(long double);// upload local counts
		return 0;
	}
	Eps2 = Eps*0.4;

	// use eps2
	long double global_sensitivity, sum = 0 ;

	// the following parallel implementation is trouble some
	// #pragma omp parallel for reduction(+:sum) // Use reduction to ensure sum is updated correctly
	for(int u=0;u<g.num_nodes();u++){

		// if(g.degree[u]==0) continue;

		if(edge_clipping && priv_deg[u]<=0) continue; 

		// when edge clipping is inplace, we can only visit at most priv_deg[u] neighbors for u
		long double s1 = 0, s2 =0, s3 = 0;
		
		unordered_map<vid_t,int> count_wedge(0);

		long double du = g.degree[u]; 

		if(edge_clipping  && (du>priv_deg[u]) ){
			du = priv_deg[u];
		}

		s3 += (du *(du-1)/2) * ( (g.is_upper(u) ? g.num_v1:g.num_v2)-1);  
		
		long double sum_deg_v = 0;
		int visited_nb = 0;
		for(auto v: g.neighbor[u]){
			if(edge_clipping && visited_nb==priv_deg[u]){
				break;
			}
			sum_deg_v += g2.degree[v]; 
			for(auto w: g2.neighbor[v]){
				if(u!=w){
					count_wedge[w]++;
				}
			}
			if(edge_clipping)visited_nb++;
		}
		
		for(auto ele:count_wedge){	
			// only execute this for butterfly counting
			if( (!count_cate) && (ele.second>=2) ){
				s1+=(ele.second-1)*ele.second/2;
			}
			s2 += ele.second * (du-1);
		}

		int deg_up = edge_clipping ? priv_deg[u] : g.degree[u]; 
		if(count_cate){
			// for caterpillar counting: we only need to count s2, and s3.
			if(g.is_upper(u)){
				global_sensitivity = 3*deg_up* g2.v2_max_degree + sum_deg_v + 2*p*deg_up*(g.num_v1-1); 
			}else{
				global_sensitivity = 3*deg_up* g2.v1_max_degree + sum_deg_v + 2*p*deg_up*(g.num_v2-1); 
			}
			sum += s2 - 2*p*s3;
		}else{
			// for buttrfly counting: 
			if(g.is_upper(u)){
				global_sensitivity = (1-p)*deg_up* g2.v2_max_degree + p*sum_deg_v + p*p*deg_up*(g.num_v1-1); 
			}else{
				global_sensitivity = (1-p)*deg_up* g2.v1_max_degree + p*sum_deg_v + p*p*deg_up*(g.num_v2-1); 
			}
			sum += s1 - p*s2 + p*p*s3; 
		}
		sum += stats::rlaplace(0.0, (global_sensitivity/Eps2), engine); // add calibrated noise
		// communication_cost += sizeof(long double);
	}
	// return sum/(4*(1-2*p)*(1-2*p));
	double t3 = omp_get_wtime();

	RR_time += t2 - t1; 
	deg_esti_time += t1 - t0; 
	local_count_time += t3-t2;

	if(count_cate){
		if(sampling_noisy_graph) sum/= p____ ; 
		return sum/(2*(1-2*p));
	}else{
		if(sampling_noisy_graph) sum/= p____*p____ ; 
		return sum/(4*(1-2*p)*(1-2*p));
	}
}

void compute_m3_m2(long double &m4, long double &m3, long double &m2, long double &m1, long double &m0, BiGraph& g2){

	long double caterpillar = 0;
	long double chopsticks = 0; 
	long double wedges=0; 

	long double g2_num_edges = g2.num_edges; 
	
	long double g2_num_v1 = g2.num_v1; 
	long double g2_num_v2 = g2.num_v2; 

	#pragma omp parallel for reduction(+:caterpillar, chopsticks)
	for(int i=0;i<g2.num_v1;i++){
		for(auto j:g2.neighbor[i]){
			// for each edge (i,j)
			caterpillar += (g2.degree[i]-1)*(g2.degree[j]-1);
			chopsticks += g2_num_edges-g2.degree[i]-g2.degree[j]+1;
		}
	}

	chopsticks/=2; // each is counted twice

	#pragma omp parallel for reduction(+:wedges)
	for(int i=0;i<g2.num_nodes();i++){
		long double deg_i = g2.degree[i]; 
		if(deg_i==0) continue;
		if(g2.is_upper(i)){
			wedges += (g2_num_v1-1) * deg_i*(deg_i-1)/2; // is this correct? I am curious. 
		}else{
			wedges += (g2_num_v2-1) * deg_i*(deg_i-1)/2; 
		}
	}

	m3 = caterpillar - 4 * m4; 

	long double m21 = wedges - 4*m4 - 2*m3;

	long double m22 = chopsticks - 2*m4 - m3;

	m2 = m21 + m22; 

	m1 = g2_num_edges*(g2_num_v1-1)*(g2_num_v2-1) -4*m4 - 3*m3 -2*m2; 

	m0 = (g2_num_v1*(g2_num_v1-1)/2 )*(g2_num_v2*(g2_num_v2-1)/2)- m4 - m3 - m2 - m1;

	cout<<"m4 = "<<m4 <<endl;
	cout<<"m3 = "<<m3 <<endl;
	cout<<"m2 = "<<m2 <<endl;
	cout<<"m1 = "<<m1 <<endl;
	cout<<"m0 = "<<m0 <<endl;
}

long double one_round_btf(BiGraph& g, unsigned long seed){

	long double t0 = omp_get_wtime();

	std::mt19937 rng(std::random_device{}());

	BiGraph g2(g);

	// user side:
	construct_noisy_graph(g, g2, seed);

	long double t1 = omp_get_wtime();

	RR_time += t1-t0; 

	// if(eva_comm) return 0; // if evaluating communication cost.

	// server side: 
	long double m4;

	if(sampling_noisy_graph){
		cout<<"sampling ratio = "<<p____ <<endl;
		long double sum__=0;
		int num_itr = 1;
		for(int xxx=0;xxx<num_itr;xxx++){
			// get a sampled subgraph: 
			BiGraph g__(g); // sampled subgraph 
			init_genrand(rng());
			for(int i=0;i<g2.num_v1;i++){
				for(auto j : g2.neighbor[i]){
					if(genrand_real2() < p____){ 
						g__.addEdge(i,j);
					}
				}
			}
			g__.computePriority();
			cout<<"num edges in sampled noisy graph = "<<g__.num_edges <<endl;
			m4 = BFC_EVP(g__);
			m4 /= (p____ * p____ * p____ * p____); 
			sum__+=m4;
		}
		m4 = sum__/num_itr;

		cout<<"estimated m4' = "<<m4<<endl;

	}else{
		m4 = BFC_EVP(g2);
	}

	cout<<"m4 is ready"<<endl;

	long double t1x = omp_get_wtime();

	// should we estimate cate and other things independently? 

	long double m3=0, m2=0, m1=0, m0=0, estimate, mu = exp(Eps); 

	// BiGraph g__(g); // sampled subgraph 
	// init_genrand(rng());
	// for(int i=0;i<g2.num_v1;i++){
	// 	for(auto j : g2.neighbor[i]){
	// 		if(genrand_real2() < p____){ 
	// 			g__.addEdge(i,j);
	// 		}
	// 	}
	// }
	// g__.computePriority();
	// cout<<"num edges in sampled noisy graph = "<<g__.num_edges <<endl;
	// long double m4___ = BFC_EVP(g__) / (p____ * p____ * p____ * p____); 
	// compute_m3_m2(m4___,m3,m2,m1,m0,g2);

	compute_m3_m2(m4,m3,m2,m1,m0,g2);

	cout<<"computing m3, m2, m1, m0"<<endl;

	if(count_cate){
		cout<<"\tcaterpillar estimation: "<<endl;
		estimate = -4*mu*mu*mu*m4 + mu*mu*(mu*mu + 3)*m3 - 2*mu*(mu*mu + 1)* m2 + (3*mu*mu + 1) * m1 -4*mu* m0; 
		estimate /= ((mu-1)*(mu-1)*(mu-1)*(mu-1));
		naive_estis[iteration] = m4*4+m3; 
	}else{
		cout<<"\tbtf estimation: "<<endl;
		estimate = mu*mu*mu*mu* m4 - mu*mu*mu * m3 + mu*mu* m2 - mu * m1 +  m0; 
		estimate /= ((mu-1)*(mu-1)*(mu-1)*(mu-1));
		naive_estis[iteration] = m4; 
		// the time needed to get m4 on G2 alone.
		naive_server_side += t1x-t1;
	}
	long double t2 = omp_get_wtime();

	if(count_cate){
		// this is for efficiency evaluations
		// cout<<"// count the number of caterpillar on g2 directly."<<endl;
		// long double cater_g2 = get_cate(g2);
		// long double t3 = omp_get_wtime();
		// naive_server_side += t3 - t2;
	}
	
	// record time elapsed.
	// RR_time += t1-t0; // compute the total time from randomized response. 
	server_side_time += t2 - t1;

	return estimate;
}

// the challenge lies in how to compute deg(u, w) = {v \in N(u), v < w } for each u < w combination.
long double BFC_EVP(BiGraph& g)
{	
	long double BTF = 0;
	#pragma omp parallel for reduction(+:BTF)
	for(int u=0;u<g.num_nodes();u++){
		if(g.degree[u]<=1) continue;
		unordered_map<vid_t,int> count_wedge(0);
		for(auto v: g.neighbor[u]){
			// u -> v 
			for(auto w: g.neighbor[v]){
				// u->v->w
				// if((w>v)&(w>u)) {

				// this step is not invoked for g3.
				if(g.com_p(w,v) & g.com_p(w,u)){// this is a lot faster.
					count_wedge[w] = count_wedge[w]+1;
				}
			}
		}
		// long double btf_u = 0;
		for(auto ele:count_wedge){
			if(ele.second>=2){
				BTF += (ele.second-1)*ele.second/2; 
			}
		}
	}
	return BTF;
}

long double get_wedges(BiGraph& g){
	long double wedges=0; 
	for(int i=0;i<g.num_nodes();i++){
		long double deg_i = g.degree[i]; 
		if(deg_i==0) continue;
		if(g.is_upper(i)){
			wedges +=  deg_i*(deg_i-1)/2; 
		}else{
			wedges += deg_i*(deg_i-1)/2; 
		}
	}
	return wedges; 
}

long double get_laplace(long double parameter){
	return stats::rlaplace(0.0, parameter, engine); 
}

long double get_cate(BiGraph& g){
	long double caterpillar = 0;
	for(int i=0;i<g.num_v1;i++){
		for(auto j:g.neighbor[i]){
			caterpillar += (g.degree[i]-1)*(g.degree[j]-1);
		}
	}
	return caterpillar; 
}