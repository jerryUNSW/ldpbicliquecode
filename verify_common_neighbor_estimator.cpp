/**
 * Verification: Unbiased estimation of Var(S'(u,v)) from G'
 *
 * === Grand Scheme ===
 *
 * We are exploring a new paradigm for (p,q)-biclique counting under Edge LDP.
 * Instead of directly modeling the transformation probability of different motifs
 * (as in the existing pqbiclique framework), we take a building-block approach:
 *
 *   1. Under the one-round framework, use the noisy graph G' to estimate S(u,v),
 *      the number of common neighbors shared by two vertices u, v in U(G).
 *   2. Use S'(u,v) to estimate the number of (2,q)-bicliques involving u and v,
 *      which equals C(S(u,v), q).
 *
 * Starting with p=2 (i.e., pairs of vertices in U), the key challenge is:
 *   - S'(u,v) is an unbiased estimator of S(u,v), easy to construct from G'.
 *   - But C(S', q) is NOT unbiased for C(S, q) due to nonlinearity.
 *   - For q=2: E[S'(S'-1)/2] = C(S,2) + Var(S')/2
 *   - So we need a good (unbiased) estimate of Var(S') to correct the bias.
 *
 * This file verifies that we CAN estimate Var(S') unbiasedly from G'.
 *
 * === Method ===
 *
 * For a pair (u,v) in U, each vertex w in L falls into one of four categories:
 *   - (a,b) = (A(u,w), A(v,w)) in {(0,0), (1,0), (0,1), (1,1)}
 * Let n_ab = number of w in each category. Then:
 *
 *   Var(S'(u,v)) = n00 * C1 + (n10 + n01) * C2 + n11 * C3
 *
 * where C1, C2, C3 are known constants depending only on epsilon:
 *   - p = 1/(e^eps + 1)           (flipping probability)
 *   - delta = 1 - 2p              (signal strength)
 *   - sigma^2 = p(1-p)
 *   - C1 = [sigma^2]^2 / delta^4                          (case 0,0)
 *   - C2 = sigma^2 * (sigma^2 + delta^2) / delta^4        (case 1,0 or 0,1)
 *   - C3 = (sigma^2 + delta^2)^2 / delta^4 - 1            (case 1,1)
 *
 * Since the variance formula is LINEAR in n00, n10, n01, n11, plugging in
 * unbiased estimators of these counts yields an unbiased estimator of Var(S').
 *
 * The n-count estimators are constructed by iterating over w in L(G):
 *   - a_hat(w) = (A'(u,w) - p) / delta       (unbiased for A(u,w))
 *   - b_hat(w) = (A'(v,w) - p) / delta       (unbiased for A(v,w))
 *   - n11_hat = sum_w a_hat * b_hat
 *   - n10_hat = sum_w a_hat * (1 - b_hat)
 *   - n01_hat = sum_w (1 - a_hat) * b_hat
 *   - n00_hat = sum_w (1 - a_hat) * (1 - b_hat)
 *
 * The butterfly (2,2)-biclique estimator for pair (u,v) is then:
 *   C(S,2)_hat = S'(S'-1)/2 - Var_hat(S')/2
 *
 * === What We Verify (Monte Carlo) ===
 *
 *   1. E[S'(u,v)] = S(u,v)                    (S' is unbiased)
 *   2. Empirical Var(S') matches formula       (variance formula is correct)
 *   3. E[Var_hat(S')] = Var(S')                (variance estimator is unbiased)
 *   4. E[C(S,2)_hat] = C(S,2)                 (butterfly estimator is unbiased)
 */

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>
#include <cassert>

using namespace std;

struct BipartiteGraph {
    int nU, nL;                          // |U|, |L|
    vector<vector<bool>> adj;            // adj[u][w] = true if edge (u,w) exists
    
    BipartiteGraph(int nU, int nL) : nU(nU), nL(nL), adj(nU, vector<bool>(nL, false)) {}
    
    void addEdge(int u, int w) {
        assert(u >= 0 && u < nU && w >= 0 && w < nL);
        adj[u][w] = true;
    }
    
    // True common neighbor count S(u,v)
    int commonNeighbors(int u, int v) const {
        int count = 0;
        for (int w = 0; w < nL; w++)
            if (adj[u][w] && adj[v][w]) count++;
        return count;
    }
    
    // Degree of u in L
    int degree(int u) const {
        int d = 0;
        for (int w = 0; w < nL; w++)
            if (adj[u][w]) d++;
        return d;
    }
    
    // Count n00, n10, n01, n11 for pair (u,v)
    void countCategories(int u, int v, int &n00, int &n10, int &n01, int &n11) const {
        n00 = n10 = n01 = n11 = 0;
        for (int w = 0; w < nL; w++) {
            bool a = adj[u][w], b = adj[v][w];
            if (!a && !b) n00++;
            else if (a && !b) n10++;
            else if (!a && b) n01++;
            else n11++;
        }
    }
};

// Generate noisy graph G' via randomized response
BipartiteGraph noisyGraph(const BipartiteGraph &G, double p, mt19937 &rng) {
    BipartiteGraph Gp(G.nU, G.nL);
    uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (int u = 0; u < G.nU; u++) {
        for (int w = 0; w < G.nL; w++) {
            bool real_edge = G.adj[u][w];
            double r = dist(rng);
            if (real_edge) {
                Gp.adj[u][w] = (r >= p);   // keep with prob 1-p, flip with prob p
            } else {
                Gp.adj[u][w] = (r < p);    // flip with prob p
            }
        }
    }
    return Gp;
}

// Compute S'(u,v) from noisy graph
double estimateCommonNeighbors(const BipartiteGraph &Gp, int u, int v, double p, double delta) {
    double S_hat = 0.0;
    for (int w = 0; w < Gp.nL; w++) {
        double a_hat = (Gp.adj[u][w] ? 1.0 : 0.0);
        double b_hat = (Gp.adj[v][w] ? 1.0 : 0.0);
        double phi = (a_hat - p) * (b_hat - p) / (delta * delta);
        S_hat += phi;
    }
    return S_hat;
}

// Compute estimated n-counts from noisy graph
void estimateNCounts(const BipartiteGraph &Gp, int u, int v, double p, double delta,
                     double &n00_hat, double &n10_hat, double &n01_hat, double &n11_hat) {
    n00_hat = n10_hat = n01_hat = n11_hat = 0.0;
    for (int w = 0; w < Gp.nL; w++) {
        double a_hat = ((Gp.adj[u][w] ? 1.0 : 0.0) - p) / delta;  // unbiased for A(u,w)
        double b_hat = ((Gp.adj[v][w] ? 1.0 : 0.0) - p) / delta;  // unbiased for A(v,w)
        
        n11_hat += a_hat * b_hat;
        n10_hat += a_hat * (1.0 - b_hat);
        n01_hat += (1.0 - a_hat) * b_hat;
        n00_hat += (1.0 - a_hat) * (1.0 - b_hat);
    }
}

// Compute true Var(S') from the formula using true n-counts
double trueVariance(int n00, int n10, int n01, int n11, double p, double delta) {
    double sigma2 = p * (1.0 - p);
    double d4 = delta * delta * delta * delta;
    double d2 = delta * delta;
    
    double C1 = (sigma2 * sigma2) / d4;                          // case (0,0)
    double C2 = sigma2 * (sigma2 + d2) / d4;                     // case (1,0) or (0,1)
    double C3 = (sigma2 + d2) * (sigma2 + d2) / d4 - 1.0;       // case (1,1)
    
    return n00 * C1 + (n10 + n01) * C2 + n11 * C3;
}

// Compute estimated Var(S') by plugging in estimated n-counts
double estimateVariance(double n00_hat, double n10_hat, double n01_hat, double n11_hat,
                        double p, double delta) {
    double sigma2 = p * (1.0 - p);
    double d4 = delta * delta * delta * delta;
    double d2 = delta * delta;
    
    double C1 = (sigma2 * sigma2) / d4;
    double C2 = sigma2 * (sigma2 + d2) / d4;
    double C3 = (sigma2 + d2) * (sigma2 + d2) / d4 - 1.0;
    
    return n00_hat * C1 + (n10_hat + n01_hat) * C2 + n11_hat * C3;
}

int main() {
    // === Build a small bipartite graph ===
    // U = {0, 1, 2}, L = {0, 1, 2, 3, 4}
    int nU = 3, nL = 5;
    BipartiteGraph G(nU, nL);
    
    // Edges for u=0: neighbors in L = {0, 1, 2, 3}
    G.addEdge(0, 0); G.addEdge(0, 1); G.addEdge(0, 2); G.addEdge(0, 3);
    // Edges for u=1: neighbors in L = {0, 1, 4}
    G.addEdge(1, 0); G.addEdge(1, 1); G.addEdge(1, 4);
    // Edges for u=2: neighbors in L = {1, 2, 3, 4}
    G.addEdge(2, 1); G.addEdge(2, 2); G.addEdge(2, 3); G.addEdge(2, 4);
    
    // Test pair: u=0, v=1
    // Common neighbors: {0, 1} → S = 2
    // n11=2, n10=2 (w=2,3), n01=1 (w=4), n00=0
    
    vector<pair<int,int>> test_pairs = {{0, 1}, {0, 2}, {1, 2}};
    vector<double> epsilons = {0.5, 1.0, 2.0, 3.0};
    
    int num_trials = 100000;
    mt19937 rng(42);
    
    cout << fixed << setprecision(6);
    cout << "=== Verification of Common Neighbor & Variance Estimators ===" << endl;
    cout << "Bipartite graph: |U|=" << nU << ", |L|=" << nL << endl;
    cout << "Number of trials: " << num_trials << endl;
    cout << endl;
    
    for (auto [ui, vi] : test_pairs) {
        int true_S = G.commonNeighbors(ui, vi);
        int n00, n10, n01, n11;
        G.countCategories(ui, vi, n00, n10, n01, n11);
        double true_C2 = true_S * (true_S - 1) / 2.0;  // C(S, 2)
        
        cout << "========================================" << endl;
        cout << "Pair (u=" << ui << ", v=" << vi << ")" << endl;
        cout << "True S(u,v) = " << true_S << endl;
        cout << "True C(S,2) = " << true_C2 << endl;
        cout << "n00=" << n00 << " n10=" << n10 << " n01=" << n01 << " n11=" << n11 << endl;
        cout << "d_u=" << G.degree(ui) << " d_v=" << G.degree(vi) << endl;
        cout << "----------------------------------------" << endl;
        
        for (double eps : epsilons) {
            double p = 1.0 / (exp(eps) + 1.0);
            double delta = 1.0 - 2.0 * p;
            double true_var = trueVariance(n00, n10, n01, n11, p, delta);
            
            // Run Monte Carlo
            double sum_S = 0, sum_S2 = 0;
            double sum_var_hat = 0;
            double sum_butterfly_hat = 0;
            
            for (int t = 0; t < num_trials; t++) {
                BipartiteGraph Gp = noisyGraph(G, p, rng);
                
                // Estimate S'(u,v)
                double S_hat = estimateCommonNeighbors(Gp, ui, vi, p, delta);
                sum_S += S_hat;
                sum_S2 += S_hat * S_hat;
                
                // Estimate n-counts
                double n00_h, n10_h, n01_h, n11_h;
                estimateNCounts(Gp, ui, vi, p, delta, n00_h, n10_h, n01_h, n11_h);
                
                // Estimate Var(S')
                double var_hat = estimateVariance(n00_h, n10_h, n01_h, n11_h, p, delta);
                sum_var_hat += var_hat;
                
                // Butterfly estimator: C(S,2)_hat = S'(S'-1)/2 - Var_hat/2
                double butterfly_hat = S_hat * (S_hat - 1.0) / 2.0 - var_hat / 2.0;
                sum_butterfly_hat += butterfly_hat;
            }
            
            double mean_S = sum_S / num_trials;
            double empirical_var = sum_S2 / num_trials - mean_S * mean_S;
            double mean_var_hat = sum_var_hat / num_trials;
            double mean_butterfly = sum_butterfly_hat / num_trials;
            
            cout << "  epsilon = " << eps 
                 << "  (p=" << setprecision(4) << p 
                 << ", delta=" << delta << ")" << setprecision(6) << endl;
            
            // Check 1: E[S'] ≈ S
            cout << "    E[S']          = " << setw(12) << mean_S 
                 << "   (true S = " << true_S << ")" 
                 << "   err=" << setprecision(4) << abs(mean_S - true_S) << setprecision(6) << endl;
            
            // Check 2: Empirical Var(S') ≈ formula Var(S')
            cout << "    EmpVar(S')     = " << setw(12) << empirical_var 
                 << "   (formula = " << true_var << ")"
                 << "   err=" << setprecision(4) << abs(empirical_var - true_var) << setprecision(6) << endl;
            
            // Check 3: E[Var_hat(S')] ≈ Var(S')
            cout << "    E[Var_hat(S')] = " << setw(12) << mean_var_hat 
                 << "   (true   = " << true_var << ")"
                 << "   err=" << setprecision(4) << abs(mean_var_hat - true_var) << setprecision(6) << endl;
            
            // Check 4: E[butterfly_hat] ≈ C(S,2)
            cout << "    E[C(S,2)_hat]  = " << setw(12) << mean_butterfly 
                 << "   (true   = " << true_C2 << ")"
                 << "   err=" << setprecision(4) << abs(mean_butterfly - true_C2) << setprecision(6) << endl;
            
            cout << endl;
        }
    }
    
    // === Larger random graph test ===
    cout << "========================================" << endl;
    cout << "=== Larger Random Bipartite Graph ===" << endl;
    int nU2 = 10, nL2 = 20;
    double edge_prob = 0.3;
    BipartiteGraph G2(nU2, nL2);
    
    uniform_real_distribution<double> dist(0.0, 1.0);
    for (int u = 0; u < nU2; u++)
        for (int w = 0; w < nL2; w++)
            if (dist(rng) < edge_prob)
                G2.addEdge(u, w);
    
    // Test a few pairs
    vector<pair<int,int>> pairs2 = {{0,1}, {2,5}, {3,7}};
    double eps = 1.0;
    double p = 1.0 / (exp(eps) + 1.0);
    double delta = 1.0 - 2.0 * p;
    
    cout << "|U|=" << nU2 << ", |L|=" << nL2 << ", edge_prob=" << edge_prob << endl;
    cout << "epsilon=" << eps << ", p=" << p << ", delta=" << delta << endl;
    cout << "Trials: " << num_trials << endl << endl;
    
    for (auto [ui, vi] : pairs2) {
        int true_S = G2.commonNeighbors(ui, vi);
        int n00, n10, n01, n11;
        G2.countCategories(ui, vi, n00, n10, n01, n11);
        double true_var = trueVariance(n00, n10, n01, n11, p, delta);
        double true_C2 = true_S * (true_S - 1) / 2.0;
        
        double sum_S = 0, sum_S2 = 0, sum_var_hat = 0, sum_butterfly_hat = 0;
        
        for (int t = 0; t < num_trials; t++) {
            BipartiteGraph Gp = noisyGraph(G2, p, rng);
            double S_hat = estimateCommonNeighbors(Gp, ui, vi, p, delta);
            sum_S += S_hat;
            sum_S2 += S_hat * S_hat;
            
            double n00_h, n10_h, n01_h, n11_h;
            estimateNCounts(Gp, ui, vi, p, delta, n00_h, n10_h, n01_h, n11_h);
            double var_hat = estimateVariance(n00_h, n10_h, n01_h, n11_h, p, delta);
            sum_var_hat += var_hat;
            
            double butterfly_hat = S_hat * (S_hat - 1.0) / 2.0 - var_hat / 2.0;
            sum_butterfly_hat += butterfly_hat;
        }
        
        double mean_S = sum_S / num_trials;
        double empirical_var = sum_S2 / num_trials - mean_S * mean_S;
        double mean_var_hat = sum_var_hat / num_trials;
        double mean_butterfly = sum_butterfly_hat / num_trials;
        
        cout << "Pair (u=" << ui << ", v=" << vi << "): S=" << true_S 
             << ", C(S,2)=" << true_C2 << endl;
        cout << "  E[S']=" << mean_S << " (true=" << true_S << ")" << endl;
        cout << "  EmpVar=" << empirical_var << " Formula=" << true_var 
             << " E[Var_hat]=" << mean_var_hat << endl;
        cout << "  E[C(S,2)_hat]=" << mean_butterfly << " (true=" << true_C2 << ")" << endl;
        cout << endl;
    }
    
    cout << "=== DONE ===" << endl;
    return 0;
}
