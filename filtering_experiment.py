#!/usr/bin/env python3
"""
Experiment: How much do various filtering/clipping strategies reduce MSE?

Modes tested:
  1. No filtering (baseline)
  2. Public degree filter (upper bound for degree-based filtering)
  3. Public cn filter (absolute upper bound — oracle knows true common neighbors)
  4. Noisy degree filter (realistic, uses eps0 budget)
  5. Positive-part estimator (clip local_res to max(0, ...))
  6. f_u_w threshold (post-filter: skip if f_u_w < K)
  7. Public degree + positive-part (combined)
  8. Winsorization (clip extreme local_res values)
"""

import numpy as np
from collections import defaultdict
import sys
from math import comb, exp

def parse_bipartite_graph(filename):
    upper_adj = defaultdict(set)
    lower_adj = defaultdict(set)
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                u, v = int(parts[0]), int(parts[1])
                upper_adj[u].add(v)
                lower_adj[v].add(u)
    return upper_adj, lower_adj

def compute_local_res_python(K, f, s2):
    f2 = f*f; f3 = f2*f; f4 = f3*f; s4 = s2*s2
    if K == 1: return f
    if K == 2: return f2/2.0 - 0.5*f - s2/2.0
    if K == 3: return f3/6.0 - 0.5*f2 - f*s2/2.0 + f/3.0 + 0.5*s2
    if K == 4: return f4/24.0 - 0.25*f3 - f2*s2/4.0 + 0.458333*f2 + 0.75*f*s2 - 0.25*f + s4/8.0 - 0.458333*s2
    return 0

def simulate_one_round(upper_adj, lower_adj, epsilon, K, mode='none', param=None):
    """
    Simulate one round. Returns estimate and metadata.
    """
    # Different eps0 allocations depending on mode
    if mode == 'noisy_deg_15' or mode == 'noisy_deg_w95':
        eps0 = epsilon * 0.15
    elif mode == 'noisy_deg_25':
        eps0 = epsilon * 0.25
    else:
        eps0 = epsilon * 0.05
    
    eps1 = epsilon * 0.6
    eps2 = epsilon - eps1 - eps0
    
    mu = exp(eps1)
    p_flip = 1.0 / (mu + 1.0)
    gamma = 1.0 / (1.0 - 2.0 * p_flip)
    
    # Noisy degree estimates
    deg_estis = {}
    true_degrees = {}
    for v in upper_adj:
        true_degrees[v] = len(upper_adj[v])
        deg_estis[v] = true_degrees[v] + np.random.laplace(0, 1.0/eps0)
    for v in lower_adj:
        true_degrees[v] = len(lower_adj[v])
        deg_estis[v] = true_degrees[v] + np.random.laplace(0, 1.0/eps0)
    
    # Pick smaller partition
    upper_verts = sorted(upper_adj.keys())
    lower_verts = sorted(lower_adj.keys())
    if len(upper_verts) <= len(lower_verts):
        enum_verts = upper_verts; adj = upper_adj
    else:
        enum_verts = lower_verts; adj = lower_adj
    
    total_estimate = 0.0
    num_pairs_total = 0
    num_pairs_used = 0
    all_local_res = []  # for winsorization
    all_pair_info = []  # for winsorization
    
    for i in range(len(enum_verts)):
        u = enum_verts[i]
        for j in range(i+1, len(enum_verts)):
            w = enum_verts[j]
            num_pairs_total += 1
            true_cn = len(adj[u] & adj[w])
            
            # === PRE-FILTERS (skip before computing f_u_w) ===
            if mode == 'public_deg':
                if min(true_degrees[u], true_degrees[w]) < param:
                    continue
            elif mode == 'public_cn':
                if true_cn < K:
                    continue
            elif mode in ('noisy_deg', 'noisy_deg_15', 'noisy_deg_25', 'noisy_deg_w95'):
                if min(deg_estis[u], deg_estis[w]) < param:
                    continue
            
            # === Compute f_u_w ===
            deg_u = true_degrees[u]
            f_u_w = 0.0
            for nb in adj[u]:
                true_edge_nb_w = 1 if nb in adj[w] else 0
                if np.random.random() < p_flip:
                    noisy_edge = 1 - true_edge_nb_w
                else:
                    noisy_edge = true_edge_nb_w
                if noisy_edge:
                    f_u_w += gamma
                else:
                    f_u_w += (-p_flip) / (1 - 2*p_flip)
            f_u_w += np.random.laplace(0, gamma / eps2)
            
            sigma2 = 2 * gamma**2 / eps2**2 + p_flip * (1-p_flip) * deg_u / (1-2*p_flip)**2
            
            # === Compute local_res ===
            local_res = compute_local_res_python(K, f_u_w, sigma2)
            
            # === Winsorization modes: defer summation ===
            if mode in ('winsorize', 'noisy_deg_w95'):
                all_local_res.append(local_res)
                continue
            
            num_pairs_used += 1
            total_estimate += local_res
    
    # Winsorization: clip at percentile
    if mode in ('winsorize', 'noisy_deg_w95') and all_local_res:
        arr = np.array(all_local_res)
        pct = param if mode == 'winsorize' else 95
        lo = np.percentile(arr, 100 - pct)
        hi = np.percentile(arr, pct)
        arr_clipped = np.clip(arr, lo, hi)
        total_estimate = np.sum(arr_clipped)
        num_pairs_used = len(arr_clipped)
    
    return {
        'estimate': total_estimate,
        'num_pairs_total': num_pairs_total,
        'num_pairs_used': num_pairs_used,
    }


def run_experiment(graph_file, K, epsilon, num_rounds=30):
    print(f"Loading graph: {graph_file}")
    upper_adj, lower_adj = parse_bipartite_graph(graph_file)
    print(f"Upper: {len(upper_adj)}, Lower: {len(lower_adj)}")
    print(f"K={K}, epsilon={epsilon}, rounds={num_rounds}")
    
    if len(upper_adj) <= len(lower_adj):
        enum_verts = sorted(upper_adj.keys()); adj = upper_adj
    else:
        enum_verts = sorted(lower_adj.keys()); adj = lower_adj
    
    true_count = 0
    cn_dist = defaultdict(int)
    for i in range(len(enum_verts)):
        for j in range(i+1, len(enum_verts)):
            cn = len(adj[enum_verts[i]] & adj[enum_verts[j]])
            cn_dist[cn] += 1
            true_count += comb(cn, K)
    
    total_pairs = sum(cn_dist.values())
    signal_pairs = sum(v for k,v in cn_dist.items() if k >= K)
    
    print(f"\nTrue (2,{K})-biclique count: {true_count}")
    print(f"Total pairs: {total_pairs}, Signal: {signal_pairs} ({100*signal_pairs/total_pairs:.1f}%), Noise: {total_pairs-signal_pairs} ({100*(total_pairs-signal_pairs)/total_pairs:.1f}%)")
    print(f"CN distribution: {dict(sorted(cn_dist.items()))}")
    
    modes = [
        ('none',           None,  'Baseline (no filtering)'),
        ('public_cn',      None,  'Oracle: public cn filter (absolute upper bound)'),
        ('public_deg',     K,     f'Public degree filter (min_deg >= {K})'),
        ('noisy_deg',      K,     f'Noisy degree filter (eps0=5%, deg_est >= {K})'),
        ('noisy_deg_15',   K,     f'Noisy degree filter (eps0=15%, deg_est >= {K})'),
        ('noisy_deg_25',   K,     f'Noisy degree filter (eps0=25%, deg_est >= {K})'),
        ('winsorize',      95,    'Winsorize at 5th/95th percentile'),
        ('noisy_deg_w95',  K,     f'Noisy deg (eps0=15%) + winsorize 95th'),
    ]
    
    print(f"\n{'='*80}")
    print(f"{'Mode':<55} {'Pairs%':>7} {'Bias':>10} {'StdDev':>10} {'RMSE':>10} {'RelRMSE':>8}")
    print(f"{'='*80}")
    
    for mode, param, description in modes:
        estimates = []
        pairs_pcts = []
        
        for r in range(num_rounds):
            np.random.seed(r * 1000 + 42)
            result = simulate_one_round(upper_adj, lower_adj, epsilon, K,
                                       mode=mode, param=param)
            estimates.append(result['estimate'])
            pairs_pcts.append(100 * result['num_pairs_used'] / result['num_pairs_total'] if result['num_pairs_total'] > 0 else 0)
        
        estimates = np.array(estimates)
        mean_est = np.mean(estimates)
        bias = mean_est - true_count
        std_est = np.std(estimates)
        rmse = np.sqrt(np.mean((estimates - true_count)**2))
        avg_pct = np.mean(pairs_pcts)
        rel_rmse = f"{rmse/true_count*100:.1f}%" if true_count > 0 else "N/A"
        
        print(f"  {description:<53} {avg_pct:>6.1f}% {bias:>10.1f} {std_est:>10.1f} {rmse:>10.1f} {rel_rmse:>8}")
    
    print(f"{'='*80}")


if __name__ == "__main__":
    graph_file = sys.argv[1] if len(sys.argv) > 1 else "/data1/yizhangh/bidata/co.e"
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    epsilon = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
    num_rounds = int(sys.argv[4]) if len(sys.argv) > 4 else 30
    
    run_experiment(graph_file, K, epsilon, num_rounds)
