#!/usr/bin/env python3
"""
Variance Analysis for Multi-Round Biclique Counting
Analyzes which (u,w) vertex pairs contribute most to variance
"""

import sys
import os
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def parse_edge_list(filename):
    """Parse bipartite graph from edge list"""
    edges = []
    upper_vertices = set()
    lower_vertices = set()
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    u, v = int(parts[0]), int(parts[1])
                    edges.append((u, v))
                    upper_vertices.add(u)
                    lower_vertices.add(v)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None, None, None
    
    return edges, upper_vertices, lower_vertices

def analyze_pair_variance(edges, upper_vertices, lower_vertices, num_rounds=10):
    """
    Analyze variance contribution of each (u,w) pair
    For each pair, simulate multiple rounds and track variance
    """
    
    print("\n=== Variance Analysis for (u,w) Pairs ===")
    print(f"Upper vertices: {len(upper_vertices)}")
    print(f"Lower vertices: {len(lower_vertices)}")
    print(f"Edges: {len(edges)}")
    print(f"Running {num_rounds} rounds\n")
    
    # Build adjacency lists
    upper_neighbors = defaultdict(set)
    lower_neighbors = defaultdict(set)
    
    for u, v in edges:
        upper_neighbors[u].add(v)
        lower_neighbors[v].add(u)
    
    # For each (u,w) pair, estimate butterflies across rounds
    pair_estimates = defaultdict(list)
    
    for round_num in range(num_rounds):
        print(f"Round {round_num + 1}/{num_rounds}...", end=' ', flush=True)
        
        # For each upper vertex u
        for u in upper_vertices:
            neighbors_u = upper_neighbors[u]
            
            # For each lower vertex w
            for w in lower_vertices:
                neighbors_w = lower_neighbors[w]
                
                # Count common neighbors (butterflies involving this pair)
                # A butterfly (2,2)-biclique is: u-v1-w-v2-u
                # Common neighbors of u and w in opposite partition
                common = len(neighbors_u & neighbors_w)
                
                # Estimate with noise (simulate DP)
                noise = np.random.laplace(0, 1.0)
                estimate = max(0, common + noise)
                
                pair_estimates[(u, w)].append(estimate)
        
        print("done")
    
    # Compute variance for each pair
    pair_variances = []
    
    for (u, w), estimates in pair_estimates.items():
        if len(estimates) < 2:
            continue
        
        mean = np.mean(estimates)
        variance = np.var(estimates, ddof=1)
        std_dev = np.std(estimates, ddof=1)
        
        pair_variances.append({
            'u': u,
            'w': w,
            'variance': variance,
            'mean': mean,
            'std_dev': std_dev,
            'estimates': estimates
        })
    
    # Sort by variance (descending)
    pair_variances.sort(key=lambda x: x['variance'], reverse=True)
    
    # Print top contributors
    print("\n=== Top 20 Vertex Pairs by Variance ===")
    print(f"{'Rank':<6} {'u':<8} {'w':<8} {'Variance':<15} {'Std Dev':<15} {'Mean Est':<15}")
    print("-" * 70)
    
    for i, pv in enumerate(pair_variances[:20]):
        print(f"{i+1:<6} {pv['u']:<8} {pv['w']:<8} {pv['variance']:<15.6f} {pv['std_dev']:<15.6f} {pv['mean']:<15.6f}")
    
    # Statistics
    print("\n=== Variance Statistics ===")
    total_variance = sum(pv['variance'] for pv in pair_variances)
    max_variance = max(pv['variance'] for pv in pair_variances)
    min_variance = min(pv['variance'] for pv in pair_variances)
    avg_variance = total_variance / len(pair_variances)
    
    print(f"Total pairs analyzed: {len(pair_variances)}")
    print(f"Average variance: {avg_variance:.6f}")
    print(f"Max variance: {max_variance:.6f}")
    print(f"Min variance: {min_variance:.6f}")
    
    # Top 10% contribution
    top_10_percent = max(1, len(pair_variances) // 10)
    top_10_variance = sum(pv['variance'] for pv in pair_variances[:top_10_percent])
    
    print(f"\nTop 10% pairs ({top_10_percent} pairs) contribute: {top_10_variance/total_variance*100:.2f}% of total variance")
    
    # Save detailed results
    with open('pair_variance_analysis.txt', 'w') as f:
        f.write("Vertex Pair Variance Analysis\n")
        f.write(f"Upper vertices: {len(upper_vertices)}, Lower vertices: {len(lower_vertices)}\n")
        f.write(f"Edges: {len(edges)}\n")
        f.write(f"Rounds: {num_rounds}\n\n")
        f.write("u\tw\tVariance\tStdDev\tMeanEst\n")
        
        for pv in pair_variances:
            f.write(f"{pv['u']}\t{pv['w']}\t{pv['variance']:.6f}\t{pv['std_dev']:.6f}\t{pv['mean']:.6f}\n")
    
    print("\nDetailed results saved to: pair_variance_analysis.txt")
    
    return pair_variances

def suggest_optimizations(pair_variances, upper_vertices, lower_vertices):
    """Suggest optimizations based on variance analysis"""
    
    print("\n=== Optimization Suggestions ===")
    
    # Identify high-variance pairs
    high_var_pairs = pair_variances[:max(1, len(pair_variances) // 10)]
    high_var_u = set(pv['u'] for pv in high_var_pairs)
    high_var_w = set(pv['w'] for pv in high_var_pairs)
    
    print(f"\n1. High-variance upper vertices ({len(high_var_u)} vertices):")
    for u in sorted(high_var_u)[:10]:
        print(f"   u={u}")
    
    print(f"\n2. High-variance lower vertices ({len(high_var_w)} vertices):")
    for w in sorted(high_var_w)[:10]:
        print(f"   w={w}")
    
    print("\n3. Recommended optimizations:")
    print("   a) Allocate more privacy budget to high-variance pairs")
    print("   b) Use stratified sampling: separate high/low variance pairs")
    print("   c) Apply variance reduction (control variates, importance sampling)")
    print("   d) Adaptive epsilon allocation based on pair variance")
    print("   e) Focus refinement rounds on high-variance pairs")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 variance_analysis.py <graph_file> [num_rounds]")
        print("  graph_file: path to bipartite graph edge list")
        print("  num_rounds: number of rounds for variance analysis (default: 10)")
        sys.exit(1)
    
    graph_file = sys.argv[1]
    num_rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    
    if not os.path.exists(graph_file):
        print(f"Error: File not found: {graph_file}")
        sys.exit(1)
    
    print(f"Loading graph: {graph_file}")
    edges, upper_vertices, lower_vertices = parse_edge_list(graph_file)
    
    if edges is None:
        print("Error: Failed to parse graph")
        sys.exit(1)
    
    # Run variance analysis
    pair_variances = analyze_pair_variance(edges, upper_vertices, lower_vertices, num_rounds)
    
    # Suggest optimizations
    suggest_optimizations(pair_variances, upper_vertices, lower_vertices)
    
    print("\n=== Analysis Complete ===")

if __name__ == "__main__":
    main()
