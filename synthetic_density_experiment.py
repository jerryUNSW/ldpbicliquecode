#!/usr/bin/env python3
"""
Synthetic Density Experiment
Creates very dense synthetic graphs and tests algorithm performance across extreme density ranges
"""

import os
import sqlite3
import subprocess
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def create_synthetic_graph(n1, n2, density, output_prefix):
    """Create a synthetic bipartite graph with specified density"""
    
    # Create output directory
    os.makedirs('density-exp/synthetic', exist_ok=True)
    
    # Calculate number of edges
    total_possible_edges = n1 * n2
    num_edges = int(total_possible_edges * density)
    
    # Generate random edges
    edges = set()
    while len(edges) < num_edges:
        u = random.randint(0, n1 - 1)
        v = random.randint(0, n2 - 1)
        edges.add((u, v))
    
    # Write edge file
    edge_file = f'density-exp/synthetic/{output_prefix}.e'
    with open(edge_file, 'w') as f:
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")
    
    # Write metadata file
    meta_file = f'density-exp/synthetic/{output_prefix}.meta'
    with open(meta_file, 'w') as f:
        f.write(f"{n1}\n")
        f.write(f"{n2}\n")
        f.write(f"{len(edges)}\n")
    
    return edge_file, meta_file, len(edges)

def create_database():
    """Create database for synthetic density experiment results"""
    conn = sqlite3.connect('density-exp/synthetic_density_experiment.db')
    cursor = conn.cursor()
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS synthetic_density_experiment (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            graph_type TEXT,
            n1 INTEGER,
            n2 INTEGER,
            algorithm TEXT,
            p INTEGER,
            q INTEGER,
            epsilon REAL,
            edge_density REAL,
            total_edges INTEGER,
            ground_truth REAL,
            estimate REAL,
            relative_error REAL,
            T_samples INTEGER,
            run_id INTEGER,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    conn.commit()
    conn.close()
    print("Synthetic database created successfully!")

def run_algorithm_on_synthetic(graph_path, algorithm, p, q, epsilon, run_id):
    """Run a single algorithm on synthetic graph"""
    
    # Map algorithm names to numbers
    alg_map = {
        'Naive': 0,
        'oneR': 1,
        'ADV': 2,
        'ADV+': 3,
        'ADV++': 4
    }
    
    if algorithm not in alg_map:
        return {'success': False, 'error': f'Unknown algorithm: {algorithm}'}
    
    alg_num = alg_map[algorithm]
    
    try:
        # Run biclique algorithm
        biclique_path = '/data/yizhangh/ldp-pq/biclique'
        cmd = [biclique_path, str(epsilon), graph_path, '1', str(alg_num), str(p), str(q)]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode != 0:
            return {'success': False, 'error': f'Algorithm failed: {result.stderr}'}
        
        # Parse output
        output_lines = result.stdout.strip().split('\n')
        
        # Extract relative error and ground truth
        rel_error = None
        ground_truth = None
        estimate = None
        
        for line in output_lines:
            if 'adv rel err' in line:
                rel_error = float(line.split('=')[-1].strip())
            elif 'real count =' in line:
                ground_truth = float(line.split('=')[-1].strip())
            elif 'estimate =' in line and 'Mean' not in line:
                estimate = float(line.split('=')[-1].strip())
        
        if rel_error is None or ground_truth is None:
            return {'success': False, 'error': 'Could not parse results'}
        
        return {
            'success': True,
            'relative_error': rel_error,
            'ground_truth': ground_truth,
            'estimate': estimate
        }
        
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Algorithm timed out'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def run_synthetic_experiment():
    """Run the synthetic density experiment for multiple graph sizes"""
    
    # Create output directory and database
    os.makedirs('density-exp/synthetic', exist_ok=True)
    create_database()
    
    # Experiment parameters - test multiple graph sizes
    graph_sizes = [
        (50, 100),    # Small
        (100, 200),   # Medium (original)
        (150, 300),   # Large
        (200, 400)    # Very large
    ]
    algorithms = ['Naive', 'oneR', 'ADV++']
    p_values = [2, 3]  # Test P=2 and P=3
    q_values = [3, 4, 5]  # Test Q=3, 4, 5
    epsilon = 1.0
    
    # Density range: 10%, 15%, 20%, ..., 90%, 95%, 100% (skip the broken low densities)
    edge_densities = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    num_graphs_per_density = 5  # Generate 5 different graphs for each density
    num_runs_per_graph = 2  # 2 runs per graph
    
    print(f"Starting synthetic density experiment")
    print(f"Graph sizes: {graph_sizes}")
    print(f"Algorithms: {algorithms}")
    print(f"P values: {p_values}")
    print(f"Q values: {q_values}")
    print(f"Epsilon: {epsilon}")
    print(f"Densities: {edge_densities}")
    print(f"Graphs per density: {num_graphs_per_density}")
    print(f"Runs per graph: {num_runs_per_graph}")
    
    conn = sqlite3.connect('density-exp/synthetic_density_experiment.db')
    cursor = conn.cursor()
    
    total_experiments = len(graph_sizes) * len(p_values) * len(q_values) * len(algorithms) * len(edge_densities) * num_graphs_per_density * num_runs_per_graph
    experiment_count = 0
    
    for n1, n2 in graph_sizes:
        print(f"\n{'='*60}")
        print(f"EXPERIMENTING WITH GRAPH SIZE: {n1} x {n2}")
        print(f"{'='*60}")
        
        for p in p_values:
            for q in q_values:
                print(f"\n--- P={p}, Q={q} ---")
                
                for density in edge_densities:
                    print(f"\n=== Density {density:.2f} ({density:.1%}) ===")
                    
                    for graph_id in range(num_graphs_per_density):
                        print(f"\n  Graph {graph_id+1}/{num_graphs_per_density}:")
                        
                        # Create synthetic graph with unique prefix
                        graph_prefix = f"synthetic_n{n1}x{n2}_P{p}_Q{q}_d{density:.2f}_g{graph_id}"
                        edge_file, meta_file, total_edges = create_synthetic_graph(n1, n2, density, graph_prefix)
                        graph_path = meta_file.replace('.meta', '')
                        
                        print(f"    Created graph: {total_edges} edges ({density:.1%} density)")
                        
                        for algorithm in algorithms:
                            print(f"    Running {algorithm}...")
                            rel_errors = []
                            
                            for run_id in range(num_runs_per_graph):
                                experiment_count += 1
                                print(f"      Run {run_id+1}/{num_runs_per_graph} ({experiment_count}/{total_experiments})")
                                
                                result = run_algorithm_on_synthetic(graph_path, algorithm, p, q, epsilon, run_id)
                                
                                if result['success']:
                                    rel_errors.append(result['relative_error'])
                                    
                                    # Store result in database
                                    cursor.execute('''
                                        INSERT INTO synthetic_density_experiment 
                                        (graph_type, n1, n2, algorithm, p, q, epsilon, edge_density, 
                                         total_edges, ground_truth, estimate, relative_error, T_samples, run_id)
                                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                                    ''', (
                                        'synthetic', n1, n2, algorithm, p, q, epsilon, density,
                                        total_edges, result['ground_truth'], result['estimate'], 
                                        result['relative_error'], -1, run_id
                                    ))
                                else:
                                    print(f"        Failed: {result.get('error', 'Unknown error')}")
                            
                            if rel_errors:
                                mean_rel_error = np.mean(rel_errors)
                                std_rel_error = np.std(rel_errors)
                                print(f"      Mean MRE: {mean_rel_error:.3f} ± {std_rel_error:.3f}")
                        
                        # Clean up temporary files for this graph
                        if os.path.exists(edge_file):
                            os.remove(edge_file)
                        if os.path.exists(meta_file):
                            os.remove(meta_file)
    
    conn.commit()
    conn.close()
    
    print(f"\nSynthetic experiment completed! Results stored in density-exp/synthetic_density_experiment.db")

def plot_synthetic_results():
    """Plot the synthetic density experiment results for each graph size"""
    
    # Load results from database
    conn = sqlite3.connect('density-exp/synthetic_density_experiment.db')
    
    query = '''
    SELECT graph_type, n1, n2, algorithm, p, q, edge_density, 
           AVG(relative_error) as mean_rel_error,
           COUNT(*) as num_runs
    FROM synthetic_density_experiment 
    GROUP BY graph_type, n1, n2, algorithm, p, q, edge_density
    ORDER BY graph_type, n1, n2, algorithm, p, q, edge_density
    '''
    
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    if df.empty:
        print("No results found in database!")
        return
    
    # Algorithm colors and styles (only Naive, oneR, ADV++)
    algorithm_styles = {
        'Naive': {'color': '#1f77b4', 'marker': 'o', 'linestyle': '-'},
        'oneR': {'color': '#ff7f0e', 'marker': 's', 'linestyle': '-'},
        'ADV++': {'color': '#9467bd', 'marker': 'v', 'linestyle': '-'}
    }
    
    # Only plot the algorithms we want (Naive, oneR, ADV++)
    algorithms_to_plot = ['Naive', 'oneR', 'ADV++']
    
    # Get unique combinations of graph sizes and P,Q values
    unique_combinations = df[['n1', 'n2', 'p', 'q']].drop_duplicates().values
    
    print(f"Creating plots for {len(unique_combinations)} different combinations...")
    
    # Create a separate plot for each combination
    for n1, n2, p, q in unique_combinations:
        print(f"\nCreating plot for graph size {n1}×{n2}, P={p}, Q={q}...")
        
        # Filter data for this combination
        size_data = df[(df['n1'] == n1) & (df['n2'] == n2) & (df['p'] == p) & (df['q'] == q)]
        
        if size_data.empty:
            print(f"  No data found for graph size {n1}×{n2}, P={p}, Q={q}")
            continue
        
        # Create plot
        plt.figure(figsize=(14, 10))
        
        for algorithm in algorithms_to_plot:
            if algorithm in size_data['algorithm'].values:
                alg_data = size_data[size_data['algorithm'] == algorithm]
                style = algorithm_styles.get(algorithm, {'color': 'gray', 'marker': 'o', 'linestyle': '-'})
                
                plt.plot(
                    alg_data['edge_density'] * 100,  # Convert to percentage
                    alg_data['mean_rel_error'],
                    label=algorithm,
                    **style,
                    linewidth=2,
                    markersize=8
                )
        
        plt.xlabel('Edge Density (%)', fontsize=24)
        plt.ylabel('Mean Relative Error', fontsize=24)
        plt.title(f'Effect of Graph Density on Algorithm Performance\\nSynthetic Graph {n1}×{n2} (P={p}, Q={q})', fontsize=28)
        plt.legend(fontsize=20)
        plt.grid(True, alpha=0.3)
        plt.yscale('log')
        
        # Set x-axis to show percentages (10%, 15%, 20%, ..., 90%, 95%, 100%)
        plt.xticks([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100], fontsize=18)
        plt.yticks(fontsize=18)
        
        plt.tight_layout()
        plt.savefig(f'density-exp/synthetic_{n1}x{n2}_P{p}_Q{q}_extreme_density.pdf', dpi=300)
        plt.show()
        
        print(f"  Plot saved as density-exp/synthetic_{n1}x{n2}_P{p}_Q{q}_extreme_density.pdf")
        
        # Print summary statistics for this combination
        print(f"  Summary for {n1}×{n2}, P={p}, Q={q}:")
        for algorithm in algorithms_to_plot:
            if algorithm in size_data['algorithm'].values:
                alg_data = size_data[size_data['algorithm'] == algorithm]
                min_error = alg_data['mean_rel_error'].min()
                max_error = alg_data['mean_rel_error'].max()
                print(f"    {algorithm}: MRE range {min_error:.2f} - {max_error:.2f}")
    
    print(f"\nAll plots completed! Created {len(unique_combinations)} plots.")

if __name__ == "__main__":
    run_synthetic_experiment()
    print("\nExperiment completed! Use plot_synthetic_results() to generate plots from stored data.")
