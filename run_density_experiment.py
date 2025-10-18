#!/usr/bin/env python3
"""
Density Experiment Runner
Tests the effect of graph density on biclique counting algorithms
"""

import os
import sqlite3
import subprocess
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def create_database():
    """Create database for density experiment results"""
    conn = sqlite3.connect('density-exp/density_experiment.db')
    cursor = conn.cursor()
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS density_experiment (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            dataset TEXT,
            algorithm TEXT,
            p INTEGER,
            q INTEGER,
            epsilon REAL,
            edge_density REAL,
            original_edges INTEGER,
            sampled_edges INTEGER,
            ground_truth REAL,
            estimate REAL,
            relative_error REAL,
            std_relative_error REAL,
            T_samples INTEGER,
            run_id INTEGER,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    conn.commit()
    conn.close()
    print("Database created successfully!")

def sample_edges(dataset, density):
    """Sample edges from dataset to create different density levels"""
    edge_file = f"/data/yizhangh/bidata/{dataset}.e"
    meta_file = f"/data/yizhangh/bidata/{dataset}.meta"
    
    if not os.path.exists(edge_file) or not os.path.exists(meta_file):
        return None, None
    
    # Read metadata
    with open(meta_file, 'r') as f:
        lines = f.readlines()
        upper_count = int(lines[0].strip())
        lower_count = int(lines[1].strip())
        total_edges = int(lines[2].strip())
    
    # Read all edges
    edges = []
    with open(edge_file, 'r') as f:
        for line in f:
            u, v = map(int, line.strip().split())
            edges.append((u, v))
    
    # Sample edges
    num_edges_to_sample = int(total_edges * density)
    sampled_edges = random.sample(edges, num_edges_to_sample)
    
    # Create temporary edge file
    temp_edge_file = f"density-exp/{dataset}_density_{density:.1f}.e"
    with open(temp_edge_file, 'w') as f:
        for u, v in sampled_edges:
            f.write(f"{u} {v}\n")
    
    # Create temporary meta file
    temp_meta_file = f"density-exp/{dataset}_density_{density:.1f}.meta"
    with open(temp_meta_file, 'w') as f:
        f.write(f"{upper_count}\n")
        f.write(f"{lower_count}\n")
        f.write(f"{num_edges_to_sample}\n")
    
    return temp_edge_file, temp_meta_file

def run_algorithm(dataset, algorithm, p, q, epsilon, density, run_id):
    """Run a single algorithm on sampled graph"""
    
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
    
    # Sample edges
    temp_edge_file, temp_meta_file = sample_edges(dataset, density)
    if not temp_edge_file:
        return {'success': False, 'error': 'Failed to sample edges'}
    
    try:
        # Run biclique algorithm
        biclique_path = '/data/yizhangh/ldp-pq/biclique'
        dataset_path = temp_meta_file.replace('.meta', '')  # Remove .meta extension
        cmd = [biclique_path, str(epsilon), dataset_path, '1', str(alg_num), str(p), str(q)]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            return {'success': False, 'error': f'Algorithm failed: {result.stderr}'}
        
        # Parse output
        output_lines = result.stdout.strip().split('\n')
        
        # Extract relative error from the last line
        for line in reversed(output_lines):
            if 'adv rel err' in line:
                rel_error = float(line.split('=')[-1].strip())
                break
        else:
            return {'success': False, 'error': 'Could not parse relative error'}
        
        # Get total edges from metadata
        with open(temp_meta_file, 'r') as f:
            lines = f.readlines()
            sampled_edges = int(lines[2].strip())
        
        return {
            'success': True,
            'relative_error': rel_error,
            'total_edges': sampled_edges,
            'sampled_edges': sampled_edges
        }
        
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Algorithm timed out'}
    except Exception as e:
        return {'success': False, 'error': str(e)}
    finally:
        # Clean up temporary files
        if os.path.exists(temp_edge_file):
            os.remove(temp_edge_file)
        if os.path.exists(temp_meta_file):
            os.remove(temp_meta_file)

def run_experiment():
    """Run the complete density experiment"""
    
    # Create output directory
    os.makedirs('density-exp', exist_ok=True)
    
    # Create database
    create_database()
    
    # Experiment parameters
    dataset = 'bpywiki'
    algorithms = ['Naive', 'oneR', 'ADV', 'ADV+', 'ADV++']
    p = 2
    q = 3
    epsilon = 1.0
    edge_densities = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    num_runs = 3  # Reduced for faster testing
    
    print(f"Starting density experiment on {dataset}")
    print(f"Algorithms: {algorithms}")
    print(f"Parameters: p={p}, q={q}, epsilon={epsilon}")
    print(f"Densities: {edge_densities}")
    print(f"Runs per configuration: {num_runs}")
    
    conn = sqlite3.connect('density-exp/density_experiment.db')
    cursor = conn.cursor()
    
    total_experiments = len(algorithms) * len(edge_densities) * num_runs
    experiment_count = 0
    
    for algorithm in algorithms:
        print(f"\nRunning {algorithm}...")
        
        for density in edge_densities:
            print(f"  Density: {density:.1f}")
            rel_errors = []
            
            for run_id in range(num_runs):
                experiment_count += 1
                print(f"    Run {run_id+1}/{num_runs} ({experiment_count}/{total_experiments})")
                
                result = run_algorithm(dataset, algorithm, p, q, epsilon, density, run_id)
                
                if result['success']:
                    rel_errors.append(result['relative_error'])
                    
                    # Store result in database
                    cursor.execute('''
                        INSERT INTO density_experiment 
                        (dataset, algorithm, p, q, epsilon, edge_density, 
                         original_edges, sampled_edges, ground_truth, estimate, 
                         relative_error, std_relative_error, T_samples, run_id)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (
                        dataset, algorithm, p, q, epsilon, density,
                        result['total_edges'], result['sampled_edges'], 
                        0.0, 0.0, result['relative_error'], 0.0, -1, run_id
                    ))
                else:
                    print(f"      Failed: {result.get('error', 'Unknown error')}")
            
            if rel_errors:
                mean_rel_error = np.mean(rel_errors)
                std_rel_error = np.std(rel_errors)
                print(f"    Mean MRE: {mean_rel_error:.3f} ± {std_rel_error:.3f}")
    
    conn.commit()
    conn.close()
    
    print(f"\nExperiment completed! Results stored in density-exp/density_experiment.db")

def plot_results():
    """Plot the density experiment results"""
    
    # Load results from database
    conn = sqlite3.connect('density-exp/density_experiment.db')
    
    query = '''
    SELECT dataset, algorithm, p, q, edge_density, 
           AVG(relative_error) as mean_rel_error,
           COUNT(*) as num_runs
    FROM density_experiment 
    GROUP BY dataset, algorithm, p, q, edge_density
    ORDER BY dataset, algorithm, p, q, edge_density
    '''
    
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    if df.empty:
        print("No results found in database!")
        return
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Algorithm colors and styles
    algorithm_styles = {
        'Naive': {'color': '#1f77b4', 'marker': 'o', 'linestyle': '-'},
        'oneR': {'color': '#ff7f0e', 'marker': 's', 'linestyle': '-'},
        'ADV': {'color': '#2ca02c', 'marker': '^', 'linestyle': '-'},
        'ADV+': {'color': '#d62728', 'marker': 'd', 'linestyle': '-'},
        'ADV++': {'color': '#9467bd', 'marker': 'v', 'linestyle': '-'}
    }
    
    dataset = df['dataset'].iloc[0]
    p = df['p'].iloc[0]
    q = df['q'].iloc[0]
    
    for algorithm in df['algorithm'].unique():
        alg_data = df[df['algorithm'] == algorithm]
        
        style = algorithm_styles.get(algorithm, {'color': 'gray', 'marker': 'o', 'linestyle': '-'})
        
        plt.plot(
            alg_data['edge_density'] * 100,  # Convert to percentage
            alg_data['mean_rel_error'],
            label=algorithm,
            **style,
            linewidth=2,
            markersize=8
        )
    
    plt.xlabel('Edge Density (%)', fontsize=14)
    plt.ylabel('Mean Relative Error', fontsize=14)
    plt.title(f'Effect of Graph Density on Algorithm Performance\\n{dataset} (P={p}, Q={q})', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    # Set x-axis to show percentages (10%, 20%, ..., 90%)
    plt.xticks([10, 20, 30, 40, 50, 60, 70, 80, 90])
    
    plt.tight_layout()
    plt.savefig(f'density-exp/{dataset}_P{p}_Q{q}_density_vs_mre.pdf', dpi=300)
    plt.show()
    
    print(f"Plot saved as density-exp/{dataset}_P{p}_Q{q}_density_vs_mre.pdf")

if __name__ == "__main__":
    run_experiment()
    plot_results()
