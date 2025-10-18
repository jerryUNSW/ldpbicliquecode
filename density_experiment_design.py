#!/usr/bin/env python3
"""
Experiment Design: Effect of Graph Density on Algorithm Performance

This script designs an experiment to evaluate how graph density affects the 
performance of naive and advanced biclique counting algorithms.

Experiment Design:
- Keep vertices fixed (n1, n2)
- Sample edges at different densities: 10%, 20%, ..., 90%
- Measure Mean Relative Error (MRE) for each algorithm
- Plot MRE vs. edge density
"""

import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import os

def create_density_experiment_table():
    """Create a new table for density experiment results"""
    
    conn = sqlite3.connect('density_experiment.db')
    cursor = conn.cursor()
    
    # Create table for density experiment results
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS density_experiment (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            dataset TEXT NOT NULL,
            algorithm TEXT NOT NULL,
            p INTEGER NOT NULL,
            q INTEGER NOT NULL,
            epsilon REAL NOT NULL,
            edge_density REAL NOT NULL,  -- 0.1 to 0.9 (10% to 90%)
            original_edges INTEGER NOT NULL,
            sampled_edges INTEGER NOT NULL,
            ground_truth REAL NOT NULL,
            estimate REAL NOT NULL,
            relative_error REAL NOT NULL,
            std_relative_error REAL NOT NULL,
            T_samples INTEGER DEFAULT -1,
            UNIQUE(dataset, algorithm, p, q, epsilon, edge_density)
        )
    ''')
    
    conn.commit()
    conn.close()
    print("Created density_experiment table")

def design_experiment_parameters():
    """Define experiment parameters"""
    
    experiment_config = {
        # Datasets to test (start with small dataset)
        'datasets': [
            'bpywiki',      # Small dataset to start with
        ],
        
        # Algorithms to compare (for p=2)
        'algorithms': [
            'Naive',
            'oneR',
            'ADV',      # MRCN
            'ADV+',     # MRCN+
            'ADV++'     # MRCN++
        ],
        
        # Graph parameters
        'p_values': [2],  # Focus on p=2
        'q_values': [3],  # Focus on q=3
        'epsilon': 1.0,
        
        # Density levels to test
        'edge_densities': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        
        # Number of runs per configuration (for statistical significance)
        'num_runs': 5,
        
        # Output settings
        'output_dir': 'density-exp',
        'plot_format': 'pdf'
    }
    
    return experiment_config

def generate_experiment_script():
    """Generate the actual experiment script"""
    
    script_content = '''#!/usr/bin/env python3
"""
Density Experiment Runner
This script runs the density experiment by sampling edges at different densities
and measuring algorithm performance.
"""

import os
import sys
import subprocess
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def sample_edges(input_file, output_file, density):
    """Sample edges from input file with given density"""
    
    # Read all edges
    with open(input_file, 'r') as f:
        edges = f.readlines()
    
    total_edges = len(edges)
    num_sample = int(total_edges * density)
    
    # Randomly sample edges
    sampled_edges = np.random.choice(edges, size=num_sample, replace=False)
    
    # Write sampled edges
    with open(output_file, 'w') as f:
        f.writelines(sampled_edges)
    
    return total_edges, num_sample

def run_algorithm_on_sampled_graph(dataset, algorithm, p, q, epsilon, density, run_id):
    """Run a single algorithm on a sampled graph"""
    
    # Create sampled graph file
    original_edge_file = f"../bidata/{dataset}.e"
    sampled_edge_file = f"temp_{dataset}_{density}_{run_id}.e"
    
    # Sample edges
    total_edges, sampled_edges = sample_edges(original_edge_file, sampled_edge_file, density)
    
    # Create temporary meta file
    meta_file = f"temp_{dataset}_{density}_{run_id}.meta"
    with open(f"../bidata/{dataset}.meta", 'r') as f:
        meta_lines = f.readlines()
    
    # Update edge count in meta file
    with open(meta_file, 'w') as f:
        f.write(meta_lines[0])  # n1
        f.write(meta_lines[1])  # n2
        f.write(f"{sampled_edges}\\n")  # Updated edge count
    
    # Map algorithm names to algorithm numbers
    alg_map = {
        'Naive': 0,
        'oneR': 1,
        'ADV': 2,
        'ADV+': 3,
        'ADV++': 4
    }
    
    alg_num = alg_map[algorithm]
    
    # Run the algorithm
    cmd = [
        "./biclique",
        str(alg_num),
        f"temp_{dataset}_{density}_{run_id}",
        str(epsilon),
        "0",  # sampling ratio
        str(p),
        str(q)
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            # Parse output for relative error
            output_lines = result.stdout.split('\\n')
            rel_error = None
            
            for line in output_lines:
                if 'adv rel err' in line:
                    rel_error = float(line.split('=')[1].strip())
                    break
            
            if rel_error is not None:
                return {
                    'relative_error': rel_error,
                    'total_edges': total_edges,
                    'sampled_edges': sampled_edges,
                    'success': True
                }
        
        return {'success': False, 'error': result.stderr}
        
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Timeout'}
    finally:
        # Clean up temporary files
        for temp_file in [sampled_edge_file, meta_file]:
            if os.path.exists(temp_file):
                os.remove(temp_file)

def run_density_experiment():
    """Run the complete density experiment"""
    
    config = {
        'datasets': ['bpywiki'],
        'algorithms': ['Naive', 'oneR', 'ADV', 'ADV+', 'ADV++'],
        'p_values': [2],
        'q_values': [3],
        'edge_densities': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        'num_runs': 3,
        'epsilon': 1.0
    }
    
    # Create results directory
    os.makedirs('density_experiment_results', exist_ok=True)
    
    # Connect to database
    conn = sqlite3.connect('density_experiment.db')
    cursor = conn.cursor()
    
    total_experiments = (len(config['datasets']) * 
                        len(config['algorithms']) * 
                        len(config['p_values']) * 
                        len(config['q_values']) * 
                        len(config['edge_densities']) * 
                        config['num_runs'])
    
    experiment_count = 0
    
    print(f"Starting density experiment with {total_experiments} total runs...")
    
    for dataset in config['datasets']:
        print(f"\\nProcessing dataset: {dataset}")
        
        for algorithm in config['algorithms']:
            print(f"  Algorithm: {algorithm}")
            
            for p in config['p_values']:
                for q in config['q_values']:
                    print(f"    P={p}, Q={q}")
                    
                    for density in config['edge_densities']:
                        rel_errors = []
                        
                        for run_id in range(config['num_runs']):
                            experiment_count += 1
                            print(f"      Density={density:.1f}, Run={run_id+1} ({experiment_count}/{total_experiments})")
                            
                            result = run_algorithm_on_sampled_graph(
                                dataset, algorithm, p, q, config['epsilon'], density, run_id
                            )
                            
                            if result['success']:
                                rel_errors.append(result['relative_error'])
                                
                                # Store individual result
                                cursor.execute('''
INSERT OR REPLACE INTO density_experiment 
(dataset, algorithm, p, q, epsilon, edge_density, 
 original_edges, sampled_edges, ground_truth, estimate, 
 relative_error, std_relative_error, T_samples)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                                ''', (
                                    dataset, algorithm, p, q, config['epsilon'], density,
                                    result['total_edges'], result['sampled_edges'], 
                                    0.0, 0.0, result['relative_error'], 0.0, -1
                                ))
                            else:
                                print(f"        Failed: {result.get('error', 'Unknown error')}")
                        
                        if rel_errors:
                            mean_rel_error = np.mean(rel_errors)
                            std_rel_error = np.std(rel_errors)
                            print(f"      Mean MRE: {mean_rel_error:.3f} ± {std_rel_error:.3f}")
    
    conn.commit()
    conn.close()
    print("\\nDensity experiment completed!")

if __name__ == "__main__":
    run_density_experiment()
'''
    
    with open('run_density_experiment.py', 'w') as f:
        f.write(script_content)
    
    print("Generated run_density_experiment.py")

def generate_analysis_script():
    """Generate script to analyze and plot results"""
    
    analysis_script = '''#!/usr/bin/env python3
"""
Density Experiment Analysis
Analyze and plot the results of the density experiment.
"""

import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def load_density_results():
    """Load results from density experiment database"""
    
    conn = sqlite3.connect('density_experiment.db')
    
    query = '''
    SELECT dataset, algorithm, p, q, edge_density, 
           AVG(relative_error) as mean_rel_error,
           STDDEV(relative_error) as std_rel_error,
           COUNT(*) as num_runs
    FROM density_experiment 
    GROUP BY dataset, algorithm, p, q, edge_density
    ORDER BY dataset, algorithm, p, q, edge_density
    '''
    
    df = pd.read_sql_query(query, conn)
    conn.close()
    
    return df

def plot_density_vs_mre(df, dataset, p=2, q=3):
    """Plot MRE vs edge density for a specific dataset and Q value"""
    
    # Filter data
    data = df[(df['dataset'] == dataset) & (df['p'] == p) & (df['q'] == q)]
    
    if data.empty:
        print(f"No data found for {dataset}, P={p}, Q={q}")
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
    
    for algorithm in data['algorithm'].unique():
        alg_data = data[data['algorithm'] == algorithm]
        
        style = algorithm_styles.get(algorithm, {'color': 'gray', 'marker': 'o', 'linestyle': '-'})
        
        plt.errorbar(
            alg_data['edge_density'] * 100,  # Convert to percentage
            alg_data['mean_rel_error'],
            yerr=alg_data['std_rel_error'],
            label=algorithm,
            **style,
            linewidth=2,
            markersize=8,
            capsize=5
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

def plot_all_density_results():
    """Plot results for all datasets and Q values"""
    
    df = load_density_results()
    
    if df.empty:
        print("No data found in density experiment database")
        return
    
    # Get unique combinations
    datasets = df['dataset'].unique()
    q_values = df['q'].unique()
    
    for dataset in datasets:
        for q in q_values:
            plot_density_vs_mre(df, dataset, p=2, q=q)

def generate_summary_table():
    """Generate a summary table of results"""
    
    df = load_density_results()
    
    if df.empty:
        print("No data found")
        return
    
    # Create pivot table
    summary = df.pivot_table(
        index=['dataset', 'algorithm', 'p', 'q'],
        columns='edge_density',
        values='mean_rel_error',
        aggfunc='mean'
    )
    
    print("\\nDensity Experiment Summary (Mean Relative Error):")
    print("=" * 60)
    print(summary.round(3))
    
    # Save to CSV
    summary.to_csv('density_experiment_results/summary_table.csv')
    print("\\nSummary saved to density_experiment_results/summary_table.csv")

if __name__ == "__main__":
    plot_all_density_results()
    generate_summary_table()
'''
    
    with open('analyze_density_experiment.py', 'w') as f:
        f.write(analysis_script)
    
    print("Generated analyze_density_experiment.py")

def main():
    """Main function to set up the density experiment"""
    
    print("Setting up Density Experiment...")
    print("=" * 50)
    
    # Create database table
    create_density_experiment_table()
    
    # Generate experiment configuration
    config = design_experiment_parameters()
    print(f"\\nExperiment Configuration:")
    print(f"  Datasets: {config['datasets']}")
    print(f"  Algorithms: {config['algorithms']}")
    print(f"  Edge Densities: {config['edge_densities']}")
    print(f"  Q Values: {config['q_values']}")
    print(f"  Number of runs per config: {config['num_runs']}")
    
    # Generate scripts
    generate_experiment_script()
    generate_analysis_script()
    
    # Create output directory
    os.makedirs('density-exp', exist_ok=True)
    
    print(f"\\nExperiment setup complete!")
    print(f"\\nTo run the experiment:")
    print(f"  1. python3 run_density_experiment.py")
    print(f"  2. python3 analyze_density_experiment.py")
    
    print(f"\\nExpected outputs:")
    print(f"  - density_experiment.db (results database)")
    print(f"  - density_experiment_results/ (plots and summary)")
    print(f"  - Plots showing MRE vs edge density for each algorithm")

if __name__ == "__main__":
    main()
