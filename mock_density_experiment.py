#!/usr/bin/env python3
"""
Mock Density Experiment
Creates sample data to demonstrate the density experiment plotting functionality
"""

import os
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def create_mock_database():
    """Create database with mock density experiment results"""
    
    # Create output directory
    os.makedirs('density-exp', exist_ok=True)
    
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
    
    # Clear existing data
    cursor.execute('DELETE FROM density_experiment')
    
    # Generate mock data
    dataset = 'bpywiki'
    algorithms = ['Naive', 'oneR', 'ADV', 'ADV+', 'ADV++']
    p = 2
    q = 3
    epsilon = 1.0
    edge_densities = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    num_runs = 3
    
    # Mock relative errors that show realistic trends
    # Naive: high error, decreases with density
    # oneR: moderate error, decreases with density  
    # ADV: low error, relatively stable
    # ADV+: very low error, relatively stable
    # ADV++: lowest error, relatively stable
    
    for algorithm in algorithms:
        for density in edge_densities:
            for run_id in range(num_runs):
                # Generate realistic relative errors with some noise
                if algorithm == 'Naive':
                    base_error = 100.0 * (1.0 - density * 0.7)  # High error, decreases with density
                    noise = np.random.normal(0, base_error * 0.1)
                elif algorithm == 'oneR':
                    base_error = 20.0 * (1.0 - density * 0.5)  # Moderate error, decreases with density
                    noise = np.random.normal(0, base_error * 0.1)
                elif algorithm == 'ADV':
                    base_error = 5.0 + np.random.normal(0, 1.0)  # Low error, relatively stable
                    noise = np.random.normal(0, 0.5)
                elif algorithm == 'ADV+':
                    base_error = 2.0 + np.random.normal(0, 0.5)  # Very low error, relatively stable
                    noise = np.random.normal(0, 0.2)
                else:  # ADV++
                    base_error = 1.0 + np.random.normal(0, 0.2)  # Lowest error, relatively stable
                    noise = np.random.normal(0, 0.1)
                
                relative_error = max(0.01, base_error + noise)  # Ensure positive
                
                cursor.execute('''
                    INSERT INTO density_experiment 
                    (dataset, algorithm, p, q, epsilon, edge_density, 
                     original_edges, sampled_edges, ground_truth, estimate, 
                     relative_error, std_relative_error, T_samples, run_id)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    dataset, algorithm, p, q, epsilon, density,
                    100000, int(100000 * density), 
                    0.0, 0.0, relative_error, 0.0, -1, run_id
                ))
    
    conn.commit()
    conn.close()
    print("Mock database created with sample data!")

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
    
    # Print summary statistics
    print("\nSummary of results:")
    for algorithm in df['algorithm'].unique():
        alg_data = df[df['algorithm'] == algorithm]
        min_error = alg_data['mean_rel_error'].min()
        max_error = alg_data['mean_rel_error'].max()
        print(f"  {algorithm}: MRE range {min_error:.2f} - {max_error:.2f}")

if __name__ == "__main__":
    create_mock_database()
    plot_results()

