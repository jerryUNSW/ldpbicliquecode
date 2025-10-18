#!/usr/bin/env python3
"""
Separate plotting script for synthetic density experiment results
This script reads from the database and creates plots for different P,Q combinations
Uses the professional style from budget allocation plots (Figure 9 style)
"""

import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import os

# Set publication-quality plot parameters (matching budget allocation plots)
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams["legend.framealpha"] = 0
plt.rcParams["legend.handletextpad"] = 0.1
plt.rcParams["legend.columnspacing"] = 0.2
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 8

def plot_synthetic_results(db_path='density-exp/synthetic_density_experiment.db'):
    """Plot the synthetic density experiment results for each P,Q combination"""
    
    # Load results from database
    conn = sqlite3.connect(db_path)
    
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
    
    # Algorithm mapping from database names to display names
    algorithm_display_names = {
        'Naive': 'Naive',
        'oneR': 'oneR', 
        'ADV++': 'MRCN++'
    }
    
    # Algorithm colors and styles (black, red, blue with distinct markers)
    algorithm_styles = {
        'Naive': {'color': 'black', 'marker': 'o', 'linestyle': '-'},
        'oneR': {'color': 'red', 'marker': 's', 'linestyle': '-'},
        'MRCN++': {'color': 'blue', 'marker': '^', 'linestyle': '-'}
    }
    
    # Only plot the algorithms we want (Naive, oneR, ADV++ from database)
    algorithms_to_plot = ['Naive', 'oneR', 'ADV++']
    
    # Get unique combinations of graph sizes and P,Q values
    unique_combinations = df[['n1', 'n2', 'p', 'q']].drop_duplicates().values
    
    print(f"Creating plots for {len(unique_combinations)} different combinations...")
    
    # Create a separate plot for each combination
    for n1, n2, p, q in unique_combinations:
        print(f"\nCreating plot for graph size {n1}×{n2}, P={p}, Q={q}...")
        
        # Filter data for this combination
        size_data = df[(df['n1'] == n1) & (df['n2'] == n2) & (df['p'] == p) & (df['q'] == q)]
        
        # Filter to only include densities from 10% to 100% (0.1 to 1.0) with 10% increments
        density_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        size_data = size_data[size_data['edge_density'].isin(density_values)]
        
        if size_data.empty:
            print(f"  No data found for graph size {n1}×{n2}, P={p}, Q={q}")
            continue
        
        # Create plot
        plt.figure(figsize=(6, 5))  # Match budget allocation plot size
        
        for algorithm in algorithms_to_plot:
            if algorithm in size_data['algorithm'].values:
                alg_data = size_data[size_data['algorithm'] == algorithm]
                display_name = algorithm_display_names.get(algorithm, algorithm)
                style = algorithm_styles.get(display_name, {'color': 'black', 'marker': 'o', 'linestyle': '-'})
                
                # Sort data by edge_density to ensure proper line connections
                alg_data = alg_data.sort_values('edge_density')
                
                # Filter out infinite values for better visualization
                alg_data_clean = alg_data[alg_data['mean_rel_error'] != float('inf')]
                
                if len(alg_data_clean) > 0:
                    # Plot with professional style (matching budget allocation plots)
                    plt.plot(
                        alg_data_clean['edge_density'] * 100,  # Convert to percentage (10, 20, etc.)
                        alg_data_clean['mean_rel_error'],
                        marker=style['marker'], 
                        markerfacecolor='white', 
                        markeredgecolor=style['color'],
                        markeredgewidth=1.5, 
                        color=style['color'], 
                        linewidth=2,
                        linestyle=style['linestyle'],
                        markersize=8,  # Match biclique plot marker size
                        label=display_name
                    )
        
        # Labels (matching budget allocation plot style)
        plt.xlabel('Density (%)', fontsize=20)
        plt.ylabel('mean relative error', fontsize=20)
        
        # No title for density plots
        
        # Legend
        plt.legend(loc='best', fontsize=14, ncol=1)
        
        # Grid
        plt.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
        plt.yscale('log')
        
        # Format x-axis to show density values nicely (every 20% from 20% to 100%)
        plt.xticks([20, 40, 60, 80, 100], 
                   ['20', '40', '60', '80', '100'], 
                   rotation=0, fontsize=20)
        plt.tick_params(axis='y', labelsize=20)
        
        plt.tight_layout()
        plt.savefig(f'density-exp/synthetic_{n1}x{n2}_P{p}_Q{q}_extreme_density.pdf', dpi=300, bbox_inches='tight')
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

def plot_specific_combinations(db_path='density-exp/synthetic_density_experiment.db', 
                              graph_sizes=None, p_values=None, q_values=None):
    """Plot specific combinations of graph sizes and P,Q values"""
    
    if graph_sizes is None:
        graph_sizes = [(50, 100), (100, 200), (150, 300), (200, 400)]
    if p_values is None:
        p_values = [2, 3]
    if q_values is None:
        q_values = [3, 4, 5]
    
    # Load results from database
    conn = sqlite3.connect(db_path)
    
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
    
    # Algorithm colors and styles (black, red, blue with distinct markers)
    algorithm_styles = {
        'Naive': {'color': 'black', 'marker': 'o', 'linestyle': '-'},
        'oneR': {'color': 'red', 'marker': 's', 'linestyle': '-'},
        'MRCN++': {'color': 'blue', 'marker': '^', 'linestyle': '-'}
    }
    
    algorithms_to_plot = ['Naive', 'oneR', 'MCN++']
    
    # Filter for specific combinations
    combinations_to_plot = []
    for n1, n2 in graph_sizes:
        for p in p_values:
            for q in q_values:
                if ((df['n1'] == n1) & (df['n2'] == n2) & (df['p'] == p) & (df['q'] == q)).any():
                    combinations_to_plot.append((n1, n2, p, q))
    
    print(f"Creating plots for {len(combinations_to_plot)} specified combinations...")
    
    for n1, n2, p, q in combinations_to_plot:
        print(f"\nCreating plot for graph size {n1}×{n2}, P={p}, Q={q}...")
        
        # Filter data for this combination
        size_data = df[(df['n1'] == n1) & (df['n2'] == n2) & (df['p'] == p) & (df['q'] == q)]
        
        if size_data.empty:
            print(f"  No data found for graph size {n1}×{n2}, P={p}, Q={q}")
            continue
        
        # Create plot
        plt.figure(figsize=(6, 5))  # Match budget allocation plot size
        
        for algorithm in algorithms_to_plot:
            if algorithm in size_data['algorithm'].values:
                alg_data = size_data[size_data['algorithm'] == algorithm]
                display_name = algorithm_display_names.get(algorithm, algorithm)
                style = algorithm_styles.get(display_name, {'color': 'black', 'marker': 'o', 'linestyle': '-'})
                
                # Sort data by edge_density to ensure proper line connections
                alg_data = alg_data.sort_values('edge_density')
                
                # Filter out infinite values for better visualization
                alg_data_clean = alg_data[alg_data['mean_rel_error'] != float('inf')]
                
                if len(alg_data_clean) > 0:
                    # Plot with professional style (matching budget allocation plots)
                    plt.plot(
                        alg_data_clean['edge_density'] * 100,  # Convert to percentage (10, 20, etc.)
                        alg_data_clean['mean_rel_error'],
                        marker=style['marker'], 
                        markerfacecolor='white', 
                        markeredgecolor=style['color'],
                        markeredgewidth=1.5, 
                        color=style['color'], 
                        linewidth=2,
                        linestyle=style['linestyle'],
                        markersize=8,  # Match biclique plot marker size
                        label=display_name
                    )
        
        plt.xlabel('Density (%)', fontsize=24)
        plt.ylabel('Mean Relative Error', fontsize=24)
        plt.title(f'Synthetic Graph, n1={n1}, n2={n2}, P={p}, Q={q}', fontsize=28)
        plt.legend(fontsize=20)
        plt.grid(True, alpha=0.3)
        plt.yscale('log')
        
        # Set x-axis to show percentages (every 10%)
        plt.xticks([10, 20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=18)
        plt.yticks(fontsize=18)
        
        plt.tight_layout()
        plt.savefig(f'density-exp/synthetic_{n1}x{n2}_P{p}_Q{q}_extreme_density.pdf', dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"  Plot saved as density-exp/synthetic_{n1}x{n2}_P{p}_Q{q}_extreme_density.pdf")
    
    print(f"\nAll plots completed! Created {len(combinations_to_plot)} plots.")

if __name__ == "__main__":
    # Create output directory if it doesn't exist
    os.makedirs('density-exp', exist_ok=True)
    
    # Plot all available combinations
    plot_synthetic_results()
    
    # Example: Plot only specific combinations
    # plot_specific_combinations(
    #     graph_sizes=[(100, 200), (200, 400)],
    #     p_values=[2],
    #     q_values=[3, 4]
    # )
