#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sqlite3
import os
import math

def plt_settings():
    """Set matplotlib settings to match budget allocation styled plotting"""
    # Set publication-quality plot parameters (matching budget allocation script)
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams["legend.framealpha"] = 0
    plt.rcParams["legend.handletextpad"] = 0.1
    plt.rcParams["legend.columnspacing"] = 0.2
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 8
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['xtick.minor.width'] = 1.0
    plt.rcParams['ytick.minor.width'] = 1.0
    plt.rcParams["figure.figsize"] = (6, 5)

def get_pos_and_labels(indices):
    """Get positions and labels for y-axis ticks"""
    pos = [10**i for i in indices]
    labels = [f'$10^{{{i}}}$' for i in indices]
    return pos, labels

def get_data_from_db(db_path, dataset, epsilon=1.0):
    """Get algorithm comparison data from SQLite database"""
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Query to get all results for the dataset
    query = """
    SELECT algorithm, q, relative_error, T_samples
    FROM algorithm_comparison_p3 
    WHERE dataset = ? AND epsilon = ? AND p = 3
    ORDER BY algorithm, q
    """
    
    cursor.execute(query, (dataset, epsilon))
    rows = cursor.fetchall()
    conn.close()
    
    # Organize data into results dictionary
    results = {}
    algorithms = set()
    
    for row in rows:
        algorithm_name, q_value, relative_error, t_samples = row
        algorithms.add(algorithm_name)
        
        if q_value not in results:
            results[q_value] = {}
        
        results[q_value][algorithm_name] = relative_error
    
    # Sort algorithms with Naive first, then others alphabetically
    # Map ADV names to MRCN names for display
    algorithm_mapping = {
        'ADV': 'MRCN',
        'ADV+': 'MRCN+', 
        'ADV++': 'MRCN++'
    }
    
    # Apply mapping to algorithm names
    mapped_algorithms = [algorithm_mapping.get(alg, alg) for alg in algorithms]
    algorithm_order = ['Naive'] + sorted([alg for alg in mapped_algorithms if alg != 'Naive'])
    
    return results, algorithm_order

def create_p3_plot(dataset, epsilon, results, algorithms):
    """Create the algorithm comparison plot with proper styling for p=3"""
    
    plt_settings()
    
    # Prepare data
    q_values = sorted(results.keys())
    
    # Define colors and styles matching the p=2 plots
    method_styles = {
        'Naive': {'color': 'white', 'hatch': '', 'edgecolor': 'black'},
        'OneR': {'color': 'lightgray', 'hatch': '', 'edgecolor': 'black'},
        'MRCN': {'color': 'grey', 'hatch': '', 'edgecolor': 'black'},
        'MRCN+': {'color': 'black', 'hatch': '////', 'edgecolor': 'black'},
        'MRCN++': {'color': 'darkgrey', 'hatch': '//', 'edgecolor': 'black'}
    }
    
    # Create figure with more bottom margin
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.15)  # Add more bottom margin
    
    # Plot bars for each algorithm
    x = np.arange(len(q_values))
    total_width, n = 0.8, len(algorithms)
    width = total_width / n
    
    # Create mapping from display names back to database names
    reverse_mapping = {
        'MRCN': 'ADV',
        'MRCN+': 'ADV+', 
        'MRCN++': 'ADV++'
    }
    
    for i, alg in enumerate(algorithms):
        values = []
        # Get the database name for data lookup
        db_alg = reverse_mapping.get(alg, alg)
        
        for j, q in enumerate(q_values):
            val = results[q].get(db_alg, float('inf'))
            # Handle infinite values by setting them to a large number for plotting
            val = val if val != float('inf') else 1e6
            values.append(val)
        
        if values:  # Only plot if there are valid values
            style = method_styles.get(alg, {'color': 'lightgrey', 'hatch': '', 'edgecolor': 'black'})
            ax.bar(x + i * width, values, width=width, 
                   label=alg, 
                   linewidth=0.5, 
                   edgecolor=style['edgecolor'],
                   fc=style['color'],
                   hatch=style['hatch'])
    
    # Customize plot
    ax.set_xlabel('q', fontsize=20)
    ax.set_ylabel('Mean Relative Error', fontsize=20)
    ax.set_yscale('log')
    
    # Set x-axis ticks - center them under the bars
    ax.set_xticks(x + (n-1) * width / 2)
    ax.set_xticklabels([str(q) for q in q_values], fontsize=20)
    
    # Set y-axis limits and ticks - always use log scale
    all_values = []
    # Create mapping from display names back to database names
    reverse_mapping = {
        'MRCN': 'ADV',
        'MRCN+': 'ADV+', 
        'MRCN++': 'ADV++'
    }
    
    for q in q_values:
        for alg in algorithms:
            # Get the database name for data lookup
            db_alg = reverse_mapping.get(alg, alg)
            val = results[q].get(db_alg, float('inf'))
            if val != float('inf'):
                all_values.append(val)
    
    if all_values:
        y_min = np.min(all_values)
        y_max = np.max(all_values)
        
        # Print debug info
        print(f"Y-axis range: min={y_min:.3f}, max={y_max:.3f}")
        print(f"All values: {sorted(all_values)}")
        
        # Use a more reasonable range that shows all algorithms clearly
        # For log scale, use powers of 10 that encompass the data
        y_min_log = math.floor(math.log10(y_min)) if y_min > 0 else -3
        y_max_log = math.ceil(math.log10(y_max)) if y_max > 0 else 3
        
        # Standardize y-axis limits for jester plots to ensure consistent visual size
        # Use fixed limits that work for all jester plots (p3 and p4, eps1 and eps2)
        y_min_adjusted = 0.001  # Fixed at 10^-3 for all jester plots
        y_max_adjusted = 10000.0  # Fixed at 10^4 for all jester plots
        
        ax.set_ylim(y_min_adjusted, y_max_adjusted)
        
        # Set y-axis ticks - always use log scale
        # Create more appropriate tick spacing based on the range
        min_tick = math.floor(math.log10(y_min_adjusted))
        max_tick = math.ceil(math.log10(y_max_adjusted))
        
        # For small ranges, use every power of 10
        if max_tick - min_tick <= 4:
            step = 1
        else:
            step = 2
            
        indices = [i for i in range(min_tick, max_tick + 1, step)]
        pos, labels = get_pos_and_labels(indices)
        ax.set_yticks(pos)
        ax.set_yticklabels(labels, fontsize=20)
    
    # Legend
    ax.legend(fontsize=14, ncol=2, loc="upper left", columnspacing=0.5, frameon=False)
    
    plt.tight_layout()
    return fig

def main():
    # Database path
    db_path = "algorithm_comparison_p3.db"
    
    if not os.path.exists(db_path):
        print(f"Database file {db_path} not found.")
        return
    
    # List of datasets to plot (modify this list as needed)
    datasets_to_plot = ["lrcwiki", "edit-stwiktionary", "nips", "csbwiki", "librec-filmtrust-ratings", "jester2-swapped"]
    
    # Epsilon values to plot
    epsilon_values = [1.0, 2.0]
    
    print(f"Plotting datasets: {datasets_to_plot}")
    
    for dataset in datasets_to_plot:
        for epsilon in epsilon_values:
            print(f"\nProcessing {dataset}, ε={epsilon}...")
            
            try:
                # Get data from database
                results, algorithms = get_data_from_db(db_path, dataset, epsilon=epsilon)
                
                if not results:
                    print(f"No results found for {dataset}, ε={epsilon}")
                    continue
                
                print(f"Found results for {dataset} dataset, ε={epsilon}")
                print(f"Q values: {sorted(results.keys())}")
                print(f"Algorithms: {algorithms}")
                
                # Print summary of results
                for q in sorted(results.keys()):
                    print(f"  Q={q}:")
                    for alg in algorithms:
                        if alg in results[q]:
                            print(f"    {alg}: {results[q][alg]:.3f}")
                
                # Create and save main plot
                fig = create_p3_plot(dataset, epsilon, results, algorithms)
                eps_str = str(int(epsilon)) if epsilon == int(epsilon) else str(epsilon).replace('.', '')
                output_filename = f'algorithm_comparison_{dataset}_eps{eps_str}_FIXED_p3.pdf'
                fig.savefig(output_filename, dpi=300, bbox_inches='tight')
                print(f"Plot saved as: {output_filename}")
                plt.close()
                
                print(f"Completed plots for {dataset}, ε={epsilon}")
                print("---")
                
            except Exception as e:
                print(f"Error processing {dataset}, ε={epsilon}: {e}")
                continue

if __name__ == "__main__":
    main()
