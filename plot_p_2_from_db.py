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
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    plt.rcParams['font.size'] = 18
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['xtick.minor.width'] = 1.0
    plt.rcParams['ytick.minor.width'] = 1.0

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
    WHERE dataset = ? AND epsilon = ? AND p = 2
    ORDER BY algorithm, q
    """
    
    cursor.execute(query, (dataset, epsilon))
    rows = cursor.fetchall()
    conn.close()
    
    if not rows:
        return None, None, None
    
    # Organize data by algorithm
    algorithms = {}
    q_values = set()
    
    for algorithm, q, relative_error, T_samples in rows:
        if algorithm not in algorithms:
            algorithms[algorithm] = {}
        algorithms[algorithm][q] = relative_error
        q_values.add(q)
    
    q_values = sorted(list(q_values))
    
    return algorithms, q_values, None

def create_p2_plot(dataset, algorithms, q_values, output_file):
    """Create the algorithm comparison bar plot for p=2"""
    
    plt_settings()
    
    # Prepare data
    q_values = sorted(q_values)
    
    # Define colors and styles matching the p=3 plots
    method_styles = {
        'Naive': {'color': 'white', 'hatch': '', 'edgecolor': 'black'},
        'OneR': {'color': 'lightgrey', 'hatch': '', 'edgecolor': 'black'},
        'MRCN': {'color': 'grey', 'hatch': '', 'edgecolor': 'black'},
        'MRCN+': {'color': 'black', 'hatch': '////', 'edgecolor': 'black'},
        'MRCN++': {'color': 'darkgrey', 'hatch': '//', 'edgecolor': 'black'}
    }
    
    # Create figure with more bottom margin
    fig, ax = plt.subplots(figsize=(10, 7))
    plt.subplots_adjust(bottom=0.15)  # Add more bottom margin
    
    # Plot bars for each algorithm
    x = np.arange(len(q_values))
    total_width, n = 0.8, len(algorithms)
    width = total_width / n
    
    for i, alg in enumerate(algorithms):
        values = []
        
        for j, q in enumerate(q_values):
            val = algorithms[alg].get(q, float('inf'))
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
    ax.set_xlabel('q', fontsize=32)
    ax.set_ylabel('Mean Relative Error', fontsize=32)
    ax.set_title('(P = 2, ε = 1)', fontsize=30)
    ax.set_yscale('log')
    
    # Set x-axis ticks - center them under the bars
    ax.set_xticks(x + (n-1) * width / 2)
    ax.set_xticklabels([str(q) for q in q_values], fontsize=26)
    
    # Set y-axis limits and ticks - always use log scale
    all_values = []
    for alg in algorithms:
        for q in q_values:
            val = algorithms[alg].get(q, float('inf'))
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
        
        # Ensure we have enough range to show all values clearly
        # For small values, don't go too far below
        y_min_adjusted = max(10 ** (y_min_log - 1), 0.0001)  # Don't go below 0.0001
        
        # Set y_max to be greater than the actual maximum plotted value
        # Use the next power of 10 above the maximum value
        if y_max >= 1.0:
            y_max_adjusted = 10 ** (y_max_log + 1)
        else:
            # For values less than 1, go to the next higher value (like 1.0)
            y_max_adjusted = 1.0
        
        # Ensure y_max is at least 10,000 to accommodate very high relative errors
        # This ensures the plot can show all values clearly, including potential Naive results
        y_max_adjusted = max(y_max_adjusted, 10000.0)
        
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
        ax.set_yticklabels(labels, fontsize=26)
    
    # Legend
    ax.legend(fontsize=26, ncol=2, loc="upper left", columnspacing=0.5, frameon=False)
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"Plot saved as: {output_file}")
    
    plt.close()

def main():
    db_path = "algorithm_comparison_p3.db"
    
    if not os.path.exists(db_path):
        print(f"Database file {db_path} not found.")
        return
    
    # List of datasets to plot (modify this list as needed)
    datasets_to_plot = ["jester2-swapped"]
    
    print(f"Plotting datasets: {datasets_to_plot}")
    
    for dataset in datasets_to_plot:
        print(f"\nProcessing {dataset}...")
        
        # Get data from database
        algorithms, q_values, _ = get_data_from_db(db_path, dataset)
        
        if algorithms is None:
            print(f"No data found for {dataset}")
            continue
        
        print(f"Found results for {dataset} dataset, ε=1.0")
        print(f"Q values: {q_values}")
        print(f"Algorithms: {list(algorithms.keys())}")
        
        # Map algorithm names for display
        algorithm_mapping = {
            'OneR': 'OneR',
            'ADV': 'MRCN',
            'ADV+': 'MRCN+', 
            'ADV++': 'MRCN++'
        }
        
        # Create mapping from display names back to database names
        reverse_mapping = {
            'OneR': 'OneR',
            'MRCN': 'ADV',
            'MRCN+': 'ADV+', 
            'MRCN++': 'ADV++'
        }
        
        # Create ordered algorithm list: Naive, OneR, MRCN, MRCN+, MRCN++
        display_algorithms = ['Naive']
        if 'OneR' in algorithms:
            display_algorithms.append('OneR')
        for alg in ['ADV', 'ADV+', 'ADV++']:
            if alg in algorithms:
                display_algorithms.append(algorithm_mapping[alg])
        
        # Create new algorithms dict with display names
        display_algorithms_dict = {}
        for display_alg in display_algorithms:
            if display_alg == 'Naive':
                display_algorithms_dict[display_alg] = algorithms['Naive']
            elif display_alg == 'OneR':
                if 'OneR' in algorithms:
                    display_algorithms_dict[display_alg] = algorithms['OneR']
            else:
                db_alg = reverse_mapping[display_alg]
                if db_alg in algorithms:
                    display_algorithms_dict[display_alg] = algorithms[db_alg]
        
        # Print some sample data
        for q in q_values[:3]:  # Show first 3 Q values
            print(f"  Q={q}:")
            for alg in display_algorithms:
                if alg in display_algorithms_dict and q in display_algorithms_dict[alg]:
                    print(f"    {alg}: {display_algorithms_dict[alg][q]:.3f}")
        
        # Create the plot
        output_file = f"algorithm_comparison_{dataset}_eps1_FIXED.pdf"
        create_p2_plot(dataset, display_algorithms_dict, q_values, output_file)
        
        print(f"Completed plots for {dataset}")
        print("---")

if __name__ == "__main__":
    main()
