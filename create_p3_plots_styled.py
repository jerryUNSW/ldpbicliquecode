#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import re
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

def parse_p3_output(output_text, filename=""):
    """Parse the output from p3 new results format"""
    
    # Extract dataset name from filename
    if filename:
        # Extract dataset name from filename like "p3_lrcwiki_new.txt"
        match = re.search(r'p3_([^_]+(?:_[^_]+)*)_new\.txt', filename)
        if match:
            dataset = match.group(1)
        else:
            dataset = "unknown"
    else:
        dataset = "unknown"
    
    epsilon = 1.0  # Default epsilon for P=3 experiments
    
    # Parse individual algorithm runs
    # Algorithm mapping: 0=Naive, 2=MRCN, 3=MRCN+, 4=MRCN++
    alg_mapping = {0: 'Naive', 2: 'MRCN', 3: 'MRCN+', 4: 'MRCN++'}
    
    results = {}
    
    # Find all algorithm runs
    sections = re.split(r'DATASET: \w+\s+ALG=(\d+)\s+Q=(\d+)', output_text)
    
    for i in range(1, len(sections), 3):
        if i + 2 < len(sections):
            alg_num = int(sections[i])
            q_value = int(sections[i + 1])
            section_content = sections[i + 2]
            
            # Extract relative error from this section
            rel_err_match = re.search(r'adv rel err = ([0-9.e+-]+)', section_content)
            if rel_err_match:
                rel_error = float(rel_err_match.group(1))
                
                if q_value not in results:
                    results[q_value] = {}
                
                if alg_num in alg_mapping:
                    alg_name = alg_mapping[alg_num]
                    results[q_value][alg_name] = rel_error
    
    return dataset, epsilon, results

def create_p3_plot(dataset, epsilon, results):
    """Create the algorithm comparison plot with proper styling for p=3"""
    
    plt_settings()
    
    # Prepare data
    q_values = sorted(results.keys())
    algorithms = ['Naive', 'MRCN', 'MRCN+', 'MRCN++']
    
    # Define colors and styles matching the reference script
    method_styles = {
        'Naive': {'color': 'white', 'hatch': '', 'edgecolor': 'black'},
        'MRCN': {'color': 'grey', 'hatch': '', 'edgecolor': 'black'},
        'MRCN+': {'color': 'black', 'hatch': '////', 'edgecolor': 'black'},
        'MRCN++': {'color': 'darkgrey', 'hatch': '//', 'edgecolor': 'black'}
    }
    
    # Create figure with more bottom margin
    fig, ax = plt.subplots(figsize=(10, 7))
    plt.subplots_adjust(bottom=0.15)  # Add more bottom margin
    
    # Plot bars for each algorithm
    x = np.arange(len(q_values))
    total_width, n = 0.9, len(algorithms)
    width = total_width / (n + 1)
    
    # Calculate positions
    positions = [x + i * width - (n * width / 2) for i in range(n)]
    
    for i, alg in enumerate(algorithms):
        values = []
        valid_positions = []
        for j, q in enumerate(q_values):
            val = results[q].get(alg, float('inf'))
            # Handle infinite values by setting them to a large number for plotting
            val = val if val != float('inf') else 1e6
            values.append(val)
            valid_positions.append(positions[i][j])
        
        if values:  # Only plot if there are valid values
            style = method_styles[alg]
            ax.bar(valid_positions, values, width=width, 
                   label=alg, 
                   linewidth=0.5, 
                   edgecolor=style['edgecolor'],
                   fc=style['color'],
                   hatch=style['hatch'])
    
    # Customize plot
    ax.set_xlabel('q', fontsize=32)
    ax.set_ylabel('Mean Relative Error', fontsize=32)
    ax.set_title('(P = 3, ε = 1)', fontsize=30)
    ax.set_yscale('log')
    
    # Set x-axis ticks
    ax.set_xticks([x[i] - width for i in range(len(q_values))])
    ax.set_xticklabels([str(q) for q in q_values], fontsize=26)
    
    # Set y-axis limits and ticks - always use log scale
    all_values = []
    for q in q_values:
        for alg in algorithms:
            val = results[q].get(alg, float('inf'))
            if val != float('inf'):
                all_values.append(val)
    
    if all_values:
        y_min = np.min(all_values)
        y_max = np.max(all_values)
        
        # Print debug info
        print(f"Y-axis range: min={y_min:.3f}, max={y_max:.3f}")
        print(f"All values: {sorted(all_values)}")
        
        # Always use log scale with appropriate range
        y_min_adjusted = y_min / 10 if y_min > 0 else y_min * 10
        y_max_adjusted = y_max * 10
        
        ax.set_ylim(y_min_adjusted, y_max_adjusted)
        
        # Set y-axis ticks - always use log scale
        indices = [i for i in range(math.floor(math.log10(y_min_adjusted)), 
                                   math.ceil(math.log10(y_max_adjusted)) + 1, 2)]
        pos, labels = get_pos_and_labels(indices)
        ax.set_yticks(pos)
        ax.set_yticklabels(labels, fontsize=26)
    
    # Legend
    ax.legend(fontsize=26, ncol=2, loc="upper center", columnspacing=0.5, frameon=False)
    
    plt.tight_layout()
    return fig

def main():
    # Results directory
    results_dir = "larger_datasets_batch_results_FIXED"
    
    # Datasets to process (p=3 results) - only include datasets that have _new.txt files
    datasets = ["csbwiki", "lrcwiki", "nips", "rmwiki"]
    
    for dataset in datasets:
        # Read the output from the p3 new results
        output_file = f"{results_dir}/p3_{dataset}_new.txt"
        
        if not os.path.exists(output_file):
            print(f"Output file {output_file} not found. Skipping {dataset}.")
            continue
            
        try:
            with open(output_file, 'r') as f:
                output_text = f.read()
        except Exception as e:
            print(f"Error reading {output_file}: {e}")
            continue
        
        # Parse the output
        dataset_name, epsilon, results = parse_p3_output(output_text, output_file)
        
        if results is None:
            print(f"Failed to parse results for {dataset}")
            continue
        
        print(f"Parsed results for {dataset_name} dataset, ε={epsilon}")
        print(f"Q values: {sorted(results.keys())}")
        print(f"Algorithms: {list(results[list(results.keys())[0]].keys())}")
        
        # Create and save main plot
        fig = create_p3_plot(dataset_name, epsilon, results)
        output_filename = f'algorithm_comparison_{dataset_name}_eps{int(epsilon)}_FIXED_p3.pdf'
        fig.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as: {output_filename}")
        plt.close()
        
        print(f"Completed plots for {dataset_name}")
        print("---")

if __name__ == "__main__":
    main()
