#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import re
import os
import math

def plt_settings():
    """Set matplotlib settings to match old plotting style"""
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2
    plt.rcParams['xtick.minor.width'] = 0.8
    plt.rcParams['ytick.minor.width'] = 0.8

def get_pos_and_labels(indices):
    """Get positions and labels for y-axis ticks"""
    pos = [10**i for i in indices]
    labels = [f'$10^{{{i}}}$' for i in indices]
    return pos, labels

def parse_output(output_text):
    """Parse the output from test_batch_with_ground_truth.cpp"""
    
    # Extract configuration
    config_match = re.search(r'Dataset: (.+)\nEpsilon: (.+)\nRounds: (.+)', output_text)
    if config_match:
        dataset = config_match.group(1).split('/')[-1]  # Extract just the dataset name
        epsilon = float(config_match.group(2))
        rounds = int(config_match.group(3))
    else:
        dataset = "unknown"
        epsilon = 1.0
        rounds = 3
    
    # Extract ground truth
    ground_truth = {}
    gt_pattern = r'\(2,(\d+)\)-bicliques: (\d+)'
    for match in re.finditer(gt_pattern, output_text):
        q = int(match.group(1))
        count = int(match.group(2))
        ground_truth[q] = count
    
    # Extract algorithm comparison summary
    summary_section = re.search(r'=== Algorithm Comparison Summary ===(.*?)(?=\n\n|\Z)', output_text, re.DOTALL)
    if not summary_section:
        print("Could not find algorithm comparison summary")
        return None, None, None, None
    
    summary_text = summary_section.group(1)
    
    # Parse the summary table
    lines = summary_text.strip().split('\n')
    header_line = lines[1]  # "Q\tNaive\t\toneR\t\tADV\t\tADV+\t\tADV++"
    data_lines = lines[2:]  # Skip header and separator line
    
    # Extract algorithm names from header - use proper names
    algorithms = ['Naive', 'OneR', 'MRCN', 'MRCN+', 'MRCN++']
    
    # Parse data
    results = {}
    for line in data_lines:
        if line.strip():
            parts = line.split('\t')
            if len(parts) >= 6:
                q = int(parts[0])
                results[q] = {}
                for i, alg in enumerate(algorithms):
                    if i + 1 < len(parts):
                        # Parse scientific notation
                        value_str = parts[i + 1].strip()
                        try:
                            results[q][alg] = float(value_str)
                        except ValueError:
                            results[q][alg] = float('inf')
    
    return dataset, epsilon, results, ground_truth

def create_plot(dataset, epsilon, results, ground_truth):
    """Create the algorithm comparison plot with old script styling"""
    
    plt_settings()
    
    # Prepare data
    q_values = sorted(results.keys())
    algorithms = ['Naive', 'OneR', 'MRCN', 'MRCN+', 'MRCN++']
    
    # Define colors and styles matching old scripts
    method_styles = {
        'Naive': {'color': 'white', 'hatch': '', 'edgecolor': 'black'},
        'OneR': {'color': 'lightgray', 'hatch': '', 'edgecolor': 'black'},
        'MRCN': {'color': 'grey', 'hatch': '', 'edgecolor': 'black'},
        'MRCN+': {'color': 'black', 'hatch': '////', 'edgecolor': 'black'},
        'MRCN++': {'color': 'darkgrey', 'hatch': '//', 'edgecolor': 'black'}
    }
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot bars for each algorithm
    x = np.arange(len(q_values))
    total_width, n = 0.9, len(algorithms)
    width = total_width / (n + 1)
    
    # Calculate positions
    positions = [x + i * width - (n * width / 2) for i in range(n)]
    
    for i, alg in enumerate(algorithms):
        values = [results[q].get(alg, float('inf')) for q in q_values]
        # Handle infinite values by setting them to a large number for plotting
        values = [v if v != float('inf') else 1e6 for v in values]
        
        style = method_styles[alg]
        ax.bar(positions[i], values, width=width, 
               label=alg, 
               linewidth=0.5, 
               edgecolor=style['edgecolor'],
               fc=style['color'],
               hatch=style['hatch'])
    
    # Customize plot
    ax.set_xlabel('q', fontsize=20)
    ax.set_ylabel('Mean Relative Error', fontsize=20)
    ax.set_title(f'{dataset.upper()} Dataset (P=2, ε={epsilon}) - FIXED', fontsize=20)
    ax.set_yscale('log')
    
    # Set x-axis ticks
    ax.set_xticks(x - width)
    ax.set_xticklabels([str(q) for q in q_values], fontsize=18)
    
    # Set y-axis limits and ticks
    all_values = []
    for q in q_values:
        for alg in algorithms:
            val = results[q].get(alg, float('inf'))
            if val != float('inf'):
                all_values.append(val)
    
    if all_values:
        y_min = np.min(all_values)
        y_max = np.max(all_values)
        y_min_adjusted = y_min / 10 if y_min > 0 else y_min * 10
        y_max_adjusted = y_max * 10
        ax.set_ylim(y_min_adjusted, y_max_adjusted)
        
        # Set y-axis ticks
        indices = [i for i in range(math.floor(math.log10(y_min_adjusted)), 
                                   math.ceil(math.log10(y_max_adjusted)) + 1, 2)]
        pos, labels = get_pos_and_labels(indices)
        ax.set_yticks(pos)
        ax.set_yticklabels(labels, fontsize=20)
    
    # Legend
    ax.legend(fontsize=20, ncol=5, loc="upper center", columnspacing=0.5, frameon=False)
    
    plt.tight_layout()
    return fig


def main():
    # Results directory
    results_dir = "larger_datasets_batch_results_FIXED"
    
    # Datasets to process
    datasets = ["co", "unicode", "librec-filmtrust-ratings", "bpywiki", "csbwiki", "rmwiki", "lrcwiki", "nips"]
    
    for dataset in datasets:
        # Read the output from the test
        output_file = f"{results_dir}/{dataset}_eps1_results.txt"
        
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
        dataset_name, epsilon, results, ground_truth = parse_output(output_text)
        
        if results is None:
            print(f"Failed to parse results for {dataset}")
            continue
        
        print(f"Parsed results for {dataset_name} dataset, ε={epsilon}")
        print(f"Q values: {sorted(results.keys())}")
        print(f"Algorithms: {list(results[list(results.keys())[0]].keys())}")
        
        # Create and save main plot
        fig = create_plot(dataset_name, epsilon, results, ground_truth)
        output_filename = f'algorithm_comparison_{dataset_name}_eps{int(epsilon)}_FIXED.pdf'
        fig.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as: {output_filename}")
        plt.close()
        
        print(f"Completed plots for {dataset_name}")
        print("---")

if __name__ == "__main__":
    main()
