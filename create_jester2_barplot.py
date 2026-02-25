#!/usr/bin/env python3
"""
Create professional bar plot for Jester2-Swapped dataset matching the exact style
from plot_p_3_from_db.py. Shows relative errors for all 4 algorithms as Q increases.
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from pathlib import Path

def plt_settings():
    """Set matplotlib settings to match Figure 12 (density plots)"""
    # Set publication-quality plot parameters
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

def extract_data_from_file(filepath):
    """Extract key metrics from a single result file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Extract Q value from filename
        q_match = re.search(r'q(\d+)_', filepath)
        q = int(q_match.group(1)) if q_match else None
        
        # Extract algorithm from filename
        if 'naive' in filepath:
            algorithm = 'Naive'
        elif 'advplusplus' in filepath:
            algorithm = 'ADV++'
        elif 'advplus' in filepath:
            algorithm = 'ADV+'
        elif 'adv' in filepath:
            algorithm = 'ADV'
        else:
            algorithm = 'Unknown'
        
        # Extract relative error
        rel_error_match = re.search(r'adv rel err = ([\d\.e\-\+]+)', content)
        rel_error = float(rel_error_match.group(1)) if rel_error_match else None
        
        return {
            'q': q,
            'algorithm': algorithm,
            'relative_error': rel_error
        }
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return None

def load_all_data(data_dir_str="jester-p-4-exp"):
    """Load data from all result files."""
    data_dir = Path(data_dir_str)
    all_data = []
    
    if not data_dir.exists():
        print(f"Directory {data_dir} not found!")
        return pd.DataFrame()
    
    for filepath in data_dir.glob("*.txt"):
        data = extract_data_from_file(str(filepath))
        if data:
            all_data.append(data)
    
    return pd.DataFrame(all_data)

def create_jester2_plot(results, algorithms):
    """Create the algorithm comparison plot with proper styling matching the reference"""
    
    plt_settings()
    
    # Prepare data
    q_values = sorted(results.keys())
    
    # Define colors and styles matching the reference plot
    method_styles = {
        'Naive': {'color': 'white', 'hatch': '', 'edgecolor': 'black'},
        'ADV': {'color': 'grey', 'hatch': '', 'edgecolor': 'black'},
        'ADV+': {'color': 'black', 'hatch': '////', 'edgecolor': 'black'},
        'ADV++': {'color': 'darkgrey', 'hatch': '//', 'edgecolor': 'black'}
    }
    
    # Create figure
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.15)  # Add more bottom margin
    
    # Plot bars for each algorithm
    x = np.arange(len(q_values))
    total_width, n = 0.8, len(algorithms)
    width = total_width / n
    
    for i, alg in enumerate(algorithms):
        values = []
        
        for j, q in enumerate(q_values):
            val = results[q].get(alg, float('inf'))
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
    """Main function to create the professional bar plot."""
    # Directories and corresponding epsilon values
    dir_epsilon_map = [
        ("jester-p-4-exp", 1.0),
        ("jester-p-4-exp-eps2", 2.0)
    ]
    
    for data_dir_str, epsilon in dir_epsilon_map:
        print(f"\n{'='*80}")
        print(f"Processing {data_dir_str}, ε={epsilon}")
        print(f"{'='*80}")
        
    print("Loading experimental data...")
        df = load_all_data(data_dir_str)
    
    if df.empty:
            print(f"No data found in {data_dir_str}!")
            continue
    
    print(f"Loaded {len(df)} data points")
    print(f"Q values: {sorted(df['q'].unique())}")
    print(f"Algorithms: {df['algorithm'].unique()}")
    
    # Organize data into results dictionary matching the reference format
    results = {}
    algorithms = ['Naive', 'ADV', 'ADV+', 'ADV++']
    
    for q in sorted(df['q'].unique()):
        results[q] = {}
        q_data = df[df['q'] == q]
        for algo in algorithms:
            algo_data = q_data[q_data['algorithm'] == algo]
            if len(algo_data) > 0:
                results[q][algo] = algo_data['relative_error'].iloc[0]
            else:
                results[q][algo] = float('inf')
    
    print("\nCreating professional bar plot...")
    fig = create_jester2_plot(results, algorithms)
    
    # Save the plot
        eps_str = str(int(epsilon)) if epsilon == int(epsilon) else str(epsilon).replace('.', '')
        output_filename = f'algorithm_comparison_jester2-swapped_eps{eps_str}_FIXED_p4.pdf'
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as: {output_filename}")
    plt.close()
    
    # Print summary
    print("\n" + "="*80)
        print(f"ALGORITHM PERFORMANCE SUMMARY (ε={epsilon})")
    print("="*80)
    
    for q in sorted(results.keys()):
        print(f"\nQ = {q}:")
        for alg in algorithms:
            val = results[q].get(alg, float('inf'))
            if val != float('inf'):
                print(f"  {alg:6s}: {val:.6f} ({val*100:.2f}%)")
            else:
                print(f"  {alg:6s}: No data")
    
    print("\nProfessional bar plot complete!")

if __name__ == "__main__":
    main()

