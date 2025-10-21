#!/usr/bin/env python3
"""
Focused visualization: Relative Error vs Q Value by Algorithm
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

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
        rel_error_match = re.search(r'relative error = ([\d\.e\-\+]+)', content)
        rel_error = float(rel_error_match.group(1)) if rel_error_match else None
        
        return {
            'q': q,
            'algorithm': algorithm,
            'relative_error': rel_error
        }
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return None

def load_all_data():
    """Load data from all result files."""
    data_dir = Path("jester-p-4-exp")
    all_data = []
    
    if not data_dir.exists():
        print(f"Directory {data_dir} not found!")
        return pd.DataFrame()
    
    for filepath in data_dir.glob("*.txt"):
        data = extract_data_from_file(str(filepath))
        if data:
            all_data.append(data)
    
    return pd.DataFrame(all_data)

def create_error_vs_q_plot(df):
    """Create the focused Relative Error vs Q Value plot."""
    if df.empty:
        print("No data to visualize!")
        return
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Define colors and markers for each algorithm
    colors = {'Naive': '#FF6B6B', 'ADV': '#4ECDC4', 'ADV+': '#45B7D1', 'ADV++': '#96CEB4'}
    markers = {'Naive': 'o', 'ADV': 's', 'ADV+': '^', 'ADV++': 'D'}
    
    # Plot each algorithm
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        plt.plot(algo_data['q'], algo_data['relative_error'], 
                marker=markers[algo], linewidth=3, label=algo, 
                markersize=10, color=colors[algo], alpha=0.8)
    
    # Customize the plot
    plt.xlabel('Q Value', fontsize=14, fontweight='bold')
    plt.ylabel('Relative Error', fontsize=14, fontweight='bold')
    plt.title('Relative Error vs Q Value by Algorithm', fontsize=16, fontweight='bold', pad=20)
    plt.yscale('log')
    plt.legend(fontsize=12, loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Set axis limits and ticks
    plt.xlim(1.5, 10.5)
    plt.xticks(range(2, 11), fontsize=12)
    plt.yticks(fontsize=12)
    
    # Add minor grid
    plt.grid(True, which='minor', alpha=0.2)
    
    # Tight layout and save
    plt.tight_layout()
    plt.savefig('error_vs_q_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summary statistics
    print("\n" + "="*60)
    print("RELATIVE ERROR vs Q VALUE SUMMARY")
    print("="*60)
    
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        print(f"\n{algo}:")
        for _, row in algo_data.iterrows():
            print(f"  Q={row['q']:2d}: {row['relative_error']:.6f}")

def main():
    """Main function to create the focused plot."""
    print("Loading experimental data...")
    df = load_all_data()
    
    if df.empty:
        print("No data found!")
        return
    
    print(f"Loaded {len(df)} data points")
    print("Creating Relative Error vs Q Value plot...")
    create_error_vs_q_plot(df)
    print("\nPlot complete! Check 'error_vs_q_plot.png'")

if __name__ == "__main__":
    main()


