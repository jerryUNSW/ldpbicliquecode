#!/usr/bin/env python3
"""
Create bar plot showing relative errors for all 4 algorithms as Q increases.
Uses log scale for better visualization of the wide range of errors.
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

def create_barplot(df):
    """Create bar plot showing relative errors vs Q for all algorithms."""
    if df.empty:
        print("No data to visualize!")
        return
    
    # Create figure
    plt.figure(figsize=(14, 10))
    
    # Define colors for each algorithm
    colors = {'Naive': '#FF6B6B', 'ADV': '#4ECDC4', 'ADV+': '#45B7D1', 'ADV++': '#96CEB4'}
    
    # Get unique Q values and algorithms
    q_values = sorted(df['q'].unique())
    algorithms = ['Naive', 'ADV', 'ADV+', 'ADV++']
    
    # Set up bar positions
    x = np.arange(len(q_values))
    width = 0.2
    
    # Create bars for each algorithm
    for i, algo in enumerate(algorithms):
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        errors = [algo_data[algo_data['q'] == q]['relative_error'].iloc[0] if len(algo_data[algo_data['q'] == q]) > 0 else 0 for q in q_values]
        
        bars = plt.bar(x + i*width, errors, width, 
                      label=algo, color=colors[algo], alpha=0.8)
        
        # Add value labels on bars
        for bar, error in zip(bars, errors):
            if error > 0:
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                        f'{error:.3f}', ha='center', va='bottom', fontsize=8, rotation=90)
    
    # Customize the plot
    plt.xlabel('Q Value', fontsize=14, fontweight='bold')
    plt.ylabel('Relative Error (Log Scale)', fontsize=14, fontweight='bold')
    plt.title('Relative Error vs Q Value by Algorithm\n(10 Rounds Each)', fontsize=16, fontweight='bold', pad=20)
    plt.yscale('log')
    plt.legend(fontsize=12, loc='upper left')
    plt.grid(True, alpha=0.3, axis='y')
    
    # Set x-axis
    plt.xticks(x + width * 1.5, [f'Q={q}' for q in q_values], fontsize=12)
    plt.yticks(fontsize=12)
    
    # Set y-axis limits for better visualization
    plt.ylim(1e-3, 1.0)
    
    # Add horizontal lines for reference
    plt.axhline(y=0.01, color='gray', linestyle='--', alpha=0.5, label='1% error')
    plt.axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, label='10% error')
    
    plt.tight_layout()
    plt.savefig('relative_error_barplot.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summary statistics
    print("\n" + "="*80)
    print("RELATIVE ERROR SUMMARY BY Q VALUE")
    print("="*80)
    
    for q in q_values:
        print(f"\nQ = {q}:")
        q_data = df[df['q'] == q]
        for algo in algorithms:
            algo_data = q_data[q_data['algorithm'] == algo]
            if len(algo_data) > 0:
                error = algo_data['relative_error'].iloc[0]
                print(f"  {algo:6s}: {error:.6f} ({error*100:.2f}%)")
            else:
                print(f"  {algo:6s}: No data")
    
    # Print best algorithm for each Q
    print(f"\n" + "="*80)
    print("BEST ALGORITHM FOR EACH Q VALUE")
    print("="*80)
    
    for q in q_values:
        q_data = df[df['q'] == q]
        if len(q_data) > 0:
            best_algo = q_data.loc[q_data['relative_error'].idxmin(), 'algorithm']
            best_error = q_data['relative_error'].min()
            print(f"Q={q:2d}: {best_algo:6s} (Error: {best_error:.6f} = {best_error*100:.2f}%)")
    
    return plt.gcf()

def main():
    """Main function to create the bar plot."""
    print("Loading experimental data...")
    df = load_all_data()
    
    if df.empty:
        print("No data found!")
        return
    
    print(f"Loaded {len(df)} data points")
    print(f"Q values: {sorted(df['q'].unique())}")
    print(f"Algorithms: {df['algorithm'].unique()}")
    
    print("\nCreating bar plot...")
    create_barplot(df)
    print("\nBar plot complete! Check 'relative_error_barplot.png'")

if __name__ == "__main__":
    main()

