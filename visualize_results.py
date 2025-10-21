#!/usr/bin/env python3
"""
Data visualization script for biclique counting experiment results.
Extracts data from jester-p-4-exp/ directory and creates comprehensive visualizations.
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
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
        
        # Extract time
        time_match = re.search(r'time:([\d\.]+)', content)
        time = float(time_match.group(1)) if time_match else None
        
        # Extract real count
        real_count_match = re.search(r'real count = ([\d\.e\-\+]+)', content)
        real_count = float(real_count_match.group(1)) if real_count_match else None
        
        # Extract estimate
        estimate_match = re.search(r'estimate = ([\d\.e\-\+]+)', content)
        estimate = float(estimate_match.group(1)) if estimate_match else None
        
        return {
            'q': q,
            'algorithm': algorithm,
            'relative_error': rel_error,
            'time': time,
            'real_count': real_count,
            'estimate': estimate
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

def create_visualizations(df):
    """Create comprehensive visualizations."""
    if df.empty:
        print("No data to visualize!")
        return
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 15))
    
    # 1. Relative Error vs Q for different algorithms
    plt.subplot(2, 3, 1)
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        plt.plot(algo_data['q'], algo_data['relative_error'], 
                marker='o', linewidth=2, label=algo, markersize=8)
    
    plt.xlabel('Q Value', fontsize=12)
    plt.ylabel('Relative Error', fontsize=12)
    plt.title('Relative Error vs Q Value by Algorithm', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    # 2. Execution Time vs Q for different algorithms
    plt.subplot(2, 3, 2)
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        plt.plot(algo_data['q'], algo_data['time'], 
                marker='s', linewidth=2, label=algo, markersize=8)
    
    plt.xlabel('Q Value', fontsize=12)
    plt.ylabel('Execution Time (seconds)', fontsize=12)
    plt.title('Execution Time vs Q Value by Algorithm', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 3. Relative Error comparison (bar chart)
    plt.subplot(2, 3, 3)
    error_by_algo = df.groupby('algorithm')['relative_error'].mean()
    bars = plt.bar(error_by_algo.index, error_by_algo.values, 
                   color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4'])
    plt.xlabel('Algorithm', fontsize=12)
    plt.ylabel('Average Relative Error', fontsize=12)
    plt.title('Average Relative Error by Algorithm', fontsize=14, fontweight='bold')
    plt.yscale('log')
    plt.xticks(rotation=45)
    
    # Add value labels on bars
    for bar, value in zip(bars, error_by_algo.values):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                f'{value:.3f}', ha='center', va='bottom', fontsize=10)
    
    # 4. Execution Time comparison (bar chart)
    plt.subplot(2, 3, 4)
    time_by_algo = df.groupby('algorithm')['time'].mean()
    bars = plt.bar(time_by_algo.index, time_by_algo.values,
                   color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4'])
    plt.xlabel('Algorithm', fontsize=12)
    plt.ylabel('Average Execution Time (seconds)', fontsize=12)
    plt.title('Average Execution Time by Algorithm', fontsize=14, fontweight='bold')
    plt.xticks(rotation=45)
    
    # Add value labels on bars
    for bar, value in zip(bars, time_by_algo.values):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                f'{value:.1f}s', ha='center', va='bottom', fontsize=10)
    
    # 5. Heatmap of relative errors
    plt.subplot(2, 3, 5)
    pivot_data = df.pivot(index='algorithm', columns='q', values='relative_error')
    sns.heatmap(pivot_data, annot=True, fmt='.3f', cmap='RdYlBu_r', 
                cbar_kws={'label': 'Relative Error'})
    plt.title('Relative Error Heatmap', fontsize=14, fontweight='bold')
    plt.xlabel('Q Value', fontsize=12)
    plt.ylabel('Algorithm', fontsize=12)
    
    # 6. Scatter plot: Time vs Relative Error
    plt.subplot(2, 3, 6)
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo]
        plt.scatter(algo_data['time'], algo_data['relative_error'], 
                   label=algo, s=100, alpha=0.7)
    
    plt.xlabel('Execution Time (seconds)', fontsize=12)
    plt.ylabel('Relative Error', fontsize=12)
    plt.title('Time vs Accuracy Trade-off', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig('biclique_experiment_results.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create a summary statistics table
    print("\n" + "="*80)
    print("EXPERIMENT SUMMARY STATISTICS")
    print("="*80)
    
    summary_stats = df.groupby('algorithm').agg({
        'relative_error': ['mean', 'std', 'min', 'max'],
        'time': ['mean', 'std', 'min', 'max']
    }).round(6)
    
    print(summary_stats)
    
    # Create detailed comparison table
    print("\n" + "="*80)
    print("DETAILED RESULTS BY Q VALUE")
    print("="*80)
    
    detailed_results = df.pivot_table(
        index='q', 
        columns='algorithm', 
        values=['relative_error', 'time'], 
        aggfunc='mean'
    ).round(6)
    
    print(detailed_results)
    
    return fig

def main():
    """Main function to run the visualization."""
    print("Loading experimental data...")
    df = load_all_data()
    
    if df.empty:
        print("No data found! Make sure the jester-p-4-exp/ directory exists and contains result files.")
        return
    
    print(f"Loaded {len(df)} data points")
    print(f"Algorithms: {df['algorithm'].unique()}")
    print(f"Q values: {sorted(df['q'].unique())}")
    
    print("\nCreating visualizations...")
    create_visualizations(df)
    
    print("\nVisualization complete! Check 'biclique_experiment_results.png' for the plots.")

if __name__ == "__main__":
    main()
