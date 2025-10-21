#!/usr/bin/env python3
"""
Detailed analysis script for biclique counting experiment results.
Creates focused visualizations and statistical analysis.
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

def create_focused_analysis(df):
    """Create focused analysis visualizations."""
    if df.empty:
        print("No data to analyze!")
        return
    
    # Create figure with focused subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Algorithm Performance Comparison (Bar Chart)
    ax1 = axes[0, 0]
    algo_performance = df.groupby('algorithm').agg({
        'relative_error': 'mean',
        'time': 'mean'
    }).reset_index()
    
    x = np.arange(len(algo_performance))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, algo_performance['relative_error'], width, 
                    label='Avg Relative Error', alpha=0.8, color='#FF6B6B')
    ax1.set_xlabel('Algorithm')
    ax1.set_ylabel('Average Relative Error')
    ax1.set_title('Algorithm Accuracy Comparison', fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(algo_performance['algorithm'])
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add value labels
    for bar, value in zip(bars1, algo_performance['relative_error']):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                f'{value:.3f}', ha='center', va='bottom', fontsize=9)
    
    # 2. Time vs Accuracy Scatter
    ax2 = axes[0, 1]
    colors = {'Naive': '#FF6B6B', 'ADV': '#4ECDC4', 'ADV+': '#45B7D1', 'ADV++': '#96CEB4'}
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo]
        ax2.scatter(algo_data['time'], algo_data['relative_error'], 
                   label=algo, s=100, color=colors[algo], alpha=0.8)
    
    ax2.set_xlabel('Execution Time (seconds)')
    ax2.set_ylabel('Relative Error')
    ax2.set_title('Time vs Accuracy Trade-off', fontweight='bold')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Q Value Impact Analysis
    ax3 = axes[1, 0]
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        ax3.plot(algo_data['q'], algo_data['relative_error'], 
                marker='o', linewidth=2, label=algo, markersize=6)
    
    ax3.set_xlabel('Q Value')
    ax3.set_ylabel('Relative Error')
    ax3.set_title('Error vs Q Value (Problem Size)', fontweight='bold')
    ax3.set_yscale('log')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Performance Ranking
    ax4 = axes[1, 1]
    
    # Calculate performance score (lower is better): error * time
    df['performance_score'] = df['relative_error'] * df['time']
    ranking = df.groupby('algorithm')['performance_score'].mean().sort_values()
    
    bars = ax4.bar(range(len(ranking)), ranking.values, 
                   color=['#96CEB4', '#45B7D1', '#4ECDC4', '#FF6B6B'])
    ax4.set_xlabel('Algorithm (Ranked by Performance)')
    ax4.set_ylabel('Performance Score (Error × Time)')
    ax4.set_title('Algorithm Performance Ranking', fontweight='bold')
    ax4.set_xticks(range(len(ranking)))
    ax4.set_xticklabels(ranking.index, rotation=45)
    ax4.set_yscale('log')
    ax4.grid(True, alpha=0.3)
    
    # Add value labels
    for i, (bar, value) in enumerate(zip(bars, ranking.values)):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                f'{value:.2f}', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('focused_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print detailed analysis
    print("\n" + "="*80)
    print("DETAILED PERFORMANCE ANALYSIS")
    print("="*80)
    
    # Best performing algorithm for each Q
    print("\nBest Algorithm for Each Q Value:")
    print("-" * 40)
    for q in sorted(df['q'].unique()):
        q_data = df[df['q'] == q]
        best_algo = q_data.loc[q_data['relative_error'].idxmin(), 'algorithm']
        best_error = q_data['relative_error'].min()
        print(f"Q={q:2d}: {best_algo:6s} (Error: {best_error:.6f})")
    
    # Statistical significance analysis
    print("\nStatistical Analysis:")
    print("-" * 40)
    naive_errors = df[df['algorithm'] == 'Naive']['relative_error']
    adv_errors = df[df['algorithm'] == 'ADV']['relative_error']
    advplus_errors = df[df['algorithm'] == 'ADV+']['relative_error']
    advplusplus_errors = df[df['algorithm'] == 'ADV++']['relative_error']
    
    print(f"Naive avg error:     {naive_errors.mean():.6f} ± {naive_errors.std():.6f}")
    print(f"ADV avg error:       {adv_errors.mean():.6f} ± {adv_errors.std():.6f}")
    print(f"ADV+ avg error:      {advplus_errors.mean():.6f} ± {advplus_errors.std():.6f}")
    print(f"ADV++ avg error:     {advplusplus_errors.mean():.6f} ± {advplusplus_errors.std():.6f}")
    
    # Improvement ratios
    print(f"\nImprovement over Naive:")
    print(f"ADV:      {naive_errors.mean() / adv_errors.mean():.2f}x better")
    print(f"ADV+:     {naive_errors.mean() / advplus_errors.mean():.2f}x better")
    print(f"ADV++:    {naive_errors.mean() / advplusplus_errors.mean():.2f}x better")
    
    return fig

def main():
    """Main function to run the detailed analysis."""
    print("Loading experimental data for detailed analysis...")
    df = load_all_data()
    
    if df.empty:
        print("No data found!")
        return
    
    print(f"Loaded {len(df)} data points")
    print("Creating focused analysis...")
    create_focused_analysis(df)
    print("\nFocused analysis complete! Check 'focused_analysis.png' for the plots.")

if __name__ == "__main__":
    main()
