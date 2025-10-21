#!/usr/bin/env python3
"""
Scalability analysis for biclique counting algorithms.
Analyzes how algorithms perform as problem size (Q) increases.
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
        
        return {
            'q': q,
            'algorithm': algorithm,
            'relative_error': rel_error,
            'time': time,
            'real_count': real_count
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

def create_scalability_analysis(df):
    """Create scalability analysis visualizations."""
    if df.empty:
        print("No data to analyze!")
        return
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Error Growth vs Q
    ax1 = axes[0, 0]
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        ax1.plot(algo_data['q'], algo_data['relative_error'], 
                marker='o', linewidth=2, label=algo, markersize=8)
    
    ax1.set_xlabel('Q Value (Problem Size)')
    ax1.set_ylabel('Relative Error')
    ax1.set_title('Error Growth vs Problem Size', fontweight='bold', fontsize=14)
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Time Complexity Analysis
    ax2 = axes[0, 1]
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        ax2.plot(algo_data['q'], algo_data['time'], 
                marker='s', linewidth=2, label=algo, markersize=8)
    
    ax2.set_xlabel('Q Value (Problem Size)')
    ax2.set_ylabel('Execution Time (seconds)')
    ax2.set_title('Time Complexity vs Problem Size', fontweight='bold', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Efficiency Analysis (Error/Time ratio)
    ax3 = axes[1, 0]
    df['efficiency'] = df['relative_error'] / df['time']
    
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        ax3.plot(algo_data['q'], algo_data['efficiency'], 
                marker='^', linewidth=2, label=algo, markersize=8)
    
    ax3.set_xlabel('Q Value (Problem Size)')
    ax3.set_ylabel('Efficiency (Error/Time)')
    ax3.set_title('Algorithm Efficiency vs Problem Size', fontweight='bold', fontsize=14)
    ax3.set_yscale('log')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Scalability Ranking
    ax4 = axes[1, 1]
    
    # Calculate scalability score (how well algorithm maintains performance as Q increases)
    scalability_scores = {}
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        if len(algo_data) > 1:
            # Calculate slope of error vs Q (lower slope = better scalability)
            q_values = algo_data['q'].values
            error_values = np.log10(algo_data['relative_error'].values)
            slope = np.polyfit(q_values, error_values, 1)[0]
            scalability_scores[algo] = -slope  # Negative slope is better
    
    algorithms = list(scalability_scores.keys())
    scores = list(scalability_scores.values())
    
    bars = ax4.bar(algorithms, scores, color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4'])
    ax4.set_xlabel('Algorithm')
    ax4.set_ylabel('Scalability Score (Higher = Better)')
    ax4.set_title('Algorithm Scalability Ranking', fontweight='bold', fontsize=14)
    ax4.grid(True, alpha=0.3)
    
    # Add value labels
    for bar, score in zip(bars, scores):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1, 
                f'{score:.3f}', ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('scalability_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print scalability analysis
    print("\n" + "="*80)
    print("SCALABILITY ANALYSIS")
    print("="*80)
    
    print("\nAlgorithm Performance by Q Value:")
    print("-" * 50)
    for q in sorted(df['q'].unique()):
        q_data = df[df['q'] == q]
        print(f"\nQ = {q}:")
        for _, row in q_data.iterrows():
            print(f"  {row['algorithm']:6s}: Error={row['relative_error']:.6f}, Time={row['time']:.2f}s")
    
    print(f"\nScalability Scores (Higher = Better):")
    print("-" * 40)
    for algo, score in sorted(scalability_scores.items(), key=lambda x: x[1], reverse=True):
        print(f"{algo:6s}: {score:.4f}")
    
    # Performance degradation analysis
    print(f"\nPerformance Degradation Analysis:")
    print("-" * 40)
    for algo in df['algorithm'].unique():
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        if len(algo_data) > 1:
            min_error = algo_data['relative_error'].min()
            max_error = algo_data['relative_error'].max()
            degradation = (max_error - min_error) / min_error * 100
            print(f"{algo:6s}: {degradation:.1f}% error increase from Q{algo_data['q'].min()} to Q{algo_data['q'].max()}")
    
    return fig

def main():
    """Main function to run the scalability analysis."""
    print("Loading experimental data for scalability analysis...")
    df = load_all_data()
    
    if df.empty:
        print("No data found!")
        return
    
    print(f"Loaded {len(df)} data points")
    print("Creating scalability analysis...")
    create_scalability_analysis(df)
    print("\nScalability analysis complete! Check 'scalability_analysis.png' for the plots.")

if __name__ == "__main__":
    main()


