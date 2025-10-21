#!/usr/bin/env python3
"""
Create professional bar plot matching the style of the reference PDF.
Shows relative errors for all 4 algorithms as Q increases, with log scale.
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

def create_professional_barplot(df):
    """Create professional bar plot matching the reference style."""
    if df.empty:
        print("No data to visualize!")
        return
    
    # Set up the plot with professional styling
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define colors matching the reference style
    colors = {
        'Naive': '#E74C3C',      # Red
        'ADV': '#3498DB',        # Blue  
        'ADV+': '#2ECC71',       # Green
        'ADV++': '#F39C12'       # Orange
    }
    
    # Get unique Q values and algorithms
    q_values = sorted(df['q'].unique())
    algorithms = ['Naive', 'ADV', 'ADV+', 'ADV++']
    
    # Set up bar positions
    x = np.arange(len(q_values))
    width = 0.2
    
    # Create bars for each algorithm
    bars_dict = {}
    for i, algo in enumerate(algorithms):
        algo_data = df[df['algorithm'] == algo].sort_values('q')
        errors = [algo_data[algo_data['q'] == q]['relative_error'].iloc[0] if len(algo_data[algo_data['q'] == q]) > 0 else 0 for q in q_values]
        
        bars = ax.bar(x + i*width, errors, width, 
                     label=algo, color=colors[algo], alpha=0.8, 
                     edgecolor='black', linewidth=0.5)
        bars_dict[algo] = bars
    
    # Customize the plot with professional styling
    ax.set_xlabel('Q Value', fontsize=14, fontweight='bold')
    ax.set_ylabel('Relative Error (Log Scale)', fontsize=14, fontweight='bold')
    ax.set_title('Algorithm Performance Comparison on Jester2-Swapped Dataset\n(ε=1, 10 Rounds Each)', 
                fontsize=16, fontweight='bold', pad=20)
    
    # Set log scale
    ax.set_yscale('log')
    
    # Set axis limits for better visualization
    ax.set_ylim(1e-3, 1.0)
    
    # Customize legend
    ax.legend(loc='upper left', fontsize=12, frameon=True, 
             fancybox=True, shadow=True, framealpha=0.9)
    
    # Set x-axis
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels([f'Q={q}' for q in q_values], fontsize=12)
    ax.tick_params(axis='y', labelsize=12)
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='y', linestyle='--')
    
    # Add reference lines
    ax.axhline(y=0.01, color='gray', linestyle=':', alpha=0.7, linewidth=1)
    ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.7, linewidth=1)
    
    # Add value labels on bars (only for significant values)
    for algo, bars in bars_dict.items():
        for bar, q in zip(bars, q_values):
            height = bar.get_height()
            if height > 0.001:  # Only label significant values
                ax.text(bar.get_x() + bar.get_width()/2, height * 1.1, 
                       f'{height:.3f}', ha='center', va='bottom', 
                       fontsize=9, rotation=90)
    
    # Professional styling
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Set background
    ax.set_facecolor('#FAFAFA')
    fig.patch.set_facecolor('white')
    
    plt.tight_layout()
    plt.savefig('algorithm_comparison_jester2_professional.png', dpi=300, bbox_inches='tight')
    plt.savefig('algorithm_comparison_jester2_professional.pdf', bbox_inches='tight')
    plt.show()
    
    # Print detailed analysis
    print("\n" + "="*80)
    print("ALGORITHM PERFORMANCE ANALYSIS")
    print("="*80)
    
    # Performance summary
    print("\nPerformance Summary by Q Value:")
    print("-" * 50)
    for q in q_values:
        q_data = df[df['q'] == q]
        print(f"\nQ = {q}:")
        for algo in algorithms:
            algo_data = q_data[q_data['algorithm'] == algo]
            if len(algo_data) > 0:
                error = algo_data['relative_error'].iloc[0]
                print(f"  {algo:6s}: {error:.6f} ({error*100:.2f}%)")
    
    # Best algorithm for each Q
    print(f"\n" + "="*80)
    print("BEST ALGORITHM FOR EACH Q VALUE")
    print("="*80)
    for q in q_values:
        q_data = df[df['q'] == q]
        if len(q_data) > 0:
            best_algo = q_data.loc[q_data['relative_error'].idxmin(), 'algorithm']
            best_error = q_data['relative_error'].min()
            print(f"Q={q:2d}: {best_algo:6s} (Error: {best_error:.6f} = {best_error*100:.2f}%)")
    
    # Overall performance ranking
    print(f"\n" + "="*80)
    print("OVERALL PERFORMANCE RANKING")
    print("="*80)
    avg_errors = df.groupby('algorithm')['relative_error'].mean().sort_values()
    for i, (algo, error) in enumerate(avg_errors.items(), 1):
        print(f"{i}. {algo:6s}: {error:.6f} ({error*100:.2f}%)")
    
    return fig

def main():
    """Main function to create the professional bar plot."""
    print("Loading experimental data...")
    df = load_all_data()
    
    if df.empty:
        print("No data found!")
        return
    
    print(f"Loaded {len(df)} data points")
    print(f"Q values: {sorted(df['q'].unique())}")
    print(f"Algorithms: {df['algorithm'].unique()}")
    
    print("\nCreating professional bar plot...")
    create_professional_barplot(df)
    print("\nProfessional bar plot complete!")
    print("Files created: algorithm_comparison_jester2_professional.png/.pdf")

if __name__ == "__main__":
    main()

