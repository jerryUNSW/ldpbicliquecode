#!/usr/bin/env python3
"""
Create bar plot comparing ADV and Naive algorithm relative errors for NIPS dataset.
"""

import matplotlib.pyplot as plt
import numpy as np

# Apply professional styling
plt.rcParams['savefig.dpi'] = 300 
plt.rcParams['figure.dpi'] = 300 
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams["legend.framealpha"] = 0
plt.rcParams["legend.handletextpad"] = 0.1
plt.rcParams["legend.columnspacing"] = 0.2
plt.rcParams['pdf.fonttype'] = 42

def plt_settings():
    """Apply styling from biclique plots."""
    plt.style.use('default')
    plt.rcParams['savefig.dpi'] = 300 
    plt.rcParams['figure.dpi'] = 300 
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams["legend.framealpha"] = 0
    plt.rcParams["legend.handletextpad"] = 0.1
    plt.rcParams["legend.columnspacing"] = 0.2
    plt.rcParams['pdf.fonttype'] = 42

# Data from the user
q_values = [4, 5, 6, 7, 8, 9, 10]
adv_errors = [37.43, 489.9, 7181.0, 138800.0, 3338000.0, 88970000.0, 2435000000.0]
naive_errors = [1557.307, 7404.205, 33927.06, 198384.0, 680666.3, 4067774.0, 16274870.0]

# Convert to percentage for better readability
adv_errors_pct = [err * 100 for err in adv_errors]
naive_errors_pct = [err * 100 for err in naive_errors]

def create_comparison_plot():
    """Create bar plot comparing ADV and Naive errors."""
    plt_settings()
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    x = np.arange(len(q_values))
    width = 0.35
    
    # Create bars
    bars1 = ax.bar(x - width/2, adv_errors_pct, width, 
                   label='ADV', color='lightblue', edgecolor='black', linewidth=0.5)
    bars2 = ax.bar(x + width/2, naive_errors_pct, width, 
                   label='Naive', color='lightcoral', edgecolor='black', linewidth=0.5)
    
    # Customize the plot
    ax.set_xlabel('Q Value', fontsize=16)
    ax.set_ylabel('Relative Error (%)', fontsize=16)
    ax.set_title('Algorithm Comparison: ADV vs Naive\nNIPS Dataset (3,Q)-Biclique Counting, ε=1.0', fontsize=18, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([f'Q={q}' for q in q_values], fontsize=14)
    ax.set_yscale('log')
    
    # Add value labels on bars
    def add_value_labels(bars, values):
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{value:.1e}%', ha='center', va='bottom', fontsize=10, rotation=90)
    
    # Add value labels (only for smaller values to avoid clutter)
    for i, (bar, value) in enumerate(zip(bars1, adv_errors_pct)):
        if value < 1e6:  # Only label smaller values
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                   f'{value:.1e}%', ha='center', va='bottom', fontsize=9, rotation=90)
    
    for i, (bar, value) in enumerate(zip(bars2, naive_errors_pct)):
        if value < 1e6:  # Only label smaller values
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                   f'{value:.1e}%', ha='center', va='bottom', fontsize=9, rotation=90)
    
    # Add legend
    ax.legend(fontsize=14, loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Set y-axis limits and ticks
    all_errors = adv_errors_pct + naive_errors_pct
    y_min = min(all_errors) / 10
    y_max = max(all_errors) * 10
    ax.set_ylim(y_min, y_max)
    
    # Create log scale ticks
    import math
    min_exp = math.floor(math.log10(y_min))
    max_exp = math.ceil(math.log10(y_max))
    tick_values = [10**i for i in range(min_exp, max_exp + 1, 2)]
    tick_labels = [f'$10^{{{i}}}$' for i in range(min_exp, max_exp + 1, 2)]
    ax.set_yticks(tick_values)
    ax.set_yticklabels(tick_labels, fontsize=12)
    
    plt.tight_layout()
    
    # Save the plot
    output_file = '/data/yizhangh/ldp-pq/algorithm_error_comparison.pdf'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Bar plot saved to: {output_file}")
    
    # Print summary statistics
    print(f"\n=== Summary Statistics ===")
    print(f"ADV Algorithm:")
    print(f"  Mean Error: {np.mean(adv_errors_pct):.2e}%")
    print(f"  Median Error: {np.median(adv_errors_pct):.2e}%")
    print(f"  Min Error: {np.min(adv_errors_pct):.2e}%")
    print(f"  Max Error: {np.max(adv_errors_pct):.2e}%")
    
    print(f"\nNaive Algorithm:")
    print(f"  Mean Error: {np.mean(naive_errors_pct):.2e}%")
    print(f"  Median Error: {np.median(naive_errors_pct):.2e}%")
    print(f"  Min Error: {np.min(naive_errors_pct):.2e}%")
    print(f"  Max Error: {np.max(naive_errors_pct):.2e}%")
    
    # Compare performance
    adv_better = sum(1 for i in range(len(q_values)) if adv_errors_pct[i] < naive_errors_pct[i])
    print(f"\nADV performs better than Naive in {adv_better}/{len(q_values)} cases ({adv_better/len(q_values)*100:.1f}%)")

if __name__ == "__main__":
    create_comparison_plot()
