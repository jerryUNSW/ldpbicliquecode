#!/usr/bin/env python3
"""
EDUCATIONAL VERSION: Plot privacy budget allocation results in the style of Figure 6.
Creates publication-quality plots with log scale and markers.

This script demonstrates:
1. Professional matplotlib styling for publication-quality figures
2. Reading and parsing CSV data files
3. Creating line plots with multiple algorithms
4. Log scale handling for wide data ranges
5. Dynamic subplot layout for multiple datasets
6. Exporting high-quality PDF figures

Author: Enhanced for educational purposes
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# ============================================================================
# SECTION 1: GLOBAL PLOT STYLING CONFIGURATION
# ============================================================================
# This section sets up publication-quality plot parameters that apply to all
# plots created by matplotlib. These settings ensure consistent, professional
# appearance suitable for academic papers and presentations.

# DPI (Dots Per Inch) settings:
# - savefig.dpi: Resolution when saving figures (300 is publication quality)
# - figure.dpi: Resolution for on-screen display
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 300

# Axis styling:
# - axes.linewidth: Thickness of plot borders (1.5 is standard for publications)
plt.rcParams['axes.linewidth'] = 1.5

# Legend styling:
# - framealpha: Transparency of legend background (0 = transparent)
# - handletextpad: Space between legend marker and text
# - columnspacing: Space between columns in multi-column legends
plt.rcParams["legend.framealpha"] = 0
plt.rcParams["legend.handletextpad"] = 0.1
plt.rcParams["legend.columnspacing"] = 0.2

# PDF export settings:
# - pdf.fonttype: 42 ensures fonts are embedded as TrueType (required for LaTeX)
plt.rcParams['pdf.fonttype'] = 42

# Line and marker styling:
# - lines.linewidth: Thickness of plot lines
# - lines.markersize: Size of markers on data points
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 8

# ============================================================================
# SECTION 2: PLOTTING FUNCTION
# ============================================================================

def plot_dataset_multi_alg(csv_files_dict, dataset_name, ax, epsilon=1.0, p=2, q=2):
    """
    Plot multiple algorithms for a single dataset on the given axis.
    
    This function demonstrates:
    - Reading CSV files with pandas
    - Creating line plots with custom markers
    - Using log scale for wide data ranges
    - Customizing axis labels and legends
    
    Args:
        csv_files_dict: Dictionary mapping algorithm names to CSV file paths
                       Example: {'MRCN': 'path/to/mrcn.csv', 'MRCN+': 'path/to/mrcn+.csv'}
        dataset_name: Name of the dataset for title (e.g., 'CO', 'UN')
        ax: Matplotlib axis object where the plot will be drawn
        epsilon: Total privacy budget (used in title)
        p: First parameter of (p,q)-biclique (used in title)
        q: Second parameter of (p,q)-biclique (used in title)
    
    Returns:
        bool: True if plot was created successfully, False otherwise
    """
    
    # Define algorithm-specific visual styles
    # Each algorithm gets a unique color and marker for easy distinction
    # Style dictionary structure: {algorithm_name: {color, marker, linestyle}}
    alg_styles = {
        'MRCN': {'color': 'blue', 'marker': 'o', 'linestyle': '-'},      # Blue circles
        'MRCN+': {'color': 'green', 'marker': 's', 'linestyle': '-'},     # Green squares
        'MRCN++': {'color': 'red', 'marker': '^', 'linestyle': '-'}        # Red triangles
    }
    
    # Track whether we successfully loaded any data
    has_data = False
    alpha1_values = None  # Will store x-axis values from the first valid dataset
    
    # Iterate through each algorithm's CSV file
    # sorted() ensures consistent ordering across runs
    for alg_name, csv_file in sorted(csv_files_dict.items()):
        # Read CSV file using pandas
        # pandas automatically handles headers and converts to DataFrame
        df = pd.read_csv(csv_file)
        
        # Check if the CSV file has data
        if len(df) == 0:
            print(f"Warning: No data for {dataset_name} - {alg_name}")
            continue
        
        # Mark that we have at least one dataset with data
        has_data = True
        
        # Extract data columns
        # .values converts pandas Series to numpy array for plotting
        alpha1_values = df['alpha1'].values      # X-axis: budget allocation parameter
        rel_errors = df['relative_error'].values  # Y-axis: relative error values
        
        # Get the style for this algorithm (default to black if not found)
        style = alg_styles.get(alg_name, {'color': 'black', 'marker': 'o', 'linestyle': '-'})
        
        # Create the line plot with custom styling
        # Key parameters explained:
        # - alpha1_values, rel_errors: X and Y data
        # - marker: Shape of data point markers ('o'=circle, 's'=square, '^'=triangle)
        # - markerfacecolor: Fill color of markers ('white' creates hollow markers)
        # - markeredgecolor: Border color of markers (matches line color)
        # - markeredgewidth: Thickness of marker borders
        # - color: Color of the line connecting markers
        # - linewidth: Thickness of the line
        # - linestyle: Style of line ('-' = solid, '--' = dashed, etc.)
        # - label: Text for legend
        ax.plot(alpha1_values, rel_errors, 
                marker=style['marker'], 
                markerfacecolor='white',           # Hollow markers (white fill)
                markeredgecolor=style['color'],     # Colored borders
                markeredgewidth=1.5, 
                color=style['color'],               # Line color
                linewidth=2,
                linestyle=style['linestyle'],
                label=alg_name)
    
    # If no data was loaded, return False to indicate failure
    if not has_data:
        return False
    
    # ========================================================================
    # Y-AXIS: LOG SCALE SETUP
    # ========================================================================
    # Log scale is essential when data spans multiple orders of magnitude
    # (e.g., from 0.001 to 1000). It compresses large values and expands
    # small values, making trends visible across the entire range.
    ax.set_yscale('log')
    
    # ========================================================================
    # AXIS LABELS AND TITLE
    # ========================================================================
    # Use LaTeX math notation for Greek letters and symbols
    # r'...' creates a raw string (backslashes are literal, not escape chars)
    # $...$ enters math mode in LaTeX
    ax.set_xlabel(r'$\alpha_1$', fontsize=20)  # Greek letter alpha with subscript
    ax.set_ylabel('relative error', fontsize=20)
    
    # Title includes dataset name and experimental parameters
    # $\\epsilon$ renders as ε (epsilon) in LaTeX
    ax.set_title(f'{dataset_name}, $\\epsilon = {epsilon}$, $p = {p}$, $q = {q}$', fontsize=16)
    
    # ========================================================================
    # LEGEND
    # ========================================================================
    # loc='best' automatically places legend where it least obscures data
    # ncol=1 means single column (vertical layout)
    ax.legend(loc='best', fontsize=14, ncol=1)
    
    # ========================================================================
    # GRID
    # ========================================================================
    # Grid helps readers estimate values between tick marks
    # alpha=0.3 makes grid semi-transparent (30% opacity)
    # linestyle='--' creates dashed lines
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    
    # ========================================================================
    # X-AXIS TICK FORMATTING
    # ========================================================================
    # Customize x-axis to show alpha1 values clearly
    if alpha1_values is not None:
        # Set exact positions of ticks (one per data point)
        ax.set_xticks(alpha1_values)
        
        # Format tick labels: show one decimal place
        # f'{v:.1f}' formats number with 1 decimal place
        ax.set_xticklabels([f'{v:.1f}' for v in alpha1_values], 
                          rotation=0,      # No rotation (horizontal text)
                          fontsize=20)
        
        # Set y-axis tick label size
        ax.tick_params(axis='y', labelsize=20)
    
    return True


# ============================================================================
# SECTION 3: MAIN FUNCTION - DATA PROCESSING AND PLOT GENERATION
# ============================================================================

def main():
    """
    Main function that orchestrates the entire plotting process.
    
    Workflow:
    1. Parse command-line arguments for input directory and parameters
    2. Find and parse CSV files from the results directory
    3. Group files by dataset and algorithm
    4. Filter to only complete datasets (all algorithms present)
    5. Generate individual plots for each dataset
    6. Generate a combined comparison plot with all datasets
    """
    
    # ========================================================================
    # COMMAND-LINE ARGUMENT PARSING
    # ========================================================================
    # sys.argv[0] is the script name, sys.argv[1:] are user arguments
    # This allows flexible usage: python script.py [dir] [epsilon] [p] [q]
    
    # Get directory from command line or use default
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = 'budget_allocation_results'
    
    # Get epsilon, P and Q values from command line (optional, with defaults)
    # float() and int() convert string arguments to numbers
    epsilon = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    p = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    q = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    
    # Validate that the directory exists
    if not os.path.exists(results_dir):
        print(f"Error: Directory {results_dir} not found")
        return
    
    # ========================================================================
    # CSV FILE DISCOVERY
    # ========================================================================
    # Find all CSV files matching the expected pattern
    # List comprehension: [f for f in ... if condition] creates a filtered list
    csv_files = [f for f in os.listdir(results_dir) 
                 if f.endswith('_budget_allocation.csv')]
    
    if not csv_files:
        print("No CSV files found")
        return
    
    # ========================================================================
    # DATASET AND ALGORITHM NAME MAPPING
    # ========================================================================
    # These dictionaries provide human-readable names for display
    # Keys are internal identifiers, values are display names
    
    # Map internal dataset keys to display names
    dataset_names = {
        'co': 'CO',
        'to': 'TO', 
        'unicode': 'UN',
        'lrcwiki': 'LR',
        'csbwiki': 'CS',
        'librec-filmtrust-ratings': 'LIBREC-FILMTRUST-RATINGS',
        'rmwiki': 'RM',
        'bag-kos': 'BAG-KOS',
        'bpywiki': 'BP'
    }
    
    # Map internal algorithm keys to display names
    algorithm_names = {
        'alg2': 'MRCN',
        'alg3': 'MRCN+',
        'alg4': 'MRCN++'
    }
    
    # ========================================================================
    # CSV FILE PARSING AND GROUPING
    # ========================================================================
    # Group CSV files by dataset and algorithm
    # Structure: datasets_dict[dataset_key][alg_name] = csv_path
    
    datasets_dict = {}  # Will store: dataset_key -> {alg_name: csv_path}
    
    for csv_file in sorted(csv_files):
        # Parse filename to extract dataset and algorithm
        # Expected format: "dataset_alg#_budget_allocation.csv"
        # Example: "co_alg2_budget_allocation.csv" -> dataset='co', alg='alg2'
        
        if '_alg' in csv_file:
            # Split on '_alg' to separate dataset from algorithm number
            parts = csv_file.replace('_budget_allocation.csv', '').split('_alg')
            dataset_key = parts[0]      # Everything before '_alg'
            alg_key = 'alg' + parts[1]   # 'alg' + number (e.g., 'alg2')
        else:
            # Handle old format without '_alg' in filename
            dataset_key = csv_file.replace('_budget_allocation.csv', '')
            alg_key = 'alg4'  # Default to MRCN++ for old format
        
        csv_path = os.path.join(results_dir, csv_file)
        
        # Validate that CSV file has data
        df = pd.read_csv(csv_path)
        if len(df) == 0:
            print(f"Skipping {dataset_key} - {alg_key}: no data")
            continue
        
        # Convert internal algorithm key to display name
        alg_name = algorithm_names.get(alg_key, alg_key.upper())
        
        # Add to datasets dictionary
        # Create nested dictionary structure if needed
        if dataset_key not in datasets_dict:
            datasets_dict[dataset_key] = {}
        datasets_dict[dataset_key][alg_name] = csv_path
    
    # ========================================================================
    # FILTER TO COMPLETE DATASETS
    # ========================================================================
    # Only plot datasets that have results for all expected algorithms
    # This ensures fair comparison (all algorithms tested on same datasets)
    
    complete_datasets = {}
    expected_algorithms = {'MRCN', 'MRCN+', 'MRCN++'}  # Set of required algorithms
    
    # Optional: exclude datasets that are incomplete or in progress
    exclude_datasets = set()  # Empty set means no exclusions
    
    for dataset_key, alg_dict in datasets_dict.items():
        # Skip excluded datasets
        if dataset_key in exclude_datasets:
            print(f"Skipping {dataset_key}: excluded (in progress or not verified)")
            continue
        
        # Check if all expected algorithms are present
        available_algs = set(alg_dict.keys())
        if available_algs == expected_algorithms:
            # All algorithms present - add to complete datasets
            complete_datasets[dataset_key] = alg_dict
        else:
            # Missing some algorithms - report and skip
            missing = expected_algorithms - available_algs
            print(f"Skipping {dataset_key}: missing algorithms {missing}")
    
    # Check if we have any complete datasets to plot
    if not complete_datasets:
        print("No complete datasets to plot (need all 3 algorithms per dataset)")
        return
    
    datasets_dict = complete_datasets  # Use only complete datasets
    
    # ========================================================================
    # GENERATE INDIVIDUAL PLOTS
    # ========================================================================
    # Create one plot per dataset, with all algorithms on the same plot
    # This allows easy comparison of algorithms for each dataset
    
    print(f"\nGenerating individual plots for {len(datasets_dict)} datasets...")
    for dataset_key, alg_csv_dict in sorted(datasets_dict.items()):
        # Get display name for dataset
        dataset_name = dataset_names.get(dataset_key, dataset_key.upper())
        
        # Create a new figure with single subplot
        # figsize=(6, 5) sets width=6 inches, height=5 inches
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        
        # Plot all algorithms for this dataset
        success = plot_dataset_multi_alg(alg_csv_dict, dataset_name, ax, 
                                        epsilon=epsilon, p=p, q=q)
        
        if success:
            # Adjust layout to prevent label clipping
            plt.tight_layout()
            
            # Generate output filename
            # Format epsilon: integer if whole number, otherwise decimal
            # Replace '.' with 'p' to avoid filesystem issues
            eps_str = f"{int(epsilon)}" if epsilon == int(epsilon) else f"{epsilon:.1f}".replace('.', 'p')
            output_pdf = os.path.join(results_dir, 
                                     f'{dataset_key}-budget-allocation-eps{eps_str}-p{p}q{q}.pdf')
            
            # Save as high-quality PDF
            # dpi=300: High resolution for publication
            # bbox_inches='tight': Remove extra whitespace around plot
            plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
            plt.close()  # Free memory by closing figure
            
            # Calculate and print statistics for each algorithm
            print(f"\n{dataset_name}:")
            for alg_name, csv_path in sorted(alg_csv_dict.items()):
                df = pd.read_csv(csv_path)
                min_error = df['relative_error'].min()  # Best (minimum) error
                best_alpha1 = df.loc[df['relative_error'].idxmin(), 'alpha1']  # Best parameter value
                print(f"  {alg_name}: Best α₁ = {best_alpha1:.2f}, RelErr = {min_error:.4f}")
            print(f"  Saved: {output_pdf}")
    
    # ========================================================================
    # GENERATE COMBINED COMPARISON PLOT
    # ========================================================================
    # Create a grid of subplots showing all datasets side-by-side
    # This allows comparison across datasets at a glance
    
    n_datasets = len(datasets_dict)
    
    # Determine grid layout based on number of datasets
    # Goal: Create a roughly square grid that fits all datasets
    if n_datasets == 1:
        n_rows, n_cols = 1, 1
        figsize = (6, 5)
    elif n_datasets == 2:
        n_rows, n_cols = 1, 2
        figsize = (12, 5)
    elif n_datasets <= 4:
        n_rows, n_cols = 2, 2
        figsize = (12, 10)
    elif n_datasets <= 6:
        n_rows, n_cols = 2, 3
        figsize = (18, 10)
    else:
        n_rows, n_cols = 3, 3
        figsize = (18, 15)
    
    # Create figure with subplot grid
    # subplots() returns figure object and array of axis objects
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    
    # Handle single subplot case (axes is not an array)
    if n_datasets == 1:
        axes = [axes]
    else:
        # Flatten 2D array of axes to 1D for easier iteration
        axes = axes.flatten()
    
    print(f"\nGenerating combined comparison plot...")
    
    # Plot each dataset in its own subplot
    for idx, (dataset_key, alg_csv_dict) in enumerate(sorted(datasets_dict.items())):
        if idx < len(axes):
            dataset_name = dataset_names.get(dataset_key, dataset_key.upper())
            plot_dataset_multi_alg(alg_csv_dict, dataset_name, axes[idx], 
                                 epsilon=epsilon, p=p, q=q)
    
    # Hide unused subplots (if grid is larger than number of datasets)
    for idx in range(len(datasets_dict), len(axes)):
        axes[idx].axis('off')  # Make subplot invisible
    
    # Adjust layout
    plt.tight_layout()
    
    # Save combined plot
    eps_str = f"{int(epsilon)}" if epsilon == int(epsilon) else f"{epsilon:.1f}".replace('.', 'p')
    output_combined_pdf = os.path.join(results_dir, 
                                      f'budget-allocation-comparison-eps{eps_str}-p{p}q{q}.pdf')
    
    plt.savefig(output_combined_pdf, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nCombined plot saved: {output_combined_pdf}")
    print("\n✓ All plots generated successfully!")


# ============================================================================
# SECTION 4: SCRIPT ENTRY POINT
# ============================================================================

if __name__ == '__main__':
    # This block runs only when script is executed directly (not imported)
    # This is a Python best practice for reusable code
    main()
