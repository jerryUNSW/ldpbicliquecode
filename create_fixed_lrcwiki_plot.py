import matplotlib.pyplot as plt
import numpy as np
import re
import os

def parse_fixed_results(filepath):
    q_values = []
    ground_truths = {}
    algorithm_errors = {
        "Naive": [],
        "oneR": [],
        "MRCN": [],
        "MRCN+": [],
        "MRCN++": []
    }
    
    current_section = None
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Extract ground truth values
            if "Ground Truth for (2,Q)-bicliques:" in line:
                # Look for the ground truth values in the output
                continue
            
            # Look for the ground truth values in the summary
            if "Ground Truth" in line and "Mean Estimate" in line:
                continue
                
            # Parse the algorithm comparison summary
            if "=== Algorithm Comparison Summary ===" in line:
                current_section = "Summary"
                continue
            
            if current_section == "Summary" and line.startswith("Q\t"):
                # Skip header line
                continue
                
            if current_section == "Summary" and line.startswith("-"):
                # Skip separator line
                continue
                
            if current_section == "Summary" and line and not line.startswith("="):
                # Parse the summary table
                parts = line.split('\t')
                if len(parts) >= 6:
                    try:
                        q = int(parts[0])
                        if q not in q_values:
                            q_values.append(q)
                        
                        # Extract relative errors for each algorithm
                        algorithm_errors["Naive"].append(float(parts[1]))
                        algorithm_errors["oneR"].append(float(parts[2]))
                        algorithm_errors["MRCN"].append(float(parts[3]))
                        algorithm_errors["MRCN+"].append(float(parts[4]))
                        algorithm_errors["MRCN++"].append(float(parts[5]))
                    except (ValueError, IndexError):
                        continue
    
    # Sort Q values
    q_values.sort()
    
    return q_values, ground_truths, algorithm_errors

def create_fixed_plot(dataset_name, epsilon, q_values, algorithm_errors, output_dir="."):
    # Filter out Q=10 if its ground truth is 0, as relative error is 1.0 by definition
    filtered_q_values = [q for q in q_values if q != 10]
    filtered_algorithm_errors = {
        algo: [err for i, err in enumerate(errors) if q_values[i] != 10]
        for algo, errors in algorithm_errors.items()
    }

    # Define colors and hatch patterns
    colors = {
        "Naive": "white",
        "oneR": "lightgrey", 
        "MRCN": "grey",
        "MRCN+": "black",
        "MRCN++": "darkgrey"
    }
    hatches = {
        "Naive": None,
        "oneR": None,
        "MRCN": None,
        "MRCN+": "//",
        "MRCN++": "x"
    }
    edgecolors = {
        "Naive": "black",
        "oneR": "black",
        "MRCN": "black",
        "MRCN+": "black",
        "MRCN++": "black"
    }

    # --- Main Comparison Plot ---
    fig, ax = plt.subplots(figsize=(12, 7))
    bar_width = 0.15
    index = np.arange(len(filtered_q_values))

    for i, (algo, errors) in enumerate(filtered_algorithm_errors.items()):
        ax.bar(index + i * bar_width, errors, bar_width, label=algo,
               color=colors[algo], edgecolor=edgecolors[algo], hatch=hatches[algo], log=True)

    ax.set_xlabel("q", fontsize=16)
    ax.set_ylabel("Relative Error (log scale)", fontsize=16)
    ax.set_title(f"Algorithm Comparison for {dataset_name.capitalize()} (P=2, ε={epsilon}) - FIXED", fontsize=18)
    ax.set_xticks(index + 2 * bar_width)
    ax.set_xticklabels([f"{q}" for q in filtered_q_values], fontsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.legend(fontsize=20, frameon=False, ncol=2, columnspacing=0.5)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    output_filename = os.path.join(output_dir, f"algorithm_comparison_{dataset_name}_eps{epsilon}_FIXED.pdf")
    plt.savefig(output_filename)
    print(f"Fixed plot saved to: {output_filename}")
    plt.close()

    # --- ADV-only Plot (for better visibility of smaller errors) ---
    fig, ax = plt.subplots(figsize=(12, 7))
    adv_algorithms = ["oneR", "MRCN", "MRCN+", "MRCN++"]
    adv_errors_filtered = {algo: filtered_algorithm_errors[algo] for algo in adv_algorithms}

    bar_width_adv = 0.18
    index_adv = np.arange(len(filtered_q_values))

    for i, (algo, errors) in enumerate(adv_errors_filtered.items()):
        ax.bar(index_adv + i * bar_width_adv, errors, bar_width_adv, label=algo,
               color=colors[algo], edgecolor=edgecolors[algo], hatch=hatches[algo], log=True)

    ax.set_xlabel("q", fontsize=16)
    ax.set_ylabel("Relative Error (log scale)", fontsize=16)
    ax.set_title(f"MRCN Algorithm Comparison for {dataset_name.capitalize()} (P=2, ε={epsilon}) - FIXED", fontsize=18)
    ax.set_xticks(index_adv + 1.5 * bar_width_adv)
    ax.set_xticklabels([f"{q}" for q in filtered_q_values], fontsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.legend(fontsize=20, frameon=False, ncol=2, columnspacing=0.5)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    output_adv_only_filename = os.path.join(output_dir, f"algorithm_comparison_{dataset_name}_adv_only_eps{epsilon}_FIXED.pdf")
    plt.savefig(output_adv_only_filename)
    print(f"Fixed ADV-only plot saved to: {output_adv_only_filename}")
    plt.close()

    # Print comparison
    print(f"\n=== Fixed Results for {dataset_name.upper()} ===")
    print("Q\tNaive\t\toneR\t\tMRCN\t\tMRCN+\t\tMRCN++")
    print("-" * 70)
    for i, q in enumerate(filtered_q_values):
        print(f"{q}\t{algorithm_errors['Naive'][i]:.2e}\t\t{algorithm_errors['oneR'][i]:.2e}\t\t"
              f"{algorithm_errors['MRCN'][i]:.2e}\t\t{algorithm_errors['MRCN+'][i]:.2e}\t\t"
              f"{algorithm_errors['MRCN++'][i]:.2e}")

if __name__ == "__main__":
    dataset_name = "lrcwiki"
    epsilon = 1
    results_filepath = f"lrcwiki_fixed_results.txt"

    print(f"Parsing {dataset_name.upper()} fixed results...")
    q_values, ground_truths, algorithm_errors = parse_fixed_results(results_filepath)

    print("Creating fixed comparison plots...")
    create_fixed_plot(dataset_name, epsilon, q_values, algorithm_errors, output_dir=".")
    print("\nFixed LRCWiki plots created successfully!")
    print("Files generated:")
    print(f"- algorithm_comparison_{dataset_name}_eps{epsilon}_FIXED.pdf")
    print(f"- algorithm_comparison_{dataset_name}_adv_only_eps{epsilon}_FIXED.pdf")
