#!/usr/bin/env python3

import pandas as pd 
import numpy as np
import os
import re
import math
import matplotlib.pyplot as plt
import subprocess
import sympy as sp
from scipy.optimize import root_scalar
from scipy.optimize import minimize
from matplotlib.pyplot import cm
import sys
import requests
import seaborn as sns
import subprocess


# hyper paremers for plotting
plt.rcParams['savefig.dpi'] = 300 
plt.rcParams['figure.dpi'] = 300 

def plt_settings():
    plt.style.use('default')
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams["legend.framealpha"] = 0
    plt.rcParams["legend.handletextpad"] = 0.1
    plt.rcParams["legend.columnspacing"] = 0.2
    # for the varying plots 
    plt.rcParams["figure.figsize"] = (6,5)
    plt.rcParams['pdf.fonttype'] = 42


def get_pos_and_labels(indices):
    positions = []
    labels = []
    for i in indices:
        positions.append(pow(10, i))
        label = f"$10^{{{'-' if i < 0 else ''}{abs(i)}}}$"
        labels.append(label)
    return tuple(positions), tuple(labels)

base_dir = "biclique-results"
methods = [
    "naive", 
    "dbe",
    "wedge-base", 
    # "wedge", 
    "wedge-avg"
]

epsilons = [round(0.2 * i, 1) for i in range(2, 13)]  # 0.4 to 2.4

pq_pairs = [(2, 2), (2, 3)]

# pq_pairs = [(3, 2), (3, 3), (3, 4)]

def extract_error(filepath):
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        for line in reversed(lines):
            if "adv rel err" in line:
                match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
                if match:
                    return float(match.group())
    except Exception as e:
        print(f"Failed to read {filepath}: {e}")
    return None


def extract_time(filepath):
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        for line in reversed(lines):
            if "time:" in line:
                match = re.search(r"time:([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
                if match:
                    return float(match.group(1))
    except Exception as e:
        print(f"Failed to read {filepath}: {e}")
    return None


def count_rounds(filepath):
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        count = sum(1 for line in lines if "relative error" in line)
        return count
    except Exception as e:
        print(f"Failed to read {filepath}: {e}")
    return 0

def collect_data():
    data = {}  

    for dataset in [
            "unicode", 
            "lrcwiki",
            # "librec-filmtrust-ratings",
            # "rmwiki",
            # "csbwiki", 
            # "bag-kos",
            # "bpywiki",
            # "nips",
            # "lastfm_band"
        ]:

        dataset_dir = os.path.join(base_dir, dataset)
        if not os.path.isdir(dataset_dir):
            continue

        for P, Q in pq_pairs:
            key = (P, Q)
            data.setdefault(key, {})
            data[key].setdefault(dataset, {m: [] for m in methods})

            for eps in epsilons:
                for method in methods:

                    if P > 2 and method == "dbe":
                        continue

                    # fname = f"{dataset}-{method}-eps{eps:.1f}-p{P}q{Q}.txt"

                    fname = f"{dataset}-{method}-eps{eps:.1f}-p{P}q{Q}-time.txt"

                    # print("fname = ", fname)
                    fpath = os.path.join(dataset_dir, fname)
                    # print("fpath = ", fpath)
                    running_time = extract_time(fpath)

                    T = count_rounds(fpath)

                    time_per_round = running_time / T 

                    print(f"Total time: {running_time:.4f} sec")
                    print(f"Rounds: {T}")
                    print(f"Time per round: {time_per_round:.6f} sec")

                    # print("err = ", err)
                    data[key][dataset][method].append(time_per_round if time_per_round is not None else float('nan'))

    return data

def plot_data(data):

    plt_settings()

    os.makedirs("plots", exist_ok=True)
    for (P, Q), datasets in data.items():
        for dataset, errors in datasets.items():
            plt.figure(figsize=(8, 5))
            for method in methods:
                plt.plot(epsilons, errors[method], label=method, marker='o')
            plt.xlabel("Epsilon")
            plt.ylabel("Relative Error (%)")
            plt.yscale('log')
            plt.title(f"{dataset} (P={P}, Q={Q})")
            plt.legend()
            plt.grid(True)
            outfile = f"plots/{dataset}-p{P}q{Q}.png"
            plt.savefig(outfile)
            plt.close()
            print(f"Saved: {outfile}")

def plot_bar_comparison(data):
    plt_settings()

    # plt.grid(False)

    os.makedirs("plots", exist_ok=True)
    
    methods_to_plot = [
        "naive", "dbe", "wedge-base", "wedge", "wedge-avg"
    ]
    # eps_min, eps_max = 0.4, 1.6
    eps_min, eps_max = 0.4, 1.6
    filtered_epsilons = [eps for eps in epsilons if eps_min <= eps <= eps_max]
    eps_indices = [i for i, eps in enumerate(epsilons) if eps_min <= eps <= eps_max]
    
    for (P, Q), datasets in data.items():
        for dataset, errors in datasets.items():

            print(f"{dataset}\t {P} \t {Q}")
            x = np.arange(len(filtered_epsilons))
            total_width, n = 0.9, len(methods_to_plot)
            width = total_width / (n + 1)
            
            # Extract values for all methods
            naive_vals = [errors["naive"][i] for i in eps_indices]
            wedge_base_vals = [errors["wedge-base"][i] for i in eps_indices]
            wedge_vals = [errors["wedge"][i] for i in eps_indices]
            wedge_avg_vals = [errors["wedge-avg"][i] for i in eps_indices]
            
            # plt.figure(figsize=(10, 6))
            
            # Plot all five methods with offset positions, edge colors, and hatches
            positions = [x + i * width - (n * width / 2) for i in range(n)]

            plt.bar(positions[0], naive_vals, width=width, 
                label='Naive', linewidth=0.5, edgecolor='black', fc='white')

            if P == 2:
                dbe_vals = [errors["dbe"][i] for i in eps_indices]
                plt.bar(positions[1], dbe_vals, width=width, label='oneround', linewidth=0.5, edgecolor='black', fc='lightgrey')

            if P in {2, 3}:
                plt.bar(positions[1 if P == 3 else 2], wedge_base_vals, width=width, label='Adv', linewidth=0.5, edgecolor='black', fc='grey')
                plt.bar(positions[2 if P == 3 else 3], wedge_avg_vals, width=width, label='Adv++', linewidth=0.5, hatch="////", edgecolor='black', fc='black')

            plt.xticks(x-width, [str(eps) if i % 2 == 0 else '' for i, eps in enumerate(filtered_epsilons)], fontsize=18)

            plt.xlabel(r"$\varepsilon$", fontsize=20)

            plt.ylabel("Mean Relative Error", fontsize=20)
            plt.yscale('log')
            
            # Dynamic y-axis range and ticks
            all_errors = naive_vals + wedge_base_vals + wedge_vals + wedge_avg_vals

            if P == 2:
                all_errors = all_errors + dbe_vals

            y_min = np.min(all_errors)
            y_max = np.max(all_errors)
            y_range = y_max - y_min
            y_min_adjusted = y_min / 10 if y_min > 0 else y_min * 10
            y_max_adjusted = y_max * 10

            print("y_min_adjusted = ", y_min_adjusted)
            print("y_max_adjusted = ", y_max_adjusted)
            plt.ylim(y_min_adjusted, y_max_adjusted)
            
            indices = [i for i in range(math.floor(math.log10(y_min_adjusted)), math.ceil(math.log10(y_max_adjusted)) + 1, 2)]
            pos, labels = get_pos_and_labels(indices)
            plt.yticks(pos, labels, fontsize=20)
            
            # plt.title(f"{dataset} (P={P}, Q={Q})", fontsize=20)

            # plt.legend(fontsize=14, ncol=2, loc="upper right", columnspacing=0.5)

            plt.legend(fontsize=14, ncol=4, loc="upper center", columnspacing=0.1)
            # plt.grid(axis='y', alpha=0.3, which='both')
            
            plt.tight_layout()
            outfile = f"plots/barplots/{dataset}-p{P}q{Q}-bar.pdf"
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved: {outfile}")

def plot_bar_comparison_running_time(data):
    plt_settings()

    os.makedirs("plots", exist_ok=True)
    
    methods_to_plot = [
        "naive", "dbe", "wedge-base", "wedge", "wedge-avg"
    ]
    
    eps_min, eps_max = 0.4, 1.6
    filtered_epsilons = [eps for eps in epsilons if eps_min <= eps <= eps_max]
    eps_indices = [i for i, eps in enumerate(epsilons) if eps_min <= eps <= eps_max]
    
    for (P, Q), datasets in data.items():
        for dataset, times in datasets.items():

            print(f"{dataset}\t {P} \t {Q}")
            x = np.arange(len(filtered_epsilons))
            total_width, n = 0.9, len(methods_to_plot)
            width = total_width / (n + 1)
            
            # Extract running times
            naive_vals = [times["naive"][i] for i in eps_indices]


            wedge_base_vals = [times["wedge-base"][i] for i in eps_indices]
            # wedge_vals = [times["wedge"][i] for i in eps_indices]
            wedge_avg_vals = [times["wedge-avg"][i] for i in eps_indices]
            
            positions = [x + i * width - (n * width / 2) for i in range(n)]

            plt.bar(positions[0], naive_vals, width=width, label='Naive', linewidth=0.5, edgecolor='black', fc='white')

            if P == 2 :
                dbe_vals = [times["dbe"][i] for i in eps_indices]
                plt.bar(positions[1], dbe_vals, width=width, label='oneround', linewidth=0.5, edgecolor='black', fc='lightgrey')
                plt.bar(positions[2], wedge_base_vals, width=width, label='Adv', linewidth=0.5, edgecolor='black', fc='grey')
                plt.bar(positions[3], wedge_avg_vals, width=width, label='Adv++', linewidth=0.5, hatch="////", edgecolor='black', fc='black')
    
            if P == 3 : 
                plt.bar(positions[1], wedge_base_vals, width=width, label='Adv', linewidth=0.5, edgecolor='black', fc='grey')
                plt.bar(positions[2], wedge_avg_vals, width=width, label='Adv++', linewidth=0.5, hatch="////", edgecolor='black', fc='black')
                
            plt.xticks(
                x - width, 
                [str(eps) if i % 2 == 0 else '' for i, eps in enumerate(filtered_epsilons)], fontsize=18)

            plt.xlabel(r"$\varepsilon$", fontsize=20)
            plt.ylabel("Running Time (sec)", fontsize=20)
            plt.yscale('log')

            all_times = naive_vals + wedge_base_vals + wedge_avg_vals

            if P == 2 :
                all_times = all_times + dbe_vals
            y_max = max(all_times)
            y_min = min(all_times)

            # print("y_max = ", y_max)
            plt.ylim(y_min, y_max * 1.2)

            plt.yticks(fontsize=20)
            plt.legend(fontsize=14, ncol=4, loc="upper center", columnspacing=0.1)
            
            plt.tight_layout()
            outfile = f"plots/time-plot/{dataset}-p{P}q{Q}-runtime.pdf"
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved: {outfile}")

if __name__ == "__main__":
    data = collect_data()

    print("# Print the collected data in table format")
    for (P, Q), datasets in data.items():
        print(f"\n=== P={P}, Q={Q} ===")

        print("datasets.items() = ", datasets.items())
        # for dataset, errors in datasets.items():
        #     print(f"\nDataset: {dataset}")
        #     df = pd.DataFrame(errors, index=epsilons)
        #     df.index.name = 'Epsilon'
        #     print(df.to_markdown(tablefmt="grid", floatfmt=".4f"))
        for dataset, errors in datasets.items():
            print(f"\nDataset: {dataset}")
            
            # Filter out 'dbe' if P > 2
            filtered_errors = {k: v for k, v in errors.items() if not (k == "dbe" and P > 2)}

            # Pad each method's list to match length of epsilons
            for method in filtered_errors:
                vals = filtered_errors[method]
                if len(vals) < len(epsilons):
                    filtered_errors[method] = vals + [np.nan] * (len(epsilons) - len(vals))

            df = pd.DataFrame(filtered_errors, index=epsilons)
            df.index.name = 'Epsilon'
            print(df.to_markdown(tablefmt="grid", floatfmt=".4f"))


    # exit(0)

    # plot_data(data)

    plot_bar_comparison_running_time(data)

    subprocess.run(["bash", "transfer.sh"])




