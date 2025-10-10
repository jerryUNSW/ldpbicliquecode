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
    "wedge", 
    "wedge-avg"
]

# epsilons = [round(0.2 * i, 1) for i in range(1, 13)]  # 0.2 to 2.4
epsilons = [round(0.2 * i, 1) for i in range(2, 13)]  # 0.4 to 2.4


naive_label = "Naive"
dbe_label = "OneR"
adv_label = "MRCN"
adv_plus_label = "MRCN++"


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

def collect_data():
    data = {}  
    # for dataset in os.listdir(base_dir):

    for dataset in [
            "unicode", 
            "lrcwiki",
            "librec-filmtrust-ratings",
            "rmwiki",
            "csbwiki", 
            "bag-kos",
            "bpywiki",
            "nips",
            # "lastfm_band"
        ]:

        dataset_dir = os.path.join(base_dir, dataset)
        if not os.path.isdir(dataset_dir):
            continue

        pq_pairs = [(2, 3), (2, 4), (2, 5)]

        for P, Q in pq_pairs:
            key = (P, Q)
            data.setdefault(key, {})
            data[key].setdefault(dataset, {m: [] for m in methods})

            for eps in epsilons:
                for method in methods:
                    fname = f"{dataset}-{method}-eps{eps:.1f}-p{P}q{Q}.txt"
                    # print("fname = ", fname)
                    fpath = os.path.join(dataset_dir, fname)
                    # print("fpath = ", fpath)
                    err = extract_error(fpath)

                    # print("err = ", err)
                    data[key][dataset][method].append(err if err is not None else float('nan'))

    # print("data = ", data.shape)

    # exit(1)

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
            dbe_vals = [errors["dbe"][i] for i in eps_indices]
            wedge_base_vals = [errors["wedge-base"][i] for i in eps_indices]
            wedge_vals = [errors["wedge"][i] for i in eps_indices]
            wedge_avg_vals = [errors["wedge-avg"][i] for i in eps_indices]
            
            # plt.figure(figsize=(10, 6))
            
            # Plot all five methods with offset positions, edge colors, and hatches
            positions = [x + i * width - (n * width / 2) for i in range(n)]

            plt.bar(positions[0], naive_vals, width=width, 
                label=naive_label, linewidth=0.5, edgecolor='black', fc='white')

            plt.bar(positions[1], dbe_vals, 
                width=width, label=dbe_label, linewidth=0.5, hatch="", edgecolor='black', fc='lightgrey')
            
            plt.bar(positions[2], wedge_base_vals, 
                width=width, label=adv_label, linewidth=0.5, edgecolor='black', fc='grey')

            # plt.bar(positions[3], wedge_vals, width=width, label='wedge', linewidth=0.5, hatch="//", edgecolor='black', fc='darkgrey')
            
            plt.bar(positions[3], wedge_avg_vals, width=width, label=adv_plus_label, 
                linewidth=0.5, hatch="////", edgecolor='black', fc='black')

            plt.xticks(
                x-width, 
                [str(eps) if i % 2 == 0 else '' for i, eps in enumerate(filtered_epsilons)], fontsize=18)

            plt.xlabel(r"$\varepsilon$", fontsize=20)

            plt.ylabel("Mean Relative Error", fontsize=20)
            plt.yscale('log')
            
            # Dynamic y-axis range and ticks
            all_errors = naive_vals + dbe_vals + wedge_base_vals + wedge_vals + wedge_avg_vals
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
            plt.title(f"(P={P}, Q={Q})", fontsize=20)

            # plt.legend(fontsize=14, ncol=2, loc="upper right", columnspacing=0.5)

            plt.legend(fontsize=14, ncol=4, loc="upper center", columnspacing=0.1)
            # plt.grid(axis='y', alpha=0.3, which='both')
            
            plt.tight_layout()
            outfile = f"plots/barplots/{dataset}-p{P}q{Q}-bar.pdf"
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved: {outfile}")


def extract_com(filepath):
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        for line in reversed(lines):
            if "com cost (MB)" in line:
                match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
                if match:
                    mb = float(match.group())
                    return int(mb * 1024 * 1024)  # Convert MB to bytes
    except Exception as e:
        print(f"Failed to read {filepath}: {e}")
    return None



def collect_com_data():
    data = {}  
    # for dataset in os.listdir(base_dir):

    for dataset in [
            "unicode", 
            "lrcwiki",
            "librec-filmtrust-ratings",
            "rmwiki",
            "csbwiki", 
            "bag-kos"
        ]:

        dataset_dir = os.path.join(base_dir, dataset)
        if not os.path.isdir(dataset_dir):
            continue

        methods__ = [
                        "naive", 
                        "dbe",
                        "wedge-base", 
                        "wedge-avg"
                    ]

        for P, Q in [(2,2), (2,3)]:
            key = (P, Q)
            data.setdefault(key, {})
            data[key].setdefault(dataset, {m: [] for m in methods__})

            for eps in [round(1.6 - 0.2 * i, 1) for i in reversed(range(7))]:

                for method in methods__:

                    # fname = f"{dataset}-{method}-eps{eps:.1f}-p{P}q{Q}-com.txt"
                    fname = f"{dataset}-{method}-eps{eps:.1f}-p{P}q{Q}-com-new.txt"
                    # print("fname = ", fname)

                    if method == "dbe":
                        # method = "naive"
                        # fname = f"{dataset}-naive-eps{eps:.1f}-p{P}q{Q}-com.txt"
                        fname = f"{dataset}-naive-eps{eps:.1f}-p{P}q{Q}-com-new.txt"

                    fpath = os.path.join(dataset_dir, fname)
                    # print("fpath = ", fpath)
                    err = extract_com(fpath)

                    # print("err = ", err)
                    data[key][dataset][method].append(err if err is not None else float('nan'))

    return data


def plot_comm_cost_per_dataset(com_data):

    plt_settings()

    os.makedirs("plots/barplots_comm", exist_ok=True)

    methods_to_plot = ["naive", "dbe", "wedge-base", "wedge-avg"]
    method_styles = {
        "naive":      {"color": "white",     "hatch": "",    "label": naive_label},
        # "dbe":        {"color": "lightgray", "hatch": "//",  "label": "OneR"},
        # "wedge-base": {"color": "gray",      "hatch": "xx",  "label": "MultiR-SS"},
        "dbe":        {"color": "lightgray", "hatch": "",  "label": dbe_label},
        "wedge-base": {"color": "gray",      "hatch": "",  "label": adv_label},
        "wedge-avg":  {"color": "black",     "hatch": "",    "label": adv_plus_label},
    }

    eps_range = [round(1.6 - 0.2 * i, 1) for i in reversed(range(7))]
    x = np.arange(len(eps_range))
    total_width, n = 0.9, len(methods_to_plot)
    width = total_width / (n + 1)

    for (P, Q), datasets in com_data.items():
        for dataset_idx, (dataset, comcosts) in enumerate(datasets.items()):

            for idx, method in enumerate(methods_to_plot):
                values_mb = comcosts[method]
                positions = x + idx * width - (n * width / 2)
                style = method_styles[method]

                plt.bar(
                    positions, values_mb, width=width,
                    color=style["color"], hatch=style["hatch"],
                    edgecolor='black', linewidth=0.5,
                    label=style["label"]
                )

            plt.xticks(x, [str(eps) for eps in eps_range], fontsize=14)
            plt.xlabel(r"$\varepsilon$", fontsize=16)
            plt.ylabel("Communication Cost (MB)", fontsize=16)
            plt.yscale('log')

            all_vals = [v for m in methods_to_plot for v in comcosts[m] if v > 0]
            y_min = min(all_vals)
            y_max = max(all_vals)
            plt.ylim(y_min / 10, y_max * 10)

            indices = list(range(
                math.floor(math.log10(y_min / 10)),
                math.ceil(math.log10(y_max * 10)) + 1,
                2
            ))
            pos, labels = get_pos_and_labels(indices)
            plt.yticks(pos, labels, fontsize=14)

            # plt.title(f"({chr(97 + dataset_idx)}) {dataset} (P={P}, Q={Q})", fontsize=16)
            plt.legend(fontsize=14, ncol=4, loc='upper right')
            # plt.grid(axis='y', linestyle='--', alpha=0.4)

            plt.tight_layout()
            outfile = f"plots/barplots_comm/{dataset}-P{P}Q{Q}-commcost.pdf"
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved: {outfile}")



if __name__ == "__main__":


    # evaluate the effect of privacy budget 
    # data = collect_data()
    # # Print the collected data in table format
    # for (P, Q), datasets in data.items():
    #     print(f"\n=== P={P}, Q={Q} ===")
    #     for dataset, errors in datasets.items():
    #         print(f"\nDataset: {dataset}")
    #         df = pd.DataFrame(errors, index=epsilons)
    #         df.index.name = 'Epsilon'
    #         print(df.to_markdown(tablefmt="grid", floatfmt=".4f"))
    # plot_bar_comparison(data)

    ## the communication cost bar plots 
    com_data = collect_com_data()
    for (P, Q), datasets in com_data.items():
        print(f"\n=== P={P}, Q={Q} ===")
        for dataset, comcosts in datasets.items():
            print(f"\nDataset: {dataset}")
            df = pd.DataFrame(comcosts, index=[round(1.6 - 0.2 * i, 1) for i in reversed(range(7))])
            # df = pd.DataFrame(comcosts)
            df.index.name = 'Epsilon'
            print(df.to_markdown(tablefmt="grid", floatfmt=".4f"))
    plot_comm_cost_per_dataset(com_data)

    subprocess.run(["bash", "transfer.sh"])