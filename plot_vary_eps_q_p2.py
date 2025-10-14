#!/usr/bin/env python3
import os
import re
import sys
from collections import defaultdict
import matplotlib.pyplot as plt


def parse_log(log_path):
    estimate = None
    real = None
    rel_err = None
    try:
        with open(log_path, 'r') as f:
            for line in f:
                if 'relative error' in line:
                    # e.g., relative error = 13.2144291293
                    try:
                        rel_err = float(re.findall(r"relative error\s*=\s*([\-0-9\.eE]+)", line)[0])
                    except Exception:
                        pass
                elif line.strip().startswith('estimate ='):
                    try:
                        estimate = float(line.split('=')[1].strip())
                    except Exception:
                        pass
                elif line.strip().startswith('real ='):
                    try:
                        real = float(line.split('=')[1].strip())
                    except Exception:
                        pass
    except FileNotFoundError:
        return None
    return {
        'estimate': estimate,
        'real': real,
        'rel_err': rel_err,
    }


def main():
    if len(sys.argv) < 4:
        print("Usage: plot_vary_eps_q_p2.py <results_root> <dataset> <eps>")
        sys.exit(1)

    results_root = sys.argv[1]
    dataset = sys.argv[2]
    eps = sys.argv[3]

    base_dir = os.path.join(results_root, dataset, f"eps{eps}")
    if not os.path.isdir(base_dir):
        print(f"Missing directory: {base_dir}")
        sys.exit(1)

    # Collect logs by q and algorithm
    # filename format: p2_q<q>_algo<algo>.log
    data = defaultdict(dict)  # data[q][algo] = metrics
    algo_labels = {0: 'Naive', 2: 'Algo2', 3: 'Algo3', 4: 'Algo4'}

    for name in os.listdir(base_dir):
        if not name.startswith('p2_q') or not name.endswith('.log'):
            continue
        m = re.match(r"p2_q(\d+)_algo(\d+)\.log", name)
        if not m:
            continue
        q = int(m.group(1))
        algo = int(m.group(2))
        log_path = os.path.join(base_dir, name)
        metrics = parse_log(log_path)
        if metrics is None:
            continue
        data[q][algo] = metrics

    if not data:
        print("No logs parsed")
        sys.exit(1)

    qs = sorted(data.keys())
    algos = [0, 2, 3, 4]

    # Plot relative error (if available)
    rel_err_matrix = []
    for algo in algos:
        rel_err_matrix.append([data[q].get(algo, {}).get('rel_err', None) for q in qs])

    # Bar plot
    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.18
    x = range(len(qs))
    for i, algo in enumerate(algos):
        vals = rel_err_matrix[i]
        # replace missing/invalid with a tiny positive to support log scale
        vals_plot = []
        for v in vals:
            if v is None or (isinstance(v, float) and (v != v)):
                vals_plot.append(1e-12)
            else:
                vals_plot.append(max(v, 1e-12))
        ax.bar([xi + (i-1.5)*width for xi in x], vals_plot, width=width, label=algo_labels.get(algo, f"Algo{algo}"))

    ax.set_title(f"Relative Error by Algorithm (dataset={dataset}, eps={eps}, p=2)")
    ax.set_xlabel("q")
    ax.set_ylabel("Relative Error (log)")
    ax.set_yscale('log')
    ax.set_xticks(list(x))
    ax.set_xticklabels([str(q) for q in qs])
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    out_png = os.path.join(base_dir, f"p2_q4to8_eps{eps}_relerr.png")
    out_pdf = os.path.join(base_dir, f"p2_q4to8_eps{eps}_relerr.pdf")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.savefig(out_pdf)
    print(f"Saved plots: {out_png}, {out_pdf}")


if __name__ == '__main__':
    main()


