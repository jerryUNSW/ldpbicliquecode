"""
Plot CN vs ADV++ comparison for (2,q)-biclique counting, q=2,3,4.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import glob

RESULT_DIR = "cn_btf_results/cn_2q_comparison"
DATASETS = ["co", "unicode", "lrcwiki", "mooc"]
DS_LABELS = {
    "co": "CO (74×257)",
    "unicode": "Unicode (255×615)",
    "lrcwiki": "lrcwiki (210×8737)",
    "mooc": "MOOC (97×7047)",
}
ALGORITHMS = ["CN", "ADVpp"]
ALG_LABELS = {"CN": "CN (new, 1-round)", "ADVpp": "ADV++ (multi-round)"}
ALG_COLORS = {"CN": "#d62728", "ADVpp": "#2ca02c"}
ALG_MARKERS = {"CN": "D", "ADVpp": "^"}
EPSILONS = [1.0, 1.5, 2.0, 2.5]
Q_VALUES = [2, 3, 4]

def load_results(dataset, alg, q, eps):
    fname = os.path.join(RESULT_DIR, f"{dataset}_{alg}_q{q}_eps{eps}.txt")
    if not os.path.exists(fname):
        return None, None
    estimates, real_val = [], None
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("REAL "):
                try:
                    real_val = float(line.split()[1])
                except:
                    pass
            else:
                try:
                    estimates.append(float(line))
                except:
                    pass
    if not estimates:
        return None, None
    return np.array(estimates), real_val

# ============================================================
# Figure: MRE comparison grid (rows=q, cols=dataset)
# ============================================================
fig, axes = plt.subplots(len(Q_VALUES), len(DATASETS), figsize=(5*len(DATASETS), 4*len(Q_VALUES)))

for qi, q in enumerate(Q_VALUES):
    for di, ds in enumerate(DATASETS):
        ax = axes[qi][di]
        for alg in ALGORITHMS:
            mre_list, eps_valid = [], []
            for eps in EPSILONS:
                estimates, real_val = load_results(ds, alg, q, eps)
                if estimates is None or real_val is None or real_val == 0:
                    continue
                mre = np.mean(np.abs(estimates - real_val) / abs(real_val))
                mre_list.append(mre)
                eps_valid.append(eps)
            if mre_list:
                ax.plot(eps_valid, mre_list, marker=ALG_MARKERS[alg], color=ALG_COLORS[alg],
                        linewidth=2, markersize=8, label=ALG_LABELS[alg])

        ax.set_xlabel("ε", fontsize=10)
        ax.set_ylabel("MRE", fontsize=10)
        ax.set_title(f"(2,{q})-biclique: {DS_LABELS[ds]}", fontsize=10)
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(EPSILONS)
        ax.legend(fontsize=8)

plt.suptitle("CN vs ADV++: MRE for (2,q)-Biclique Counting", fontsize=14, y=1.01)
plt.tight_layout()
outpath = os.path.join(RESULT_DIR, "cn_vs_advpp_2q_mre.pdf")
fig.savefig(outpath, bbox_inches='tight', dpi=150)
fig.savefig(outpath.replace('.pdf', '.png'), bbox_inches='tight', dpi=100)
print(f"Saved {outpath}")
plt.close(fig)

# ============================================================
# Summary table
# ============================================================
print("\n=== CN vs ADV++ MRE Table ===")
for q in Q_VALUES:
    print(f"\n--- (2,{q})-bicliques ---")
    header = f"{'Dataset':<15} {'ε':>4}"
    for alg in ALGORITHMS:
        header += f" {ALG_LABELS[alg]:>22}"
    print(header)
    print("-" * 70)

    for ds in DATASETS:
        for eps in EPSILONS:
            row = f"{ds:<15} {eps:>4.1f}"
            for alg in ALGORITHMS:
                estimates, real_val = load_results(ds, alg, q, eps)
                if estimates is None or real_val is None or real_val == 0:
                    row += f" {'N/A':>22}"
                else:
                    mre = np.mean(np.abs(estimates - real_val) / abs(real_val))
                    row += f" {mre:>22.4f}"
            print(row)
    print()
