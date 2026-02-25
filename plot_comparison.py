"""
Compare algorithms: OneR, ADV, ADV++, CN (common-neighbor)
1. Per-dataset: MRE vs epsilon for all 4 algorithms
2. Summary table
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

RESULT_DIR = "cn_btf_results/comparison"
DATASETS = {
    "to": "TO (33×217, btf=12)",
    "co": "CO (74×257, btf=273)",
    "unicode": "Unicode (255×615, btf=1662)",
    "lrcwiki": "lrcwiki (210×8737, btf=17.4M)",
    "librec-filmtrust-ratings": "filmtrust (1509×2072, btf=99.3M)",
    "rmwiki": "rmwiki (1226×8088, btf=97.8M)",
    "mooc": "MOOC (97×7047, btf=2.97B)",
}
ALGORITHMS = ["OneR", "ADV", "ADVpp", "CN", "PBC"]
ALG_LABELS = {"OneR": "OneR (1-round)", "ADV": "ADV (multi-round)", "ADVpp": "ADV++ (multi-round)", "CN": "CN (new, 1-round)", "PBC": "PBC (Li et al.)"}
EPSILONS = [0.5, 1.0, 1.5, 2.0, 2.5]

def load_results(dataset, alg, eps):
    fname = os.path.join(RESULT_DIR, f"{dataset}_{alg}_eps{eps}.txt")
    if not os.path.exists(fname):
        return None, None
    estimates, real_val = [], None
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("REAL"):
                real_val = float(line.split()[1])
            elif line:
                try:
                    estimates.append(float(line))
                except:
                    pass
    if len(estimates) == 0:
        return None, None
    return np.array(estimates), real_val

# ============================================================
# Figure: Per-dataset MRE comparison (2x4 grid)
# ============================================================
fig, axes = plt.subplots(2, 4, figsize=(24, 10))
axes = axes.flatten()

markers = {"OneR": "o", "ADV": "s", "ADVpp": "^", "CN": "D", "PBC": "P"}
colors = {"OneR": "#1f77b4", "ADV": "#ff7f0e", "ADVpp": "#2ca02c", "CN": "#d62728", "PBC": "#9467bd"}
linestyles = {"OneR": "--", "ADV": "-.", "ADVpp": ":", "CN": "-", "PBC": (0, (3, 1, 1, 1))}

for idx, (ds_key, ds_label) in enumerate(DATASETS.items()):
    ax = axes[idx]
    
    for alg in ALGORITHMS:
        mre_list, eps_valid = [], []
        for eps in EPSILONS:
            estimates, real_val = load_results(ds_key, alg, eps)
            if estimates is None or real_val is None or real_val == 0:
                continue
            mre = np.mean(np.abs(estimates - real_val) / abs(real_val))
            mre_list.append(mre)
            eps_valid.append(eps)
        
        if mre_list:
            ax.plot(eps_valid, mre_list, marker=markers[alg], color=colors[alg],
                    linewidth=2, markersize=8, label=ALG_LABELS[alg], linestyle=linestyles[alg])
    
    ax.set_xlabel("ε", fontsize=11)
    ax.set_ylabel("Mean Relative Error", fontsize=11)
    ax.set_title(ds_label, fontsize=12)
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(EPSILONS)
    ax.legend(fontsize=8)

# Hide unused subplot(s)
for j in range(len(DATASETS), len(axes)):
    axes[j].set_visible(False)

plt.suptitle("Algorithm Comparison: MRE vs ε for (2,2)-Biclique Counting", fontsize=14, y=1.01)
plt.tight_layout()
outpath = os.path.join(RESULT_DIR, "comparison_mre.pdf")
fig.savefig(outpath, bbox_inches='tight', dpi=150)
fig.savefig(outpath.replace('.pdf', '.png'), bbox_inches='tight', dpi=100)
print(f"Saved {outpath}")
plt.close(fig)

# ============================================================
# Summary table
# ============================================================
print("\n=== MRE Comparison Table ===")
header = f"{'Dataset':<25} {'ε':>4}"
for alg in ALGORITHMS:
    header += f" {ALG_LABELS[alg]:>20}"
print(header)
print("-" * (25 + 5 + 21 * len(ALGORITHMS)))

for ds_key, ds_label in DATASETS.items():
    short_label = ds_label.split("(")[0].strip()
    for eps in EPSILONS:
        row = f"{short_label:<25} {eps:>4.1f}"
        for alg in ALGORITHMS:
            estimates, real_val = load_results(ds_key, alg, eps)
            if estimates is None or real_val is None or real_val == 0:
                row += f" {'N/A':>20}"
            else:
                mre = np.mean(np.abs(estimates - real_val) / abs(real_val))
                row += f" {mre:>20.4f}"
        print(row)
    print()
