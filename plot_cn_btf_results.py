"""
Plot distribution of common-neighbor butterfly estimates and mean relative error vs epsilon.
Two types of figures:
  1. Per-dataset distribution plots (one per dataset, subplots for each epsilon)
  2. Mean relative error vs epsilon across all datasets
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

RESULT_DIR = "cn_btf_results"
DATASETS = {
    "to": "TO (33×217)",
    "co": "CO (74×257)",
    "unicode": "Unicode (255×615)",
    "lrcwiki": "lrcwiki (210×8737)",
    "librec-filmtrust-ratings": "filmtrust (1509×2072)",
    "rmwiki": "rmwiki (1226×8088)",
}
EPSILONS = [0.5, 1.0, 1.5, 2.0, 2.5]

def load_results(dataset, eps):
    fname = os.path.join(RESULT_DIR, f"{dataset}_eps{eps}.txt")
    if not os.path.exists(fname):
        return None, None
    estimates = []
    real_val = None
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("REAL"):
                real_val = float(line.split()[1])
            elif line:
                try:
                    estimates.append(float(line))
                except ValueError:
                    pass
    return np.array(estimates), real_val

# ============================================================
# Figure 1: Distribution plots (one figure per dataset)
# ============================================================
for ds_key, ds_label in DATASETS.items():
    fig, axes = plt.subplots(1, 5, figsize=(22, 4), sharey=False)
    fig.suptitle(f"Estimate Distribution — {ds_label}", fontsize=14, y=1.02)
    
    for i, eps in enumerate(EPSILONS):
        ax = axes[i]
        estimates, real_val = load_results(ds_key, eps)
        if estimates is None or len(estimates) == 0:
            ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f"ε = {eps}")
            continue
        
        ax.hist(estimates, bins=30, alpha=0.7, color='steelblue', edgecolor='white', density=True)
        if real_val is not None:
            ax.axvline(x=real_val, color='red', linewidth=2, linestyle='--', label=f'True = {real_val:.1f}')
        mean_est = np.mean(estimates)
        ax.axvline(x=mean_est, color='orange', linewidth=2, linestyle='-', label=f'Mean = {mean_est:.1f}')
        ax.set_title(f"ε = {eps}", fontsize=12)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_xlabel("Estimate")
        if i == 0:
            ax.set_ylabel("Density")
    
    plt.tight_layout()
    outpath = os.path.join(RESULT_DIR, f"dist_{ds_key}.pdf")
    fig.savefig(outpath, bbox_inches='tight', dpi=150)
    print(f"Saved {outpath}")
    plt.close(fig)

# ============================================================
# Figure 2: Mean relative error vs epsilon
# ============================================================
fig, ax = plt.subplots(figsize=(8, 5))

markers = ['o', 's', '^', 'D', 'v', 'P']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

for idx, (ds_key, ds_label) in enumerate(DATASETS.items()):
    mre_list = []
    eps_valid = []
    for eps in EPSILONS:
        estimates, real_val = load_results(ds_key, eps)
        if estimates is None or real_val is None or real_val == 0 or len(estimates) == 0:
            continue
        rel_errors = np.abs(estimates - real_val) / abs(real_val)
        mre = np.mean(rel_errors)
        mre_list.append(mre)
        eps_valid.append(eps)
    
    if mre_list:
        ax.plot(eps_valid, mre_list, marker=markers[idx % len(markers)], 
                color=colors[idx % len(colors)], linewidth=2, markersize=8, label=ds_label)

ax.set_xlabel("ε (Privacy Budget)", fontsize=12)
ax.set_ylabel("Mean Relative Error", fontsize=12)
ax.set_title("Common-Neighbor Butterfly Estimator: MRE vs ε", fontsize=13)
ax.legend(fontsize=9)
ax.set_yscale('log')
ax.grid(True, alpha=0.3)
ax.set_xticks(EPSILONS)

plt.tight_layout()
outpath = os.path.join(RESULT_DIR, "mre_vs_epsilon.pdf")
fig.savefig(outpath, bbox_inches='tight', dpi=150)
print(f"Saved {outpath}")
plt.close(fig)

# ============================================================
# Print summary table
# ============================================================
print("\n=== Summary Table ===")
print(f"{'Dataset':<30} {'ε':>5} {'True':>15} {'Mean':>15} {'MRE':>10} {'Rounds':>7}")
print("-" * 90)
for ds_key, ds_label in DATASETS.items():
    for eps in EPSILONS:
        estimates, real_val = load_results(ds_key, eps)
        if estimates is None or len(estimates) == 0:
            print(f"{ds_label:<30} {eps:>5.1f} {'N/A':>15} {'N/A':>15} {'N/A':>10} {'N/A':>7}")
            continue
        mean_est = np.mean(estimates)
        if real_val and real_val != 0:
            mre = np.mean(np.abs(estimates - real_val) / abs(real_val))
        else:
            mre = float('nan')
        print(f"{ds_label:<30} {eps:>5.1f} {real_val:>15.1f} {mean_est:>15.1f} {mre:>10.4f} {len(estimates):>7}")
