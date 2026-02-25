"""
Analyze CN normality and higher-q moment correction results.
1. Normality tests (Shapiro-Wilk) + Q-Q plots for S'(u,v)
2. Unbiasedness verification for C(S,q) estimators, q=2,3,4
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as sp_stats
import os
import glob

RESULT_DIR = "cn_btf_results"

def load_samples(csv_path):
    """Load per-pair S' samples."""
    df = pd.read_csv(csv_path)
    return df

def load_summary(csv_path):
    """Load moment correction summary."""
    df = pd.read_csv(csv_path)
    return df

# ============================================================
# Find all result files
# ============================================================
sample_files = sorted(glob.glob(os.path.join(RESULT_DIR, "cn_normality_*_sprime_samples.csv")))
summary_files = sorted(glob.glob(os.path.join(RESULT_DIR, "cn_normality_*_moment_summary.csv")))

if not sample_files:
    print("No sample files found. Run run_cn_normality_tests.sh first.")
    exit(1)

print(f"Found {len(sample_files)} sample files, {len(summary_files)} summary files")

# ============================================================
# 1. Normality analysis: Q-Q plots + Shapiro-Wilk
# ============================================================
outdir_norm = os.path.join(RESULT_DIR, "normality")
os.makedirs(outdir_norm, exist_ok=True)

all_normality_results = []

for sf in sample_files:
    # Parse dataset and eps from filename
    basename = os.path.basename(sf).replace("cn_normality_", "").replace("_sprime_samples.csv", "")
    # e.g. "co_eps1.0" -> dataset=co, eps=1.0
    parts = basename.rsplit("_eps", 1)
    ds_name = parts[0]
    eps_str = parts[1] if len(parts) > 1 else "?"

    df = load_samples(sf)
    pairs = df.groupby(['u', 'v', 'true_S'])

    # Pick up to 6 pairs with largest true_S for Q-Q plots
    pair_keys = sorted(pairs.groups.keys(), key=lambda x: x[2], reverse=True)
    plot_pairs = pair_keys[:6]

    if len(plot_pairs) == 0:
        continue

    ncols = min(3, len(plot_pairs))
    nrows = (len(plot_pairs) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    if nrows * ncols == 1:
        axes = np.array([axes])
    axes = np.array(axes).flatten()

    for idx, (u, v, true_S) in enumerate(plot_pairs):
        ax = axes[idx]
        group = pairs.get_group((u, v, true_S))
        samples = group['S_prime'].values

        # Shapiro-Wilk (use at most 5000 samples)
        test_samples = samples[:5000]
        if len(test_samples) >= 20:
            stat, pval = sp_stats.shapiro(test_samples)
        else:
            stat, pval = np.nan, np.nan

        all_normality_results.append({
            'dataset': ds_name, 'eps': eps_str,
            'u': u, 'v': v, 'true_S': true_S,
            'n_samples': len(samples),
            'mean': np.mean(samples), 'std': np.std(samples),
            'shapiro_stat': stat, 'shapiro_pval': pval
        })

        # Q-Q plot
        sp_stats.probplot(samples, dist="norm", plot=ax)
        ax.set_title(f"S={true_S}, W={stat:.4f}, p={pval:.4f}", fontsize=9)
        ax.get_lines()[0].set_markersize(2)

    for j in range(len(plot_pairs), len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Q-Q plots: {ds_name}, ε={eps_str}", fontsize=12)
    plt.tight_layout()
    outpath = os.path.join(outdir_norm, f"qq_{ds_name}_eps{eps_str}.pdf")
    fig.savefig(outpath, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"  Saved {outpath}")

# Print normality summary
print("\n=== Shapiro-Wilk Normality Test Summary ===")
print(f"{'Dataset':<12} {'ε':>4} {'(u,v)':>8} {'S':>5} {'W':>8} {'p-value':>10} {'Normal?':>8}")
print("-" * 65)
for r in all_normality_results:
    normal = "YES" if r['shapiro_pval'] > 0.05 else "no"
    print(f"{r['dataset']:<12} {r['eps']:>4} ({r['u']},{r['v']}){r['true_S']:>5} "
          f"{r['shapiro_stat']:>8.4f} {r['shapiro_pval']:>10.4f} {normal:>8}")

# ============================================================
# 2. Unbiasedness verification for q=2,3,4
# ============================================================
print("\n=== Moment Correction Unbiasedness (q=2,3,4) ===")
print(f"{'Dataset':<10} {'ε':>4} {'S':>5} "
      f"{'E[C2]':>12} {'true_C2':>12} {'relErr':>8} "
      f"{'E[C3]':>12} {'true_C3':>12} {'relErr':>8} "
      f"{'E[C4]':>12} {'true_C4':>12} {'relErr':>8}")
print("-" * 130)

all_bias_results = []

for sf in summary_files:
    basename = os.path.basename(sf).replace("cn_normality_", "").replace("_moment_summary.csv", "")
    parts = basename.rsplit("_eps", 1)
    ds_name = parts[0]
    eps_str = parts[1] if len(parts) > 1 else "?"

    df = load_summary(sf)
    for _, row in df.iterrows():
        true_S = row['true_S']
        if true_S < 3:
            continue  # skip pairs with very small S (C3, C4 = 0)

        for q in [2, 3, 4]:
            est = row[f'mean_C{q}_hat']
            true_val = row[f'true_C{q}']
            rel_err = abs(est - true_val) / max(abs(true_val), 1e-10)
            all_bias_results.append({
                'dataset': ds_name, 'eps': eps_str, 'true_S': true_S,
                'q': q, 'estimate': est, 'true_val': true_val, 'rel_err': rel_err
            })

        re2 = abs(row['mean_C2_hat'] - row['true_C2']) / max(abs(row['true_C2']), 1e-10)
        re3 = abs(row['mean_C3_hat'] - row['true_C3']) / max(abs(row['true_C3']), 1e-10)
        re4 = abs(row['mean_C4_hat'] - row['true_C4']) / max(abs(row['true_C4']), 1e-10)
        print(f"{ds_name:<10} {eps_str:>4} {int(true_S):>5} "
              f"{row['mean_C2_hat']:>12.2f} {row['true_C2']:>12.2f} {re2:>8.4f} "
              f"{row['mean_C3_hat']:>12.2f} {row['true_C3']:>12.2f} {re3:>8.4f} "
              f"{row['mean_C4_hat']:>12.2f} {row['true_C4']:>12.2f} {re4:>8.4f}")

# ============================================================
# 3. Summary plot: relative error of moment correction vs true_S
# ============================================================
if all_bias_results:
    bias_df = pd.DataFrame(all_bias_results)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for idx, q in enumerate([2, 3, 4]):
        ax = axes[idx]
        sub = bias_df[bias_df['q'] == q]
        for (ds, eps), grp in sub.groupby(['dataset', 'eps']):
            ax.scatter(grp['true_S'], grp['rel_err'], label=f"{ds} ε={eps}", s=20, alpha=0.7)
        ax.set_xlabel("True S(u,v)")
        ax.set_ylabel("Relative Error")
        ax.set_title(f"C(S, {q}) estimator")
        ax.set_yscale('log')
        ax.axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, label='10% error')
        ax.legend(fontsize=6, ncol=2)
        ax.grid(True, alpha=0.3)

    plt.suptitle("Moment Correction: Relative Error vs True Common Neighbors", fontsize=13)
    plt.tight_layout()
    outpath = os.path.join(outdir_norm, "moment_correction_accuracy.pdf")
    fig.savefig(outpath, bbox_inches='tight', dpi=150)
    print(f"\nSaved {outpath}")
    plt.close(fig)

# ============================================================
# 4. Distribution plot: overlay S' histogram with normal fit
# ============================================================
for sf in sample_files:
    basename = os.path.basename(sf).replace("cn_normality_", "").replace("_sprime_samples.csv", "")
    parts = basename.rsplit("_eps", 1)
    ds_name = parts[0]
    eps_str = parts[1] if len(parts) > 1 else "?"

    df = load_samples(sf)
    pairs = df.groupby(['u', 'v', 'true_S'])
    pair_keys = sorted(pairs.groups.keys(), key=lambda x: x[2], reverse=True)
    plot_pairs = pair_keys[:4]

    if len(plot_pairs) == 0:
        continue

    fig, axes = plt.subplots(1, len(plot_pairs), figsize=(5*len(plot_pairs), 4))
    if len(plot_pairs) == 1:
        axes = [axes]

    for idx, (u, v, true_S) in enumerate(plot_pairs):
        ax = axes[idx]
        samples = pairs.get_group((u, v, true_S))['S_prime'].values
        mu, sigma = np.mean(samples), np.std(samples)

        ax.hist(samples, bins=50, density=True, alpha=0.6, color='steelblue', label="S' samples")
        x = np.linspace(mu - 4*sigma, mu + 4*sigma, 200)
        ax.plot(x, sp_stats.norm.pdf(x, mu, sigma), 'r-', lw=2, label=f'N({mu:.1f}, {sigma:.1f}²)')
        ax.axvline(true_S, color='green', linestyle='--', lw=2, label=f'True S={true_S}')
        ax.set_title(f"(u={u},v={v}), S={true_S}", fontsize=10)
        ax.legend(fontsize=7)

    fig.suptitle(f"S'(u,v) distribution: {ds_name}, ε={eps_str}", fontsize=12)
    plt.tight_layout()
    outpath = os.path.join(outdir_norm, f"dist_{ds_name}_eps{eps_str}.pdf")
    fig.savefig(outpath, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"  Saved {outpath}")

print("\nDone!")
