# ldp-pq (Local Differential Privacy (p,q)-Biclique Counting)

## Overview

ldp-pq is a C++ project focused on (p,q)-biclique counting in bipartite graphs with edge-local differential privacy.

## Project Structure

The project consists of the following key files and directories:

- `main.cpp`: Entry point of the program.
- `biclique.cpp` / `biclique.h`: Implementation and declarations of biclique counting algorithms under edge LDP.
- `bigraph.cpp` / `bigraph.h`: Functionality and data structures for handling bipartite graphs.
- `utility.cpp` / `utility.h`: Common utility functions shared across modules.
- `exactcounting/`: Directory containing exact biclique counting experiment code.
- `include/`: Additional header files.
- `makefile`: Build script to compile the project using `make`.
- `README.md`: Documentation with project overview and usage instructions.

## Build Instructions

To build the project, use the following command:

```bash
make clean && make
```

## Running the Program

To run the ldp-pq program, use the following command:

```bash
./biclique <epsilon> <data_directory> <num_iterations> <algorithm_switch> <p> <q>
```

### Parameters

- `<epsilon>`: Privacy budget for edge-local differential privacy (e.g., 1). Note: When epsilon is too large, the naive algorithm may perform unrealistically well.
- `<data_directory>`: Path to the dataset directory (e.g., `../bidata/<dataset>`).
- `<num_iterations>`: Number of rounds to run the algorithm (e.g., 10). For the naive algorithm, one round is sufficient.
- `<algorithm_switch>`: Specifies the algorithm to use (see options below).
- `<p> <q>`: Parameters defining the size of the biclique to count (e.g., (p,q)-bicliques).

## Algorithm Switch Options

- **0**: Naive algorithm (single round recommended).
- **1**: One-round algorithm (feasible only on smaller datasets like "to" and "co").
- **2**: The MRCN algorithm.
- **3**: The MRCN algorithm + Multi-center optimization.
- **4**: The MRCN algorithm + Multi-center optimization + Refined Noisy Graph Construction.

## Data Format for Bipartite Graphs

The program processes bipartite graph data, which consists of two files: an edge list file (`<datafile>.e`) and a metadata file (`<datafile>.meta`).

### Edge List File (`<datafile>.e`)

The edge list file represents the connections between upper and lower vertices in the bipartite graph. Each line describes an edge in the format:

```
<upper_vertex> <lower_vertex>
```

### Metadata File (`<datafile>.meta`)

The metadata file provides essential information about the bipartite graph in the following format:

```
<upper_vertices_count>
<lower_vertices_count>
<edges_count>
```

- `Upper Vertices Count`: The number of upper vertices.
- `Lower Vertices Count`: The number of lower vertices.
- `Edges Count`: The total number of edges.

## Example Usage

1. Run the naive algorithm to count (2,3)-bicliques for 1 round with a privacy budget epsilon = 1 on the dataset `to`:

```bash
./biclique 1 ../bidata/to 1 0 2 3
```

2. Run the one-round algorithm to count (3,2)-bicliques for 10 rounds with a privacy budget epsilon = 1 on the dataset `to`:

```bash
./biclique 1 ../bidata/to 10 1 3 2
```

3. Run the advanced++ algorithm to count (3,2)-bicliques for 10 rounds with a privacy budget epsilon = 1 on the dataset `unicode`:

```bash
./biclique 1 ../bidata/unicode 10 4 3 2
```

## Important: Ground Truth Overflow Issue

### Problem Description

When computing ground truth for large Q values (typically Q ≥ 6), biclique counts can exceed 64-bit integer limits. The SQLite database correctly stores these large numbers as `REAL` (floating point) values in scientific notation, but the original C++ code was using `sqlite3_column_int64()` to retrieve them, causing integer overflow.

### Symptoms

- Ground truth values for Q ≥ 6 appear as the maximum 64-bit integer: `9223372036854775807`
- Relative error calculations become meaningless (artificially inflated)
- Algorithm performance appears much worse than it actually is
- All algorithms show similar poor performance for large Q values

### Example

For LRCWiki dataset:
- **Q=6**: True count = 4.11×10^19, but retrieved as 9.22×10^18 (overflowed)
- **Q=7**: True count = 3.27×10^22, but retrieved as 9.22×10^18 (overflowed)
- **Q=8**: True count = 2.27×10^25, but retrieved as 9.22×10^18 (overflowed)
- **Q=9**: True count = 1.40×10^28, but retrieved as 9.22×10^18 (overflowed)

### Fix Applied

The code has been updated in `biclique.cpp` to properly handle both `INTEGER` and `REAL` database values:

```cpp
// Get the column type to handle both INTEGER and REAL values
int column_type = sqlite3_column_type(stmt, 1);
long double count_value;

if (column_type == SQLITE_INTEGER) {
    count_value = static_cast<long double>(sqlite3_column_int64(stmt, 1));
} else if (column_type == SQLITE_FLOAT) {
    count_value = static_cast<long double>(sqlite3_column_double(stmt, 1));
}
```

### Impact

- **Before fix**: Relative errors were artificially inflated due to wrong ground truth
- **After fix**: Relative errors reflect actual algorithm performance
- **Result**: Algorithms show much better performance for large Q values than previously reported

### Verification

To verify the fix is working, look for these indicators in the output:
- Ground truth values are displayed in scientific notation for large Q values
- Warning messages about precision loss (if any) are shown
- Relative errors are much smaller and more reasonable for Q ≥ 6

### Database Schema

The ground truth is stored in `../biclq_counts.db` with the following schema:
```sql
CREATE TABLE pqbiclique_counts (
    dataset TEXT NOT NULL,
    p INTEGER NOT NULL,
    q INTEGER NOT NULL,
    count UNSIGNED BIG INT NOT NULL,  -- Can store as REAL for large values
    PRIMARY KEY (dataset, p, q)
);
```

Large values are automatically stored as `REAL` (scientific notation) by SQLite when they exceed integer limits.

## Variance Reduction Experiments: Why Filtering and Truncation Don't Work

We investigated several post-processing strategies to reduce variance in the ADV++ (multi-round) biclique estimator. All were tested on `unicode` (255×615, 1255 edges) and `co` (74×257, 328 edges) with ε=2.0, P=2, K=2.

### Approaches Tested

| Mode | Description | Result |
|------|-------------|--------|
| Baseline | ADV++ with no modification | rel err ~1.50 |
| Noisy degree pre-filter | Skip pairs where `min(deg_est(u), deg_est(w)) < K` | rel err ~4.00 (worse) |
| Winsorize 95th pct | Clip `local_res` at 5th/95th percentile | rel err ~3.79 (worse) |
| Deg filter + Winsorize | Combined | rel err ~5.91 (worse) |
| Truncation (auto) | Clip `|local_res|` at 3× median | rel err ~8.24 (worse) |
| Deg filter + Truncation | Combined | rel err ~6.34 (worse) |

### Why They All Fail

The core estimator `compute_local_res(K, f_u_w, σ²)` is designed to be **unbiased**. For K=2:

```
local_res = f²/2 - f/2 - σ²/2
```

The `-σ²/2` term corrects for noise. When summed over all pairs, positive and negative `local_res` values cancel out to produce an unbiased total. This cancellation is essential.

**Filtering** removes pairs from the sum. Pairs with `cn(u,w) < K` contribute `E[local_res] ≈ 0`, so removing them doesn't help the expectation — but the noisy degree filter (`eps0 = 0.05ε`) is too inaccurate to distinguish signal from noise pairs. At ε=2.0, `noise_std = 1/0.1 = 10`, making degree estimates useless for vertices with degree < 20. The filter randomly drops signal pairs, introducing negative bias.

**Truncation/Winsorization** clips extreme `local_res` values. But the distribution is asymmetric: the positive tail carries real signal (pairs with many common neighbors produce large positive `local_res`), while negative values are needed for bias correction. Symmetric clipping destroys more signal than noise, causing massive negative bias.

**Post-filtering on `f_u_w`** (e.g., skip if `f_u_w < K`) introduces selection bias because the filtering decision depends on the same noisy quantity used in the estimator.

### Upper Bound (Public Degree Filter)

Using true (non-private) degrees as a filter gives a 25× RMSE improvement on `co` (RMSE: 7359 → 292). This shows the potential is huge, but we cannot access true degrees in the LDP setting, and noisy degrees at practical ε values are too noisy to replicate this.

### Conclusion

Post-processing approaches that modify or filter the output of `compute_local_res` break the estimator's unbiasedness. Variance reduction must happen **upstream** — either by reducing the noise in `f_u_w` itself, or by designing a fundamentally different estimator with lower variance.

### Code Reference

The post-processing modes are implemented in `wedge_based_two_round_2_K_biclique()` in `biclique.cpp` (controlled by `post_processing_mode` global, set via 13th command-line argument). They remain in the code for reference but should not be used in experiments.