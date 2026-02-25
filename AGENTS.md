## Cursor Cloud specific instructions

### Overview

ldp-pq is a C++ research project for (p,q)-biclique counting under edge-local differential privacy. It is a single-binary CLI tool (not a web app or multi-service system). See `README.md` for algorithm descriptions, usage examples, and experiment results.

### System dependency

`libsqlite3-dev` must be installed (`sudo apt-get install -y libsqlite3-dev`). The runtime library `libsqlite3-0` is typically present but the dev headers are not.

### Build

```
make clean && make
```

Additional targets: `make variance_analysis` (builds successfully). The `make testp3` target has a pre-existing linker error (missing globals defined only in `main.cpp`).

### Running

The binary expects two external resources:
- **Dataset files** at a path like `../bidata/<name>.meta` and `../bidata/<name>.e` (relative to the working directory).
- **Ground truth SQLite database** at `../biclq_counts.db` (relative to working directory). Table: `pqbiclique_counts(dataset, p, q, count)`. If the dataset/p/q combo is not found in the DB, the program computes the exact count ad-hoc.

Example (from `/workspace`):
```
./biclique 1.0 ../bidata/testgraph 1 0 2 2
```

### No lint/test framework

This project has no linter configuration or automated test suite. Verification is done by building and running the binary with datasets.
