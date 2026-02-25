# AGENTS.md

## Cursor Cloud specific instructions

### Overview

This is a C++ research project (`ldp-pq`) for (p,q)-biclique counting in bipartite graphs under edge-local differential privacy. It produces command-line binaries — no web server, no containerization.

### Build

See `README.md` for full usage. Quick reference:

```bash
make clean && make          # builds ./biclique (main binary)
make variance_analysis      # builds ./variance_analysis
```

The `make testp3` target has a pre-existing linker error (missing globals from `main.cpp` not linked in that target). Use the main `biclique` binary for testing.

### System dependencies

- `g++` (C++14/17), `make`, `libsqlite3-dev`, `libgomp1` (OpenMP). These are installed via the update script.

### Running the binary

```bash
./biclique <epsilon> <data_dir> <num_iterations> <algorithm_switch> <p> <q>
```

- **Data files**: The binary expects `<data_dir>.e` (edge list) and `<data_dir>.meta` (vertex counts) files. These are not tracked in git (`.gitignore` excludes non-source files).
- **Ground truth DB**: The code opens `../biclq_counts.db` (relative to CWD) via SQLite. If the table/row is missing, it computes the exact count and inserts it. For the cloud VM, this file is at `/biclq_counts.db` when running from `/workspace`.
- **Test dataset**: A minimal synthetic dataset is at `testdata/small` (3 upper × 3 lower vertices, 6 edges, 3 butterflies). Example: `./biclique 1.0 ./testdata/small 2 0 2 2`

### Gotchas

- The `.gitignore` uses a whitelist pattern (`*` then `!*.cpp` etc.), so new non-source files (like test data, databases) won't be tracked by default.
- Algorithm switch values: 0=Naive, 1=OneR, 2=ADV, 3=ADV+, 4=ADV++, 6=PBC, 7=CN normality test, 8=CN general q.
- No linter or formatter is configured for this C++ project. Compilation with warnings is the primary code quality check.
