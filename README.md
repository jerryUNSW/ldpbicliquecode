# ldp-pq (Local Differential Privacy Biclique Counting)

## Overview

ldp-pq is a C++ project focused on biclique counting in bipartite graphs with edge-local differential privacy.

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
- **2**: Advanced base algorithm (adv-base).
- **3**: Advanced+ algorithm (adv+).
- **4**: Advanced++ algorithm (adv++).

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

## Supported Datasets

The program supports the following datasets (located in `../bidata/`):

- `to`
- `co`
- `big`
- `crime`
- `M_PL_030`
- `unicode`
- `lrcwiki`
- `librec-filmtrust-ratings`
- `rmwiki`
- `amazon-ratings`
- `edit-iowiktionary`
- `opsahl-collaboration`
- `csbwiki`
- `bag-kos`
- `bpywiki`
- `nips`
- `lastfm_band`
- `discogs_lstyle_lstyle`
- `digg-votes`
- `movielens-10m_rating`

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

## Notes

- The naive algorithm (switch 0) is optimized for a single round, as additional rounds do not improve its performance.
- The one-round algorithm (switch 1) is only feasible for smaller datasets like `to` and `co`.
- For large epsilon values, the naive algorithm may produce overly optimistic results, so use with caution.
- Output can be filtered for specific algorithm results using `grep adv` (e.g., for advanced algorithms).