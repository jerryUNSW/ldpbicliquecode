# ldp-pq (Local Differential Privacy Biclique Counting)

## Overview

ldp-pq is a C++ project focused on biclique counting in bipartite graphs with edge local differential privacy.


## Project Structure

The project consists of the following key files and directories:

- `main.cpp`: Entry point of the program.
- `biclique.cpp` / `biclique.h`: Implementation and declarations of biclique counting algorithms unde edge LDP. 
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


## Algorithm Switch Options

- **0**: Naive algorithm
- **1**: One-round algorithm
- **3**: Advanced algorithm


## Data Format for Bipartite Graphs

This program processes bipartite graph data, which consists of two files: an edge list file (`<datafile>.e`) and a metadata file (`<datafile>.meta`).

### Edge List File (`<datafile>.e`)

The edge list file represents the connections between upper and lower vertices in the bipartite graph. Each line in the file describes an edge between an upper vertex and a lower vertex. The format is as follows:

<upper_vertex> <lower_vertex>

### Metadata File (<datafile>.meta)
The metadata file provides essential information about the bipartite graph:

Upper Vertices Count: The number of upper vertices 
Lower Vertices Count: The number of lower vertices 
Edges Count: The total number of edges 

<upper_vertices_count>
<lower_vertices_count>
<edges_count>

## Example usage

1. Run the naive algorithm to count the numebr of (2,3)-bicliques for 10 rounds with a privacy budget epsilon = 2, on the dataset unicode:  
```bash
./biclique 2 unicode 10 2 3 
```

2. Run the ADV algorithm to count the numebr of (3,2)-bicliques for 10 rounds with a privacy budget epsilon = 1, on the dataset unicode:  
```bash
./biclique 1 unicode 10 3 2 
```