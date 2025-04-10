#pragma once
#ifndef __BIGRAPH_H
#define __BIGRAPH_H
#include "utility.h"

using namespace std;
using namespace std::chrono;

class BiGraph {
   public:
    BiGraph(std::string dir);
    BiGraph(const BiGraph& other);
    BiGraph();
    ~BiGraph() {
    }
    // note that addEdgeRaw is processing raw ids.
    void addEdgeRaw(vid_t u, vid_t v);

    // addEdge is processing real ids.
    void addEdge(vid_t u, vid_t v);
    // deleteEdge is processing real ids.
    void deleteEdge(vid_t u, vid_t v);

   public:
    void computePriority();
    void init(unsigned int num_v1, unsigned int num_v2);
    void loadGraph(std::string dir);
    void print_graph();
    void show_nodes();
    // vector<LinearHeap> make_heaps();
    bool same_layer(vid_t u, vid_t v) {
        if (is_upper(u) && is_upper(v)) return true;
        if (is_lower(u) && is_lower(v)) return true;
        return false;
    }

    bool diff_layer(vid_t u, vid_t v) { return !same_layer(u, v); }

    bool is_upper(vid_t u) { return (u < this->num_v1); }
    bool is_lower(vid_t u) { return (u >= this->num_v1); }
    // void generalGraph(string dir);

    bool com_p(vid_t u, vid_t v) {
        // return prio[u] > prio[v];
        return prio[u] < prio[v];
    }

    bool is_active(int vertex) { return degree[vertex] > 0; }
    int num_nodes() { return this->num_v1 + this->num_v2; }

    // checks the neighbors 
    // check whether i and j are connected 
    bool has(int i, int j) {
        return std::find(neighbor[i].begin(), neighbor[i].end(), j) !=
               neighbor[i].end();
    }

    // for each pair of vertices u, v, it checks whether this pair has been computed. 
	bool has_computed(int u, int v){
		assert(u!=v);
		int smaller = (u < v) ? u : v;
		int larger = (u < v) ? v : u;
		return edge_vector[smaller].find(larger) != edge_vector[smaller].end();
	}



    bool hasEdge(int u, int v) {
        if (u < 0 || u >= num_nodes() || v < 0 || v >= num_nodes()) {
            cout << "Invalid vertex index\n";
            // return false;
            exit(1);
        }
        // Check if u and v belong to different layers
        if (same_layer(u, v)) {
            // Invalid edge between upper and lower vertices
            cout << "Invalid edge between upper and lower vertices" << endl;
            // return false;
            exit(1);
        }
        // Get the degrees of vertices u and v
        int degreeU = degree[u];
        int degreeV = degree[v];

        // Determine which vertex to scan first based on their degrees
        int vertexToScanFirst = (degreeU < degreeV) ? u : v;
        int vertexToScanSecond = (vertexToScanFirst == u) ? v : u;

        // Check if vertexToScanFirst is connected to vertexToScanSecond
        for (const auto& neighbor : neighbor[vertexToScanFirst]) {
            if (neighbor == vertexToScanSecond) {
                // Edge found
                return true;
            }
        }
        // Edge not found
        return false;
    }

    std::string dir;
    num_t num_v1;
    num_t num_v2;
    long num_edges;
    vector<vector<vid_t>> neighbor;
    vector<int> degree;
    vector<int> prio;
    vector<int> U, V;
    vector<int> nodes;
   public:
    // max and min degrees
    int v1_max_degree;
    int v2_max_degree;

    vector<map<int, bool>> edge_vector;
};

#endif /* __BIGRAPH_H */