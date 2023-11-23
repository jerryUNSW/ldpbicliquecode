#pragma once
#ifndef __BIGRAPH_H
#define __BIGRAPH_H
#include "utility.h"

using namespace std;
using namespace std::chrono; 

// std::mt19937 rng(std::random_device{}()); // for seeding

// for each root, store the other vertices in the CC of root. 

class UnionFind 
{

public:
	int *parent, *size;
	int n;

	int *next;

	UnionFind(int n) {
		this->n = n;
		parent = new int[n + 1];
		size = new int[n + 1];

		// this is new
		next = new int[n + 1];

		for (int i = 0; i <= n; i++) {
			parent[i] = i;
			size[i] = 1;
			//
			next[i]=i;
		}
	}

	int get_size_of_root(int xxx){
	//   if(xxx==0){
	// 	  printf("cannot visit the 0th element of size\n");
	// 	  exit(1);
	//   }
		return size[xxx];
	}

	int find(int x) {
		if (parent[x] == x) // path compression
			return x;

		return parent[x] = find(parent[x]);
	}

	void unite(int x, int y) {

		// assert(x != 0 && y!=0);
		int root_x = find(x);
		int root_y = find(y);
		if (root_x == root_y) 
			return;

		// at this point we know x and y belong to diff groups
		if (size[root_x] < size[root_y]) {
			// union by size, always manipulate the roots
			parent[root_x] = root_y;
			size[root_y] += size[root_x];
		} else {
			parent[root_y] = root_x;
			size[root_x] += size[root_y];
		}

		// swap next[x] and next[y]
		int tmp = next[x]; 
		next[x] = next[y];
		next[y] = tmp; 
	}
};

class BiGraph{
public:
	BiGraph(std::string dir);
	BiGraph(const BiGraph& other); 
	BiGraph();
	~BiGraph() {
		// destroy();
		// cout<<"destorying the current bigraph "<<endl;
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
	bool same_layer(vid_t u, vid_t v ){
		if(is_upper(u) && is_upper(v)) return true; 
		if(is_lower(u) && is_lower(v)) return true; 
		return false ; 
	}

	bool diff_layer(vid_t u, vid_t v ){
		return !same_layer(u,v);
	}

	bool is_upper(vid_t u){ 
		return (u < this->num_v1);
	}
	bool is_lower(vid_t u){
		return (u>= this->num_v1);
	}
	// void generalGraph(string dir);

	bool com_p(vid_t u, vid_t v){
		// return prio[u] > prio[v]; 
		return prio[u] < prio[v]; 
	}

	bool is_active(int vertex){ 
		return degree[vertex]>0;
	}
	int num_nodes(){ 
		return this->num_v1+this->num_v2;
	}

	bool has(int i, int j){
		return std::find(neighbor[i].begin(), neighbor[i].end(), j) != neighbor[i].end(); 
	}
	bool hasEdge(int u, int v) {
		if (u < 0 || u >= num_nodes() || v < 0 || v >= num_nodes()) {
			cout<<"Invalid vertex index\n";
			// return false;
			exit(1);
		}
		// Check if u and v belong to different layers
		if (same_layer(u,v)) {
			// Invalid edge between upper and lower vertices
			cout<<"Invalid edge between upper and lower vertices"<<endl;
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

	// void show_amat(); 
	
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
};


#endif  /* __BIGRAPH_H */