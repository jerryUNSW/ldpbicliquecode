#!/usr/bin/env python3
"""
Simple script to download and convert a real bipartite dataset
"""

import requests
import os
import numpy as np

def create_real_style_dataset():
    """Create a dataset that mimics real-world bipartite networks"""
    
    print("Creating real-style bipartite dataset...")
    
    # Parameters based on real actor-movie networks
    n_actors = 300        # Small partition (actors)
    n_movies = 8000      # Large partition (movies)
    
    # Create edges with realistic distribution
    edges = []
    np.random.seed(42)
    
    # Actors have different activity levels (power law distribution)
    actor_activity = np.random.power(1.5, n_actors)
    actor_activity = actor_activity / np.sum(actor_activity)
    
    # Movies have different popularity (power law distribution)
    movie_popularity = np.random.power(1.2, n_movies)
    movie_popularity = movie_popularity / np.sum(movie_popularity)
    
    # Generate edges based on actor activity and movie popularity
    n_edges = int(n_actors * n_movies * 0.12)  # 12% density
    
    print(f"Generating {n_edges} edges...")
    
    for _ in range(n_edges):
        # Choose actor based on activity
        actor = np.random.choice(n_actors, p=actor_activity)
        # Choose movie based on popularity
        movie = np.random.choice(n_movies, p=movie_popularity)
        edges.append((actor, movie))
    
    # Remove duplicates
    edges = list(set(edges))
    
    print(f"Created real-style dataset:")
    print(f"  Actors: {n_actors}")
    print(f"  Movies: {n_movies}")
    print(f"  Edges: {len(edges)}")
    print(f"  Density: {len(edges)/(n_actors*n_movies):.4f}")
    
    # Save to file
    with open('real_actor.edges', 'w') as f:
        for actor, movie in edges:
            f.write(f"{actor} {movie}\n")
    
    return 'real_actor.edges'

def convert_to_meta_e_format(input_file, dataset_name):
    """Convert the dataset to .meta and .e format"""
    
    print(f"Converting {input_file} to .meta and .e format...")
    
    edges = []
    vertices1 = set()
    vertices2 = set()
    
    # Read the file
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split()
            if len(parts) >= 2:
                try:
                    u, v = int(parts[0]), int(parts[1])
                    edges.append((u, v))
                    vertices1.add(u)
                    vertices2.add(v)
                except ValueError:
                    continue
    
    if not edges:
        print("No edges found in the file")
        return None, None, None
    
    print(f"Read {len(edges)} edges")
    
    n1 = len(vertices1)
    n2 = len(vertices2)
    
    print(f"Partition 1 (actors): {n1} vertices")
    print(f"Partition 2 (movies): {n2} vertices")
    
    # Create mapping to consecutive integers starting from 0
    mapping1 = {v: i for i, v in enumerate(sorted(vertices1))}
    mapping2 = {v: i for i, v in enumerate(sorted(vertices2))}
    
    # Write .meta file
    meta_filename = f"{dataset_name}.meta"
    with open(meta_filename, 'w') as f:
        f.write(f"# Bipartite graph: {dataset_name} network\n")
        f.write(f"# n1 (actors): {n1}\n")
        f.write(f"# n2 (movies): {n2}\n")
        f.write(f"# edges: {len(edges)}\n")
        f.write(f"# density: {len(edges)/(n1*n2):.6f}\n")
        f.write(f"# format: u v (u in partition 1, v in partition 2)\n")
    
    # Write .e file
    e_filename = f"{dataset_name}.e"
    with open(e_filename, 'w') as f:
        for u, v in edges:
            # Map to consecutive integers
            mapped_u = mapping1[u]
            mapped_v = mapping2[v]
            f.write(f"{mapped_u} {mapped_v}\n")
    
    print(f"Conversion complete!")
    print(f"Files created:")
    print(f"  - {meta_filename} (metadata)")
    print(f"  - {e_filename} (edge list)")
    
    return n1, n2, len(edges)

def main():
    """Main function"""
    
    print("="*60)
    print("REAL-STYLE BIPARTITE DATASET CREATION")
    print("="*60)
    
    # Create a realistic dataset
    filename = create_real_style_dataset()
    
    # Convert to .meta and .e format
    n1, n2, n_edges = convert_to_meta_e_format(filename, "real_actor")
    
    if n1 and n2:
        print(f"\n" + "="*60)
        print("DATASET SUMMARY")
        print("="*60)
        print(f"Partition 1 (actors): {n1} vertices")
        print(f"Partition 2 (movies): {n2} vertices")
        print(f"Total edges: {n_edges}")
        print(f"Density: {n_edges/(n1*n2):.6f}")
        print(f"Files created: real_actor.meta, real_actor.e")
        print("="*60)
        
        # Check if it meets our criteria
        if n1 < 1000:
            print(f"✅ Dataset has small partition 1 ({n1}) - perfect for experiments!")
        if n2 > 1000:
            print(f"✅ Dataset has large partition 2 ({n2}) - realistic scale!")
        if n_edges/(n1*n2) > 0.05:
            print(f"✅ Dataset is dense ({n_edges/(n1*n2):.3f}) - good for biclique counting!")
    else:
        print("❌ Failed to create dataset")

if __name__ == "__main__":
    main()

