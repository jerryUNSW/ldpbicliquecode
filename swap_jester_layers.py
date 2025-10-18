#!/usr/bin/env python3
"""
Script to create swapped versions of jester1 and jester2 datasets where upper and lower layers are reversed.
"""

import os

def swap_dataset_layers(original_meta, original_edges, new_meta, new_edges):
    """Swap upper and lower layers of a bipartite dataset"""
    
    # Read original meta file
    with open(original_meta, 'r') as f:
        lines = f.readlines()
    
    original_upper = int(lines[0].strip())
    original_lower = int(lines[1].strip())
    num_edges = int(lines[2].strip())
    
    print(f"Original: {original_upper} upper, {original_lower} lower, {num_edges} edges")
    print(f"Swapped:  {original_lower} upper, {original_upper} lower, {num_edges} edges")
    
    # Create new meta file
    with open(new_meta, 'w') as f:
        f.write(f"{original_lower}\n")  # New upper size
        f.write(f"{original_upper}\n")  # New lower size  
        f.write(f"{num_edges}\n")       # Same number of edges
        f.write("\n")                   # Empty line
    
    # Read and swap edge file
    with open(original_edges, 'r') as f_in, open(new_edges, 'w') as f_out:
        for line in f_in:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    upper_node = int(parts[0])
                    lower_node = int(parts[1])
                    # Swap: write lower_node as new upper, upper_node as new lower
                    f_out.write(f"{lower_node} {upper_node}\n")
    
    print(f"Created swapped dataset: {new_meta} and {new_edges}")

def main():
    bidata_dir = "/data/yizhangh/bidata"
    
    # Process jester1
    print("Processing jester1...")
    swap_dataset_layers(
        os.path.join(bidata_dir, "jester1.meta"),
        os.path.join(bidata_dir, "jester1.e"),
        os.path.join(bidata_dir, "jester1-swapped.meta"),
        os.path.join(bidata_dir, "jester1-swapped.e")
    )
    
    print("\nProcessing jester2...")
    # Process jester2
    swap_dataset_layers(
        os.path.join(bidata_dir, "jester2.meta"),
        os.path.join(bidata_dir, "jester2.e"),
        os.path.join(bidata_dir, "jester2-swapped.meta"),
        os.path.join(bidata_dir, "jester2-swapped.e")
    )
    
    print("\nVerification:")
    # Verify the swapped datasets
    for dataset in ["jester1-swapped", "jester2-swapped"]:
        meta_file = os.path.join(bidata_dir, f"{dataset}.meta")
        if os.path.exists(meta_file):
            with open(meta_file, 'r') as f:
                lines = f.readlines()
            upper = int(lines[0].strip())
            lower = int(lines[1].strip())
            edges = int(lines[2].strip())
            fill_rate = edges / (upper * lower)
            print(f"{dataset}: {upper}×{lower}, {edges} edges, fill_rate={fill_rate:.6f}")

if __name__ == "__main__":
    main()
