#!/usr/bin/env python3
"""
Script to swap upper and lower layers for Discogs dataset.
Creates a new dataset where the original lower vertices become upper vertices
and the original upper vertices become lower vertices.
"""

import os

def swap_discogs_layers():
    # Input files
    input_edge_file = "../bidata/discogs.e"
    input_meta_file = "../bidata/discogs.meta"
    
    # Output files
    output_edge_file = "../bidata/discogs-dl.e"
    output_meta_file = "../bidata/discogs-dl.meta"
    
    print("Reading original discogs metadata...")
    with open(input_meta_file, 'r') as f:
        lines = f.read().strip().split('\n')
        original_upper = int(lines[0])
        original_lower = int(lines[1])
        original_edges = int(lines[2])
    
    print(f"Original dimensions:")
    print(f"  Upper vertices: {original_upper:,}")
    print(f"  Lower vertices: {original_lower:,}")
    print(f"  Edges: {original_edges:,}")
    
    print(f"\nSwapped dimensions:")
    print(f"  Upper vertices: {original_lower:,} (was lower)")
    print(f"  Lower vertices: {original_upper:,} (was upper)")
    print(f"  Edges: {original_edges:,} (unchanged)")
    
    print("\nProcessing edge file...")
    edge_count = 0
    
    with open(input_edge_file, 'r') as infile, open(output_edge_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split()
            if len(parts) >= 2:
                original_upper_vertex = int(parts[0])
                original_lower_vertex = int(parts[1])
                
                # Swap the vertices
                new_upper_vertex = original_lower_vertex
                new_lower_vertex = original_upper_vertex
                
                outfile.write(f"{new_upper_vertex} {new_lower_vertex}\n")
                edge_count += 1
                
                if edge_count % 1000000 == 0:
                    print(f"  Processed {edge_count:,} edges...")
    
    print(f"Processed {edge_count:,} edges total")
    
    # Write new metadata file
    print("Writing new metadata file...")
    with open(output_meta_file, 'w') as f:
        f.write(f"{original_lower}\n")  # New upper count (was lower)
        f.write(f"{original_upper}\n")  # New lower count (was upper)
        f.write(f"{original_edges}\n")  # Edge count unchanged
    
    print(f"\nCreated new dataset:")
    print(f"  Edge file: {output_edge_file}")
    print(f"  Meta file: {output_meta_file}")
    print(f"  New dimensions: {original_lower} upper, {original_upper} lower, {original_edges} edges")

if __name__ == "__main__":
    swap_discogs_layers()
