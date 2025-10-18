#!/usr/bin/env python3
"""
Helper script to calculate and print fill rates for all datasets in bidata directory.
Fill rate = |E| / (|U| * |L|) where |E| is number of edges, |U| is upper layer size, |L| is lower layer size.
"""

import os
import glob
import sys

def read_meta_file(meta_path):
    """Read meta file and return (upper_size, lower_size, num_edges)"""
    try:
        with open(meta_path, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 3:
            return None, None, None
            
        upper_size = int(lines[0].strip())
        lower_size = int(lines[1].strip())
        num_edges = int(lines[2].strip())
        
        return upper_size, lower_size, num_edges
    except (ValueError, FileNotFoundError, IndexError) as e:
        print(f"Error reading {meta_path}: {e}")
        return None, None, None

def calculate_fill_rate(upper_size, lower_size, num_edges):
    """Calculate fill rate = |E| / (|U| * |L|)"""
    if upper_size is None or lower_size is None or num_edges is None:
        return None
    if upper_size <= 0 or lower_size <= 0:
        return None
    
    return num_edges / (upper_size * lower_size)

def main():
    bidata_dir = "/data/yizhangh/bidata"
    
    if not os.path.exists(bidata_dir):
        print(f"Error: bidata directory {bidata_dir} does not exist")
        sys.exit(1)
    
    # Find all .meta files
    meta_files = glob.glob(os.path.join(bidata_dir, "*.meta"))
    meta_files.sort()
    
    print("Dataset Fill Rate Analysis")
    print("=" * 80)
    print(f"{'Dataset':<30} {'|U|':<8} {'|L|':<8} {'|E|':<10} {'Fill Rate':<12} {'Density'}")
    print("-" * 80)
    
    fill_rates = []
    
    for meta_file in meta_files:
        dataset_name = os.path.basename(meta_file).replace('.meta', '')
        
        upper_size, lower_size, num_edges = read_meta_file(meta_file)
        
        if upper_size is None:
            print(f"{dataset_name:<30} {'ERROR':<8} {'ERROR':<8} {'ERROR':<10} {'ERROR':<12} {'ERROR'}")
            continue
            
        fill_rate = calculate_fill_rate(upper_size, lower_size, num_edges)
        
        if fill_rate is None:
            print(f"{dataset_name:<30} {upper_size:<8} {lower_size:<8} {num_edges:<10} {'ERROR':<12} {'ERROR'}")
            continue
            
        # Categorize density
        if fill_rate >= 0.1:
            density = "Very High"
        elif fill_rate >= 0.01:
            density = "High"
        elif fill_rate >= 0.001:
            density = "Medium"
        else:
            density = "Low"
            
        print(f"{dataset_name:<30} {upper_size:<8} {lower_size:<8} {num_edges:<10} {fill_rate:<12.6f} {density}")
        
        fill_rates.append((dataset_name, fill_rate, upper_size, lower_size, num_edges))
    
    # Sort by fill rate and show top datasets
    fill_rates.sort(key=lambda x: x[1], reverse=True)
    
    print("\n" + "=" * 80)
    print("TOP 20 DATASETS BY FILL RATE:")
    print("=" * 80)
    print(f"{'Rank':<4} {'Dataset':<30} {'Fill Rate':<12} {'|U|':<8} {'|L|':<8} {'|E|':<10}")
    print("-" * 80)
    
    for i, (dataset_name, fill_rate, upper_size, lower_size, num_edges) in enumerate(fill_rates[:20], 1):
        print(f"{i:<4} {dataset_name:<30} {fill_rate:<12.6f} {upper_size:<8} {lower_size:<8} {num_edges:<10}")
    
    # Summary statistics
    valid_fill_rates = [fr[1] for fr in fill_rates if fr[1] is not None]
    if valid_fill_rates:
        print(f"\nSUMMARY STATISTICS:")
        print(f"Total datasets: {len(valid_fill_rates)}")
        print(f"Highest fill rate: {max(valid_fill_rates):.6f}")
        print(f"Lowest fill rate: {min(valid_fill_rates):.6f}")
        print(f"Average fill rate: {sum(valid_fill_rates)/len(valid_fill_rates):.6f}")
        
        # Count by density category
        very_high = sum(1 for fr in valid_fill_rates if fr >= 0.1)
        high = sum(1 for fr in valid_fill_rates if 0.01 <= fr < 0.1)
        medium = sum(1 for fr in valid_fill_rates if 0.001 <= fr < 0.01)
        low = sum(1 for fr in valid_fill_rates if fr < 0.001)
        
        print(f"\nDENSITY DISTRIBUTION:")
        print(f"Very High (≥0.1): {very_high} datasets")
        print(f"High (0.01-0.1): {high} datasets")
        print(f"Medium (0.001-0.01): {medium} datasets")
        print(f"Low (<0.001): {low} datasets")

if __name__ == "__main__":
    main()
