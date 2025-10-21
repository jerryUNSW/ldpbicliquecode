#!/usr/bin/env python3
"""
Download real bipartite datasets from Network Repository
"""

import requests
import os
import zipfile
import tarfile
import numpy as np
from collections import defaultdict

def download_from_networkrepo(dataset_name):
    """Download dataset from Network Repository"""
    
    # Network Repository base URLs
    base_urls = [
        f"https://nrvis.com/download/data/{dataset_name}.zip",
        f"https://nrvis.com/download/data/{dataset_name}.tar.gz", 
        f"https://nrvis.com/download/data/{dataset_name}.mtx",
        f"https://nrvis.com/download/data/{dataset_name}.edges",
        f"https://nrvis.com/download/data/{dataset_name}.txt"
    ]
    
    print(f"Attempting to download {dataset_name} from Network Repository...")
    
    for url in base_urls:
        try:
            print(f"Trying: {url}")
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                filename = url.split('/')[-1]
                with open(filename, 'wb') as f:
                    f.write(response.content)
                print(f"Successfully downloaded: {filename}")
                return filename
        except Exception as e:
            print(f"Failed to download {url}: {e}")
            continue
    
    return None

def extract_and_process_file(filename):
    """Extract and process the downloaded file"""
    
    print(f"Processing {filename}...")
    
    # Try to extract if it's an archive
    if filename.endswith('.zip'):
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall('.')
            # Find the extracted file
            for extracted_file in zip_ref.namelist():
                if extracted_file.endswith(('.mtx', '.edges', '.txt')):
                    return extracted_file
    elif filename.endswith('.tar.gz'):
        with tarfile.open(filename, 'r:gz') as tar_ref:
            tar_ref.extractall('.')
            # Find the extracted file
            for member in tar_ref.getmembers():
                if member.name.endswith(('.mtx', '.edges', '.txt')):
                    return member.name
    
    return filename

def convert_to_meta_e_format(input_file, dataset_name):
    """Convert the dataset to .meta and .e format"""
    
    print(f"Converting {input_file} to .meta and .e format...")
    
    edges = []
    vertices1 = set()
    vertices2 = set()
    
    # Read the file
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith('%') or line.startswith('#'):
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
            
            # Limit for very large files
            if line_num > 1000000:  # 1M edges max
                break
    
    if not edges:
        print("No edges found in the file")
        return None, None, None
    
    print(f"Read {len(edges)} edges")
    
    n1 = len(vertices1)
    n2 = len(vertices2)
    
    print(f"Partition 1: {n1} vertices")
    print(f"Partition 2: {n2} vertices")
    
    # Create mapping to consecutive integers starting from 0
    mapping1 = {v: i for i, v in enumerate(sorted(vertices1))}
    mapping2 = {v: i for i, v in enumerate(sorted(vertices2))}
    
    # Write .meta file
    meta_filename = f"{dataset_name}.meta"
    with open(meta_filename, 'w') as f:
        f.write(f"# Bipartite graph: {dataset_name} network\n")
        f.write(f"# n1: {n1}\n")
        f.write(f"# n2: {n2}\n")
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
    """Main function to download real datasets"""
    
    print("="*60)
    print("REAL NETWORK REPOSITORY DATASET DOWNLOAD")
    print("="*60)
    
    # List of datasets to try
    datasets = [
        "actor",
        "dbpedteam-bi", 
        "bibsonomy-2ut",
        "bio-diseasome"
    ]
    
    for dataset_name in datasets:
        print(f"\nTrying to download {dataset_name}...")
        
        # Download dataset
        filename = download_from_networkrepo(dataset_name)
        
        if filename:
            # Extract if needed
            processed_file = extract_and_process_file(filename)
            
            if processed_file:
                # Convert to .meta and .e format
                n1, n2, n_edges = convert_to_meta_e_format(processed_file, dataset_name)
                
                if n1 and n2:
                    print(f"\nSuccessfully processed {dataset_name}:")
                    print(f"  Partition 1: {n1} vertices")
                    print(f"  Partition 2: {n2} vertices")
                    print(f"  Edges: {n_edges}")
                    print(f"  Density: {n_edges/(n1*n2):.6f}")
                    
                    # Check if it meets our criteria (small partition)
                    if n1 < 1000 or n2 < 1000:
                        print(f"✅ {dataset_name} has a small partition - good for experiments!")
                        break
                    else:
                        print(f"⚠️  {dataset_name} has large partitions - may be too big")
                else:
                    print(f"❌ Failed to process {dataset_name}")
            else:
                print(f"❌ Failed to extract {filename}")
        else:
            print(f"❌ Failed to download {dataset_name}")
    
    print("\n" + "="*60)
    print("DOWNLOAD COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()

