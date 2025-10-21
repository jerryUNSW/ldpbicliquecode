#!/usr/bin/env python3
"""
Download datasets from Network Repository using curl and web scraping
"""

import subprocess
import os
import re
import requests
from bs4 import BeautifulSoup

def download_with_curl(dataset_name):
    """Download dataset using curl"""
    
    # Try different curl approaches
    curl_commands = [
        f"curl -L -o {dataset_name}.mtx 'https://nrvis.com/download/data/{dataset_name}.mtx'",
        f"curl -L -o {dataset_name}.edges 'https://nrvis.com/download/data/{dataset_name}.edges'",
        f"curl -L -o {dataset_name}.txt 'https://nrvis.com/download/data/{dataset_name}.txt'",
        f"curl -L -o {dataset_name}.zip 'https://nrvis.com/download/data/{dataset_name}.zip'"
    ]
    
    for cmd in curl_commands:
        try:
            print(f"Running: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
            if result.returncode == 0 and os.path.exists(cmd.split()[-1]):
                filename = cmd.split()[-1]
                file_size = os.path.getsize(filename)
                if file_size > 100:  # File has content
                    print(f"Successfully downloaded: {filename} ({file_size} bytes)")
                    return filename
        except Exception as e:
            print(f"Curl command failed: {e}")
            continue
    
    return None

def scrape_networkrepo_page(dataset_name):
    """Scrape the Network Repository page to find download links"""
    
    try:
        # Try to access the dataset page
        url = f"https://networkrepository.com/{dataset_name}.php"
        print(f"Scraping page: {url}")
        
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Look for download links
            download_links = []
            for link in soup.find_all('a', href=True):
                href = link['href']
                if any(ext in href for ext in ['.mtx', '.edges', '.txt', '.zip']):
                    download_links.append(href)
            
            print(f"Found {len(download_links)} potential download links")
            
            # Try to download from found links
            for link in download_links:
                if link.startswith('http'):
                    download_url = link
                else:
                    download_url = f"https://networkrepository.com/{link}"
                
                try:
                    print(f"Trying to download: {download_url}")
                    response = requests.get(download_url, timeout=30)
                    if response.status_code == 200:
                        filename = f"{dataset_name}_{link.split('/')[-1]}"
                        with open(filename, 'wb') as f:
                            f.write(response.content)
                        print(f"Successfully downloaded: {filename}")
                        return filename
                except Exception as e:
                    print(f"Failed to download {download_url}: {e}")
                    continue
        
    except Exception as e:
        print(f"Failed to scrape page: {e}")
    
    return None

def create_realistic_actor_dataset():
    """Create a more realistic actor dataset based on real-world statistics"""
    
    print("Creating realistic actor-movie dataset...")
    
    # Realistic parameters based on real actor-movie networks
    n_actors = 500        # More actors
    n_movies = 5000      # Fewer movies (more realistic)
    
    # Create edges with realistic distribution
    edges = []
    np.random.seed(42)
    
    # Actors have different activity levels (power law)
    actor_activity = np.random.power(1.5, n_actors)
    actor_activity = actor_activity / np.sum(actor_activity)
    
    # Movies have different popularity (power law)
    movie_popularity = np.random.power(1.2, n_movies)
    movie_popularity = movie_popularity / np.sum(movie_popularity)
    
    # Generate edges based on actor activity and movie popularity
    n_edges = int(n_actors * n_movies * 0.08)  # 8% density
    
    for _ in range(n_edges):
        # Choose actor based on activity
        actor = np.random.choice(n_actors, p=actor_activity)
        # Choose movie based on popularity
        movie = np.random.choice(n_movies, p=movie_popularity)
        edges.append((actor, movie))
    
    # Remove duplicates
    edges = list(set(edges))
    
    print(f"Created realistic dataset:")
    print(f"  Actors: {n_actors}")
    print(f"  Movies: {n_movies}")
    print(f"  Edges: {len(edges)}")
    print(f"  Density: {len(edges)/(n_actors*n_movies):.4f}")
    
    # Save to file
    with open('realistic_actor.edges', 'w') as f:
        for actor, movie in edges:
            f.write(f"{actor} {movie}\n")
    
    return 'realistic_actor.edges'

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
    """Main function"""
    
    print("="*60)
    print("REAL NETWORK REPOSITORY DOWNLOAD (WITH CURL)")
    print("="*60)
    
    # Try to download real datasets
    datasets = ["actor", "dbpedteam-bi", "bibsonomy-2ut"]
    
    for dataset_name in datasets:
        print(f"\nTrying to download {dataset_name}...")
        
        # Try curl download
        filename = download_with_curl(dataset_name)
        
        if not filename:
            # Try web scraping
            filename = scrape_networkrepo_page(dataset_name)
        
        if filename:
            # Convert to .meta and .e format
            n1, n2, n_edges = convert_to_meta_e_format(filename, dataset_name)
            
            if n1 and n2:
                print(f"\nSuccessfully processed {dataset_name}:")
                print(f"  Partition 1: {n1} vertices")
                print(f"  Partition 2: {n2} vertices")
                print(f"  Edges: {n_edges}")
                print(f"  Density: {n_edges/(n1*n2):.6f}")
                
                if n1 < 1000 or n2 < 1000:
                    print(f"✅ {dataset_name} has a small partition - good for experiments!")
                    return
                else:
                    print(f"⚠️  {dataset_name} has large partitions - may be too big")
            else:
                print(f"❌ Failed to process {dataset_name}")
        else:
            print(f"❌ Failed to download {dataset_name}")
    
    # If all downloads failed, create a realistic synthetic dataset
    print("\nAll downloads failed. Creating realistic synthetic dataset...")
    import numpy as np
    filename = create_realistic_actor_dataset()
    n1, n2, n_edges = convert_to_meta_e_format(filename, "realistic_actor")
    
    print("\n" + "="*60)
    print("DOWNLOAD COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()

