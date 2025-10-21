#!/usr/bin/env python3
"""
Print NIPS dataset results for Naive, ADV, and ADV+ algorithms.
"""

import sqlite3

def print_nips_results():
    """Print NIPS dataset results"""
    
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    print("="*100)
    print("📊 NIPS DATASET RESULTS (p=3)")
    print("="*100)
    
    # Get NIPS results for each algorithm
    algorithms = ['Naive', 'ADV', 'ADV+']
    
    for algorithm in algorithms:
        print(f"\n🔹 {algorithm} ALGORITHM:")
        print("-" * 80)
        
        cursor.execute('''
            SELECT q, ground_truth, estimate, relative_error 
            FROM algorithm_comparison_p3 
            WHERE algorithm = ? AND p = 3 AND dataset = 'nips'
            ORDER BY q
        ''', (algorithm,))
        
        results = cursor.fetchall()
        
        if results:
            print(f"{'Q':<3} {'Ground Truth':<20} {'Estimate':<20} {'Rel Error (%)':<15}")
            print("-" * 80)
            
            for q, gt, estimate, rel_error in results:
                print(f"{q:<3} {gt:<20.2e} {estimate:<20.2e} {rel_error:<15.2f}")
        else:
            print(f"❌ No results found for {algorithm}")
    
    # Summary for NIPS
    print(f"\n📊 NIPS SUMMARY:")
    print("-" * 50)
    
    for algorithm in algorithms:
        cursor.execute('''
            SELECT AVG(relative_error) as avg_error,
                   MIN(relative_error) as min_error,
                   MAX(relative_error) as max_error
            FROM algorithm_comparison_p3 
            WHERE algorithm = ? AND p = 3 AND dataset = 'nips'
        ''', (algorithm,))
        
        avg_error, min_error, max_error = cursor.fetchone()
        print(f"{algorithm}: Avg={avg_error:.2f}%, Min={min_error:.2f}%, Max={max_error:.2f}%")
    
    conn.close()

if __name__ == "__main__":
    print_nips_results()
