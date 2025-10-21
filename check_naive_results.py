#!/usr/bin/env python3
"""
Check if we have Naive results in the database for p=3.
"""

import sqlite3

def check_naive_results():
    """Check for Naive results in the database"""
    
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Check for Naive results with p=3
    cursor.execute('''
        SELECT algorithm, q, ground_truth, estimate, relative_error 
        FROM algorithm_comparison_p3 
        WHERE algorithm = 'Naive' AND p = 3
        ORDER BY q
    ''')
    
    results = cursor.fetchall()
    
    if results:
        print("✅ Naive results found for p=3:")
        print("=" * 80)
        print(f"{'Q':<3} {'Ground Truth':<15} {'Estimate':<15} {'Rel Error (%)':<15}")
        print("-" * 80)
        
        for algorithm, q, gt, estimate, rel_error in results:
            print(f"{q:<3} {gt:<15.2e} {estimate:<15.2e} {rel_error:<15.2f}")
    else:
        print("❌ No Naive results found for p=3")
    
    # Check what algorithms we have for p=3
    print(f"\n📊 All algorithms in database for p=3:")
    cursor.execute('''
        SELECT DISTINCT algorithm, COUNT(*) as count
        FROM algorithm_comparison_p3 
        WHERE p = 3
        GROUP BY algorithm
        ORDER BY algorithm
    ''')
    
    algorithms = cursor.fetchall()
    for algorithm, count in algorithms:
        print(f"  {algorithm}: {count} entries")
    
    # Check what datasets we have
    print(f"\n📊 All datasets in database:")
    cursor.execute('''
        SELECT DISTINCT dataset, COUNT(*) as count
        FROM algorithm_comparison_p3 
        GROUP BY dataset
        ORDER BY dataset
    ''')
    
    datasets = cursor.fetchall()
    for dataset, count in datasets:
        print(f"  {dataset}: {count} entries")
    
    conn.close()

if __name__ == "__main__":
    check_naive_results()
