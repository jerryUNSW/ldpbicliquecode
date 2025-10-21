#!/usr/bin/env python3
"""
Print Naive, ADV, and ADV+ results from the database.
"""

import sqlite3

def print_algorithm_results():
    """Print results for Naive, ADV, and ADV+ algorithms"""
    
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Get results for each algorithm
    algorithms = ['Naive', 'ADV', 'ADV+']
    
    for algorithm in algorithms:
        print(f"\n{'='*80}")
        print(f"📊 {algorithm} ALGORITHM RESULTS")
        print(f"{'='*80}")
        
        cursor.execute('''
            SELECT dataset, q, ground_truth, estimate, relative_error 
            FROM algorithm_comparison_p3 
            WHERE algorithm = ? AND p = 3
            ORDER BY dataset, q
        ''', (algorithm,))
        
        results = cursor.fetchall()
        
        if results:
            print(f"{'Dataset':<25} {'Q':<3} {'Ground Truth':<15} {'Estimate':<15} {'Rel Error (%)':<15}")
            print("-" * 100)
            
            current_dataset = None
            for dataset, q, gt, estimate, rel_error in results:
                if dataset != current_dataset:
                    print(f"\n{dataset}:")
                    current_dataset = dataset
                
                print(f"{'':<25} {q:<3} {gt:<15.2e} {estimate:<15.2e} {rel_error:<15.2f}")
        else:
            print(f"❌ No results found for {algorithm}")
    
    # Summary statistics
    print(f"\n{'='*80}")
    print(f"📊 SUMMARY STATISTICS")
    print(f"{'='*80}")
    
    for algorithm in algorithms:
        cursor.execute('''
            SELECT COUNT(*) as count, 
                   AVG(relative_error) as avg_error,
                   MIN(relative_error) as min_error,
                   MAX(relative_error) as max_error
            FROM algorithm_comparison_p3 
            WHERE algorithm = ? AND p = 3
        ''', (algorithm,))
        
        count, avg_error, min_error, max_error = cursor.fetchone()
        print(f"{algorithm}: {count} entries, Avg Error: {avg_error:.2f}%, Min: {min_error:.2f}%, Max: {max_error:.2f}%")
    
    conn.close()

if __name__ == "__main__":
    print_algorithm_results()
