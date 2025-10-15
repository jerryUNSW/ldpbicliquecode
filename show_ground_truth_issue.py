#!/usr/bin/env python3
"""
Script to demonstrate the ground truth overflow issue
"""

import sqlite3
import sys

def show_ground_truth_issue():
    # Connect to database
    conn = sqlite3.connect('../biclq_counts.db')
    cursor = conn.cursor()
    
    print("=== Ground Truth Overflow Analysis ===")
    print()
    
    # Get all datasets
    cursor.execute("SELECT DISTINCT dataset FROM pqbiclique_counts WHERE p = 2 ORDER BY dataset")
    datasets = [row[0] for row in cursor.fetchall()]
    
    max_64bit = 2**63 - 1  # 9223372036854775807
    
    for dataset in datasets:
        print(f"Dataset: {dataset}")
        print("Q\tGround Truth\t\t\tOverflowed?")
        print("-" * 60)
        
        cursor.execute("""
            SELECT q, count 
            FROM pqbiclique_counts 
            WHERE dataset = ? AND p = 2 AND q >= 4 
            ORDER BY q
        """, (dataset,))
        
        overflowed_count = 0
        for q, count in cursor.fetchall():
            is_overflowed = (count == max_64bit)
            if is_overflowed:
                overflowed_count += 1
            
            status = "YES" if is_overflowed else "NO"
            print(f"{q}\t{count:>20}\t{status}")
        
        if overflowed_count > 0:
            print(f"⚠️  {overflowed_count} Q values have overflowed ground truth!")
            print("   This means relative error calculations are meaningless for these Q values.")
        else:
            print("✅ No overflow issues detected.")
        
        print()
    
    conn.close()
    
    print("=== Summary ===")
    print("The issue is that when biclique counts exceed 64-bit integer limits,")
    print("they get stored as the maximum 64-bit integer value (9223372036854775807).")
    print("This makes relative error calculations meaningless because:")
    print("1. The 'ground truth' is artificially capped")
    print("2. Algorithm estimates are compared against wrong baseline")
    print("3. Relative errors appear much larger than they actually are")
    print()
    print("=== Solutions ===")
    print("1. Use arbitrary precision arithmetic for ground truth computation")
    print("2. Store ground truth as TEXT in database (scientific notation)")
    print("3. Update relative error calculation to handle large numbers")
    print("4. For Q values with overflowed ground truth, use absolute error instead")

if __name__ == "__main__":
    show_ground_truth_issue()
