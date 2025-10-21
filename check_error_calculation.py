#!/usr/bin/env python3
"""
Check how relative error is calculated and stored in the database.
"""

import sqlite3

def check_error_calculation():
    """Check the actual relative error calculation"""
    
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # Get one example from ADV Q=4
    cursor.execute('''
        SELECT q, ground_truth, estimate, relative_error 
        FROM algorithm_comparison_p3 
        WHERE algorithm = 'ADV' AND p = 3 AND dataset = 'nips' AND q = 4
    ''')
    
    result = cursor.fetchone()
    if result:
        q, gt, estimate, stored_error = result
        print(f"Database Example (ADV Q=4):")
        print(f"  Ground Truth: {gt:.2e}")
        print(f"  Estimate: {estimate:.2e}")
        print(f"  Stored Error: {stored_error}")
        
        # Calculate relative error manually
        raw_relative_error = abs(estimate - gt) / abs(gt)
        percentage_relative_error = raw_relative_error * 100
        
        print(f"\nManual Calculation:")
        print(f"  Raw relative error: {raw_relative_error:.6f}")
        print(f"  Percentage relative error: {percentage_relative_error:.2f}%")
        print(f"  Stored value matches raw: {abs(stored_error - raw_relative_error) < 1e-6}")
        print(f"  Stored value matches percentage: {abs(stored_error - percentage_relative_error) < 1e-6}")
    
    # Check a few more examples
    print(f"\n" + "="*60)
    print("Checking multiple examples:")
    print("="*60)
    
    cursor.execute('''
        SELECT algorithm, q, ground_truth, estimate, relative_error 
        FROM algorithm_comparison_p3 
        WHERE algorithm IN ('ADV', 'ADV+') AND p = 3 AND dataset = 'nips' AND q IN (4, 5)
        ORDER BY algorithm, q
    ''')
    
    results = cursor.fetchall()
    for algorithm, q, gt, estimate, stored_error in results:
        raw_error = abs(estimate - gt) / abs(gt)
        print(f"{algorithm} Q={q}: GT={gt:.2e}, Est={estimate:.2e}")
        print(f"  Raw error: {raw_error:.6f}, Stored: {stored_error:.6f}, Match: {abs(stored_error - raw_error) < 1e-6}")
        print()
    
    conn.close()

if __name__ == "__main__":
    check_error_calculation()
