#!/usr/bin/env python3
"""
Script to add jester2-swapped ADV, ADV+, and ADV++ results to the database.
"""

import sqlite3

def add_adv_results():
    # Connect to database
    conn = sqlite3.connect('algorithm_comparison_p3.db')
    cursor = conn.cursor()
    
    # ADV results for jester2-swapped (p=3, q=4-10)
    adv_results = [
        (4, 2.132051e+19, 2.1282142302957746548e+19, 0.01793, 0.01143),
        (5, 1.757e+23, 1.75165385309267381485568e+23, 0.02742, 0.01627),
        (6, 1.380e+27, 1.374676338568640554392879104e+27, 0.03517, 0.02053),
        (7, 9.726e+30, 9.688590147124980868032214597632e+30, 0.04202, 0.02441),
        (8, 6.098e+34, 6.0743130255762072724593175034855424e+34, 0.04847, 0.02819),
        (9, 3.420e+38, 3.40637528842903381231585298195821887488e+38, 0.05477, 0.03198),
        (10, 1.730e+42, 1.723585883467118343859341769038468506189824e+42, 0.06101, 0.03582)
    ]
    
    # ADV+ results for jester2-swapped (p=3, q=4-10)
    adv_plus_results = [
        (4, 2.132e+19, 2.1109794646425613532e+19, 0.01945, 0.008274),
        (5, 1.757e+23, 1.73112145095264967081984e+23, 0.02681, 0.009728),
        (6, 1.380e+27, 1.353784097832979721818734592e+27, 0.03311, 0.01188),
        (7, 9.726e+30, 9.506828451418673309988207722496e+30, 0.03902, 0.01408),
        (8, 6.098e+34, 5.9377160946629160018598338779676672e+34, 0.04483, 0.01631),
        (9, 3.420e+38, 3.31640223364390378256285604549073305600e+38, 0.05063, 0.01858),
        (10, 1.730e+42, 1.670924643349007878285909758310002687213568e+42, 0.05645, 0.02092)
    ]
    
    # ADV++ results for jester2-swapped (p=3, q=4-10)
    adv_plus_plus_results = [
        (4, 2.132e+19, 2.1117426766944560764e+19, 0.01109, 0.007783),
        (5, 1.757e+23, 1.73554107445074272305152e+23, 0.01423, 0.008824),
        (6, 1.380e+27, 1.359297096213463212805324800e+27, 0.01721, 0.01013),
        (7, 9.726e+30, 9.557319582752233360414509891584e+30, 0.02019, 0.01156),
        (8, 6.098e+34, 5.9762613769784840468240565715075072e+34, 0.02321, 0.01306),
        (9, 3.420e+38, 3.34195100491136796862999475605907439616e+38, 0.02628, 0.01459),
        (10, 1.730e+42, 1.685931787466407075871665773348994756378624e+42, 0.02939, 0.01614)
    ]
    
    # Insert ADV results
    print("Inserting ADV results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                'jester2-swapped',
                'ADV',
                q,
                ground_truth,
                estimate,
                relative_error,
                std_relative_error,
                -1,   # T_samples (no sampling)
                1.0,  # epsilon
                3     # p
            ))
            print(f"Inserted jester2-swapped ADV q={q}: rel_err={relative_error:.4f}")
        except Exception as e:
            print(f"Error inserting ADV q={q}: {e}")
    
    # Insert ADV+ results
    print("\nInserting ADV+ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                'jester2-swapped',
                'ADV+',
                q,
                ground_truth,
                estimate,
                relative_error,
                std_relative_error,
                -1,   # T_samples (no sampling)
                1.0,  # epsilon
                3     # p
            ))
            print(f"Inserted jester2-swapped ADV+ q={q}: rel_err={relative_error:.4f}")
        except Exception as e:
            print(f"Error inserting ADV+ q={q}: {e}")
    
    # Insert ADV++ results
    print("\nInserting ADV++ results...")
    for q, ground_truth, estimate, relative_error, std_relative_error in adv_plus_plus_results:
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO algorithm_comparison_p3 
                (dataset, algorithm, q, ground_truth, estimate, relative_error, std_relative_error, T_samples, epsilon, p)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                'jester2-swapped',
                'ADV++',
                q,
                ground_truth,
                estimate,
                relative_error,
                std_relative_error,
                -1,   # T_samples (no sampling)
                1.0,  # epsilon
                3     # p
            ))
            print(f"Inserted jester2-swapped ADV++ q={q}: rel_err={relative_error:.4f}")
        except Exception as e:
            print(f"Error inserting ADV++ q={q}: {e}")
    
    # Commit changes
    conn.commit()
    conn.close()
    print("\nSuccessfully added jester2-swapped ADV, ADV+, and ADV++ results to database")

if __name__ == "__main__":
    add_adv_results()
