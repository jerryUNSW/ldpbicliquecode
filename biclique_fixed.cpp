// This is a patch file showing the exact changes needed to fix the overflow issue
// Apply these changes to biclique.cpp

/*
CHANGE 1: Replace the overflow warning sections in biclique.cpp

FIND (around line 2732):
        // Also store as unsigned long long for backward compatibility (with truncation warning)
        if (count_value > static_cast<long double>(ULLONG_MAX)) {
            std::cout << "WARNING: Ground truth exceeds unsigned long long limit for dataset: " << dataset 
                      << ", p = " << P___ << ", q = " << K___ << std::endl;
            std::cout << "Original value: " << std::scientific << std::setprecision(6) << count_value << std::endl;
            std::cout << "Will use long double precision for relative error calculation" << std::endl;
            real = ULLONG_MAX; // Set to max for backward compatibility
        } else {
            real = static_cast<unsigned long long>(count_value);
        }

REPLACE WITH:
        // Use long double consistently for ground truth
        real = count_value;  // Store as long double (no truncation)
*/

/*
CHANGE 2: Replace the second overflow warning section (around line 2760):

FIND:
        // Also store as unsigned long long for backward compatibility (with truncation warning)
        if (computed_count > static_cast<long double>(ULLONG_MAX)) {
            std::cout << "WARNING: Computed count exceeds unsigned long long limit for dataset: " << dataset_to_find 
                      << ", p = " << P___ << ", q = " << K___ << std::endl;
            std::cout << "Original value: " << std::scientific << std::setprecision(6) << computed_count << std::endl;
            std::cout << "Will use long double precision for relative error calculation" << std::endl;
            real = ULLONG_MAX; // Set to max for backward compatibility
        } else {
            real = static_cast<unsigned long long>(computed_count);
        }

REPLACE WITH:
        // Use long double consistently for ground truth
        real = computed_count;  // Store as long double (no truncation)
*/

/*
CHANGE 3: Update the variable declarations

FIND (around line 35):
extern unsigned long long int real; 
extern long double real_ld;

REPLACE WITH:
extern long double real;  // Use long double instead of unsigned long long
extern long double real_ld;
*/

/*
CHANGE 4: Update the variable declarations in main.cpp and test files

FIND:
unsigned long long int real;
long double real_ld = 0.0;

REPLACE WITH:
long double real = 0.0;  // Use long double instead of unsigned long long
long double real_ld = 0.0;
*/

/*
CHANGE 5: Update the get_ground_truth function in test_p3_batch_with_ground_truth.cpp

FIND (around line 49):
    return (real_ld > 0) ? real_ld : static_cast<long double>(real);

REPLACE WITH:
    return real;  // Now both real and real_ld are long double, so just use real
*/
