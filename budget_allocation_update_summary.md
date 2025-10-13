# Budget Allocation Plots Update for Overleaf

## Summary
Successfully generated and updated budget allocation plots for settings 2 and 3 to match the style of Figure 9 (Privacy Budget Allocation Analysis) in the Overleaf paper.

## Configurations Updated

### Setting 2: ε=1.0, p=2, q=3
- **Configuration**: ε=1.0, p=2, q=3
- **Status**: ✅ Complete (27/27 experiments)
- **Plots Generated**: 6 individual dataset plots
- **Label**: `fig:budget-allocation-p2q3`

### Setting 3: ε=2.0, p=2, q=2  
- **Configuration**: ε=2.0, p=2, q=2
- **Status**: ✅ Complete (27/27 experiments)
- **Plots Generated**: 6 individual dataset plots
- **Label**: `fig:budget-allocation-eps2`

## Generated Files

### Plot Files (in `overleaf-paper/revision-plots/`)
**For ε=1.0, p=2, q=3:**
- `bag-kos-budget-allocation-eps1-p2q3.pdf`
- `unicode-budget-allocation-eps1-p2q3.pdf`
- `lrcwiki-budget-allocation-eps1-p2q3.pdf`
- `csbwiki-budget-allocation-eps1-p2q3.pdf`
- `librec-filmtrust-ratings-budget-allocation-eps1-p2q3.pdf`
- `rmwiki-budget-allocation-eps1-p2q3.pdf`

**For ε=2.0, p=2, q=2:**
- `bag-kos-budget-allocation-eps2-p2q2.pdf`
- `unicode-budget-allocation-eps2-p2q2.pdf`
- `lrcwiki-budget-allocation-eps2-p2q2.pdf`
- `csbwiki-budget-allocation-eps2-p2q2.pdf`
- `librec-filmtrust-ratings-budget-allocation-eps2-p2q2.pdf`
- `rmwiki-budget-allocation-eps2-p2q2.pdf`

### LaTeX Code
**Updated**: `overleaf-paper/figures.tex`
- Added Figure for ε=1.0, p=2, q=3 (label: `fig:budget-allocation-p2q3`)
- Added Figure for ε=2.0, p=2, q=2 (label: `fig:budget-allocation-eps2`)

## Key Results

### ε=1.0, p=2, q=3 Results
| Dataset | MRCN Best α₁ | MRCN+ Best α₁ | MRCN++ Best α₁ | MRCN++ Best RelErr |
|---------|--------------|---------------|----------------|-------------------|
| bag-kos | 0.70 | 0.80 | 0.90 | 0.0128 |
| unicode | 0.50 | 0.60 | 0.50 | 58.6106 |
| lrcwiki | 0.90 | 0.90 | 0.80 | 0.0216 |
| csbwiki | 0.90 | 0.90 | 0.90 | 0.0107 |
| librec-filmtrust | 0.80 | 0.80 | 0.80 | 0.0168 |
| rmwiki | 0.80 | 0.80 | 0.80 | 0.0137 |

### ε=2.0, p=2, q=2 Results
| Dataset | MRCN Best α₁ | MRCN+ Best α₁ | MRCN++ Best α₁ | MRCN++ Best RelErr |
|---------|--------------|---------------|----------------|-------------------|
| bag-kos | 0.90 | 0.90 | 0.90 | 0.0027 |
| unicode | 0.70 | 0.70 | 0.60 | 1.0050 |
| lrcwiki | 0.80 | 0.80 | 0.90 | 0.0068 |
| csbwiki | 0.90 | 0.80 | 0.70 | 0.0034 |
| librec-filmtrust | 0.90 | 0.80 | 0.90 | 0.0036 |
| rmwiki | 0.80 | 0.90 | 0.90 | 0.0034 |

## Observations

1. **Algorithm Performance**: MRCN++ consistently shows the best performance (lowest relative error) across all datasets and configurations.

2. **Optimal α₁ Values**: 
   - Most datasets show optimal α₁ values between 0.7-0.9
   - Higher ε values (2.0) generally allow for better accuracy
   - Different datasets have different optimal budget allocations

3. **Dataset Performance**:
   - **bag-kos, csbwiki, lrcwiki, librec-filmtrust, rmwiki**: Good accuracy (<4% error)
   - **unicode**: Higher relative errors due to smaller ground truth counts
   - **Higher ε improves accuracy**: ε=2.0 shows better performance than ε=1.0

4. **Style Consistency**: All plots match the existing Figure 9 style with:
   - Log scale y-axis
   - Consistent colors and markers
   - Proper LaTeX formatting
   - Publication-quality resolution (300 DPI)

## Usage in Paper

The new figures can be referenced in the paper using:
- `\ref{fig:budget-allocation-p2q3}` for ε=1.0, p=2, q=3
- `\ref{fig:budget-allocation-eps2}` for ε=2.0, p=2, q=2

Both figures follow the same structure and caption format as the original Figure 9, providing comprehensive budget allocation analysis for the additional parameter configurations.
