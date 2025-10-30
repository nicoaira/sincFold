# Vienna Integration Optimization Summary

## Baseline Performance

| Method | F1 Score | Notes |
|--------|----------|-------|
| **sincFold alone** | **0.7793** | Best overall |
| sincFold + Vienna (default) | 0.7509 | weight=1.0, log scaling, auto constraints |
| Vienna RNAfold alone | 0.4771 | Worst performer |

## Optimization Experiments

### Experiment 1: Weight and Scaling Method
**Parameters tested:**
- `vienna_weight`: [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]
- `vienna_linear`: [True, False]
- `vienna_temp`: [37.0]

**Key Findings:**
1. **Log scaling significantly outperforms linear scaling**
   - Log: mean F1 = 0.7443, best = 0.7519
   - Linear: mean F1 = 0.6913, best = 0.7495

2. **Best weight: 0.7** (F1 = 0.7519)
   - Only 0.13% improvement over default (weight=1.0)

3. **Weight sensitivity:**
   - Very low weights (0.1-0.3) hurt performance with linear scaling
   - High weights (2.0-5.0) perform well with log scaling
   - Optimal range: 0.5-3.0 with log scaling

### Experiment 2: Probability Threshold and Max Constraints
**Parameters tested:**
- `prob_threshold`: [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]
- `max_constraints`: [50, 100, 200, None (auto)]
- `vienna_weight`: [0.7, 3.0] (best from Exp 1)
- `vienna_linear`: [False] (log scaling only)

**Key Findings:**
1. **Probability threshold has minimal impact**
   - All thresholds (0.01-0.5) give similar mean F1 ≈ 0.739
   - Suggests that constraint quality matters less than quantity

2. **Max constraints is critical:**
   - 50 constraints: F1 = 0.708 ❌
   - 100 constraints: F1 = 0.747 ⚠️
   - 200 constraints: F1 = 0.752 ✓
   - auto (~550 for avg length 275): F1 = 0.752 ✓

3. **Sweet spot: 200 constraints** (sufficient without overhead)

## Best Configuration Found

```python
vienna_weight = 0.7
vienna_temp = 37.0
vienna_linear = False  # Use log scaling
prob_threshold = 0.3   # Doesn't matter much
max_constraints = 200  # Or auto
```

**Performance: F1 = 0.7519**

## Critical Analysis

### Why Vienna Integration Fails to Improve Performance

Despite extensive optimization, **sincFold+Vienna never outperforms sincFold alone** (0.7519 vs 0.7793).

**Hypotheses:**

1. **Model Incompatibility**
   - Vienna optimizes for minimum free energy (MFE) based on thermodynamics
   - sincFold optimizes for patterns learned from real structural data
   - These objectives are fundamentally different and sometimes contradictory

2. **Information Loss**
   - Vienna's MFE folding is deterministic and doesn't capture uncertainty
   - sincFold's probability matrix contains rich uncertainty information
   - Converting probabilities → soft constraints → MFE structure loses information

3. **Constraint Interference**
   - Vienna's thermodynamic model has strong priors (e.g., GC stacking)
   - sincFold's learned priors may conflict with thermodynamic priors
   - Soft constraints can't resolve these conflicts effectively

4. **Post-processing Gap**
   - sincFold uses sophisticated post-processing (conflict resolution)
   - Vienna applies constraints during folding, potentially creating incompatibilities

### Performance by RNA Family

Comparing sincFold vs sincFold+Vienna (best config):

| Family | sincFold F1 | sincFold+Vienna F1 | Difference |
|--------|-------------|-------------------|------------|
| tRNA | 0.9889 | 0.9777 | -0.0112 |
| 5S_rRNA | 0.9861 | 0.9313 | -0.0548 ⬇️ |
| 23S_rRNA | 0.9756 | 0.9278 | -0.0478 ⬇️ |
| SRP | 0.8985 | 0.8514 | -0.0471 ⬇️ |
| 16S_rRNA | 0.8699 | **0.8689** | -0.0010 (≈ same) |
| group_I_intron | 0.7867 | 0.7736 | -0.0131 |
| RNaseP | 0.7256 | 0.6871 | -0.0385 ⬇️ |
| telomerase | 0.3787 | 0.3718 | -0.0069 |
| tmRNA | 0.4035 | 0.3679 | -0.0356 ⬇️ |

**Vienna hurts performance across ALL families**, with largest impact on well-structured RNAs (rRNAs).

## Recommendations

### For Production Use

**Use sincFold alone** (F1 = 0.7793)

Vienna integration provides no benefit and consistently degrades performance by ~3.6% even with optimal parameters.

### For Research

If investigating Vienna integration further, consider:

1. **Selective application:** Only use Vienna on low-confidence regions
2. **Hard constraints:** Use sincFold's high-confidence predictions as hard constraints
3. **Ensemble methods:** Average sincFold and Vienna predictions
4. **Hybrid architecture:** Train sincFold to directly predict MFE-compatible structures

### Alternative Approaches

Instead of Vienna soft constraints:

1. **Use sincFold's native post-processing** (already optimal)
2. **Ensemble with other ML models** (not thermodynamic models)
3. **Fine-tune sincFold** on specific RNA families
4. **Add uncertainty estimation** to sincFold outputs

## Statistical Significance

All comparisons used **Wilcoxon signed-rank test** (paired, two-tailed):

- sincFold vs sincFold+Vienna: **p < 0.0001** (sincFold significantly better)
- sincFold vs Vienna alone: **p < 0.0001** (sincFold significantly better)
- sincFold+Vienna vs Vienna alone: **p < 0.0001** (sincFold+Vienna significantly better)

## Conclusion

**Vienna integration is not recommended.** Despite testing 66 parameter combinations:
- Best Vienna config: F1 = 0.7519
- sincFold alone: F1 = 0.7793
- **Gap: -0.0274 (-3.5%)**

The thermodynamic and learned approaches appear fundamentally incompatible when combined via soft constraints. sincFold's native predictions and post-processing are superior.

---

## Files Generated

1. `benchmark_vienna_results.csv` - Full benchmark with all 3 methods
2. `vienna_optimization_results_*.csv` - Weight & scaling optimization
3. `vienna_advanced_optimization_results.csv` - Threshold & constraint optimization
4. `optimize_vienna_params.py` - Optimization script #1
5. `optimize_vienna_advanced.py` - Optimization script #2

## Dataset

- **Source:** ArchiveII.512 test set (HuggingFace)
- **Sampling:** 10 sequences per family, random seed 42
- **Total sequences:** 90
- **Families:** 9 (16S/23S/5S rRNA, tRNA, RNaseP, SRP, group I intron, telomerase, tmRNA)
- **Metric:** Strict F1 (exact base pair matching)
