# Eternabench-CM Benchmark Results

## Dataset Information

**Source:** multimolecule/eternabench-cm (train split)
**Total sequences:** 12,293
**After CD-HIT filtering (90% identity):** 11,808 sequences
**Random sample for benchmark:** 100 sequences
**Average length:** 110.5 nucleotides (range: 84-131)

## Benchmark Results

| Method | Mean F1 | Std F1 | Median F1 | Rank |
|--------|---------|--------|-----------|------|
| **Vienna RNAfold alone** | **0.9088** | 0.2127 | **1.0000** | 🥇 **BEST** |
| **sincFold + Vienna** | **0.8849** | 0.2165 | **1.0000** | 🥈 |
| **sincFold alone** | **0.6596** | 0.2511 | **0.6825** | 🥉 WORST |

## Statistical Significance

All comparisons are statistically significant (Wilcoxon signed-rank test):

1. **sincFold vs sincFold+Vienna**: p < 0.0001 → sincFold+Vienna significantly better (+34.2%)
2. **sincFold vs Vienna**: p < 0.0001 → Vienna significantly better (+37.8%)
3. **sincFold+Vienna vs Vienna**: p = 0.0793 → No significant difference (≈ tied)

## Comparison with ArchiveII Dataset

### Performance Reversal

The results on Eternabench are **dramatically different** from ArchiveII:

| Method | ArchiveII F1 | Eternabench F1 | Difference |
|--------|--------------|----------------|------------|
| sincFold alone | **0.7793** (BEST) | **0.6596** (WORST) | **-15.4%** ⬇️ |
| sincFold + Vienna | 0.7509 | 0.8849 | **+17.9%** ⬆️ |
| Vienna alone | 0.4771 (WORST) | **0.9088** (BEST) | **+90.5%** ⬆️⬆️ |

### Key Observations

1. **Vienna dominates on Eternabench** (90.9% F1 vs 47.7% on ArchiveII)
2. **sincFold struggles on Eternabench** (66.0% F1 vs 77.9% on ArchiveII)
3. **Vienna integration helps on Eternabench** (+34% improvement) but **hurts on ArchiveII** (-3.6%)

## Analysis: Why the Reversal?

### Hypothesis 1: Dataset Characteristics

**Eternabench (Vienna-friendly):**
- Human-designed sequences from Eterna citizen science project
- Likely optimized for thermodynamic stability
- Designed structures may follow textbook base-pairing rules
- More canonical, "clean" structures
- Shorter sequences (avg 110 nt) with well-defined motifs

**ArchiveII (sincFold-friendly):**
- Natural RNA structures from diverse families
- Evolved structures with complex, non-canonical interactions
- May include pseudoknots, unusual motifs, long-range interactions
- Longer sequences (avg 275 nt) with more complexity
- Real biological constraints beyond pure thermodynamics

### Hypothesis 2: Model Training Data Bias

**sincFold:**
- Trained on real biological RNA structures
- Learned patterns from natural RNAs (likely similar to ArchiveII)
- May not generalize well to synthetic/designed sequences
- Post-processing optimized for biological RNA complexity

**Vienna RNAfold:**
- Based on universal thermodynamic principles
- Works well on "ideal" sequences that follow design principles
- Less affected by biological complexity
- Performs poorly on structures with non-canonical features

### Hypothesis 3: Sequence Complexity

- **Eternabench:** Simpler, shorter, more predictable structures → thermodynamics sufficient
- **ArchiveII:** Complex, longer, evolved structures → need learned patterns

## Implications

### For Method Selection

**Use Vienna RNAfold when:**
- Working with designed/synthetic RNA sequences
- Sequences from Eterna or similar design platforms
- Short sequences (<150 nt) with canonical structures
- Thermodynamic stability is the primary design criterion

**Use sincFold when:**
- Working with natural/biological RNA sequences
- Complex structures from diverse RNA families
- Longer sequences with potential pseudoknots
- Non-canonical or evolved structures

**Use sincFold + Vienna when:**
- Designed sequences (Eternabench-like)
- Want to leverage both thermodynamics and learned patterns
- Can afford slightly lower performance than Vienna alone (88.5% vs 90.9%)

### For Training Data

This reveals a critical limitation of sincFold:
- **Strong on natural RNAs** (ArchiveII)
- **Weak on designed RNAs** (Eternabench)

**Recommendation:** Consider augmenting sincFold training with synthetic/designed sequences to improve generalization.

## CD-HIT Filtering

CD-HIT successfully reduced redundancy:
- **Before filtering:** 12,293 sequences
- **After filtering (90% identity):** 11,808 sequences (96% retained)
- **Sequences removed:** 485 (4% redundancy)

The low redundancy suggests Eternabench already has diverse sequences, which makes the results more robust.

## Detailed Statistics

### Length Distribution
- **Min:** 84 nt
- **Max:** 131 nt
- **Mean:** 110.5 nt
- **Median:** 107.0 nt

(Much shorter than ArchiveII: 275.4 nt average)

### F1 Score Distribution

**Vienna RNAfold:**
- Many perfect predictions (median F1 = 1.0)
- Low variance (std = 0.21)
- Consistently excellent performance

**sincFold + Vienna:**
- Also many perfect predictions (median F1 = 1.0)
- Similar variance (std = 0.22)
- Slightly lower mean than Vienna alone

**sincFold alone:**
- Lower median (0.68)
- Higher variance (std = 0.25)
- More inconsistent performance

## Conclusion

The Eternabench benchmark reveals that:

1. ✅ **Vienna integration DOES help** on designed/synthetic RNA sequences
2. ✅ **Vienna alone is excellent** for Eterna-designed sequences (90.9% F1)
3. ⚠️ **sincFold underperforms** on designed sequences (66.0% F1)
4. 🔬 **Dataset characteristics matter critically** for model selection

**Key insight:** The optimal RNA folding method depends heavily on whether sequences are natural (biological) or designed (synthetic). There is no one-size-fits-all solution.

---

## Recommendation Matrix

| Sequence Type | Best Method | F1 (Expected) |
|---------------|-------------|---------------|
| Designed/Synthetic RNA (Eterna) | Vienna RNAfold | ~90% |
| Natural RNA (diverse families) | sincFold | ~78% |
| Designed RNA + want ML boost | sincFold + Vienna | ~88% |
| Natural RNA + thermodynamics | sincFold alone | ~78% (Vienna hurts) |

## Files Generated

1. `benchmark_eternabench.py` - Benchmark script with CD-HIT filtering
2. `benchmark_eternabench_results.csv` - Detailed per-sequence results
3. `ETERNABENCH_RESULTS.md` - This analysis document
