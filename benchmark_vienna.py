import pandas as pd
import numpy as np
import torch as tr
from tqdm import tqdm
import warnings
import tempfile
import os
warnings.filterwarnings('ignore')

# Import sincfold modules
from sincfold.model import sincfold
from sincfold.metrics import f1_shift, f1_strict
from sincfold.utils import dot2bp, bp2matrix, fold_with_vienna, VIENNA_AVAILABLE
from sincfold.dataset import SeqDataset, pad_batch
from torch.utils.data import DataLoader

if not VIENNA_AVAILABLE:
    raise ImportError(
        "ViennaRNA is not installed. Please install it to run this benchmark. "
        "Installation: conda install -c bioconda viennarna OR pip install ViennaRNA"
    )

import RNA

def fold_vienna_only(sequence):
    """
    Fold RNA sequence using only ViennaRNA (no sincFold).

    Args:
        sequence (str): RNA sequence

    Returns:
        list: List of 1-indexed base pairs [[i, j], ...]
    """
    sequence_rna = sequence.upper().replace('T', 'U')
    fc = RNA.fold_compound(sequence_rna)
    structure, mfe = fc.mfe()
    base_pairs = dot2bp(structure)

    if base_pairs is False:
        base_pairs = []

    return base_pairs

def fold_sincfold_only(model, sequence, seq_id, config):
    """
    Fold RNA sequence using only sincFold.

    Args:
        model: sincfold model instance
        sequence (str): RNA sequence
        seq_id (str): Sequence ID
        config (dict): Config dictionary

    Returns:
        list: List of 1-indexed base pairs [[i, j], ...]
    """
    # Create temporary CSV file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
        tmp.write('id,sequence\n')
        tmp.write(f'{seq_id},{sequence}\n')
        tmp_file = tmp.name

    try:
        # Create DataLoader
        pred_loader = DataLoader(
            SeqDataset(tmp_file, for_prediction=True, **config),
            batch_size=1,
            shuffle=False,
            num_workers=0,
            collate_fn=pad_batch,
        )

        # Predict
        predictions, _ = model.pred(pred_loader, logits=False, use_vienna=False)
        base_pairs = predictions.iloc[0]['base_pairs']
    finally:
        os.remove(tmp_file)

    return base_pairs

def fold_sincfold_vienna(model, sequence, seq_id, config, vienna_weight=1.0, vienna_temp=37.0, vienna_linear=False):
    """
    Fold RNA sequence using sincFold + ViennaRNA integration.

    Args:
        model: sincfold model instance
        sequence (str): RNA sequence
        seq_id (str): Sequence ID
        config (dict): Config dictionary
        vienna_weight (float): Weight for Vienna constraints
        vienna_temp (float): Temperature for Vienna
        vienna_linear (bool): Use linear scaling

    Returns:
        list: List of 1-indexed base pairs [[i, j], ...]
    """
    # Create temporary CSV file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
        tmp.write('id,sequence\n')
        tmp.write(f'{seq_id},{sequence}\n')
        tmp_file = tmp.name

    try:
        # Create DataLoader
        pred_loader = DataLoader(
            SeqDataset(tmp_file, for_prediction=True, **config),
            batch_size=1,
            shuffle=False,
            num_workers=0,
            collate_fn=pad_batch,
        )

        # Predict with Vienna
        predictions, _ = model.pred(
            pred_loader,
            logits=False,
            use_vienna=True,
            vienna_weight=vienna_weight,
            vienna_temp=vienna_temp,
            vienna_linear=vienna_linear
        )
        base_pairs = predictions.iloc[0]['base_pairs']
    finally:
        os.remove(tmp_file)

    return base_pairs

def calculate_f1(pred_bp, ref_bp):
    """
    Calculate F1 score with exact matching (strict).

    Precision = TP / (TP + FP)
    Recall = TP / (TP + FN)
    F1 = 2 * (Precision * Recall) / (Precision + Recall)

    Args:
        pred_bp (list): Predicted base pairs (1-indexed)
        ref_bp (list): Reference base pairs (1-indexed)

    Returns:
        tuple: (recall, precision, f1)
    """
    return f1_strict(ref_bp, pred_bp)

def main():
    print("=" * 80)
    print("RNA FOLDING BENCHMARK: sincFold vs sincFold+Vienna vs Vienna")
    print("=" * 80)
    print()
    print("F1 Score Definition (Strict Matching):")
    print("  - TP: Exact base pair matches between prediction and reference")
    print("  - FP: Predicted base pairs not in reference")
    print("  - FN: Reference base pairs not predicted")
    print("  - Precision = TP / (TP + FP)")
    print("  - Recall = TP / (TP + FN)")
    print("  - F1 = 2 * (Precision * Recall) / (Precision + Recall)")
    print()

    # Set random seed for reproducibility
    RANDOM_SEED = 42
    np.random.seed(RANDOM_SEED)
    tr.manual_seed(RANDOM_SEED)
    print(f"Random seed: {RANDOM_SEED}")
    print()

    # Load dataset
    print("Loading dataset from HuggingFace...")
    df = pd.read_parquet("hf://datasets/multimolecule/archiveii.512/test.parquet")
    print(f"Total sequences in dataset: {len(df)}")
    print(f"Unique families: {df['family'].nunique()}")
    print()

    # Sample 10 sequences from each family
    print("Sampling 10 sequences from each family...")
    sampled_df = df.groupby('family', group_keys=False).apply(
        lambda x: x.sample(n=min(10, len(x)), random_state=RANDOM_SEED)
    ).reset_index(drop=True)
    print(f"Total sequences sampled: {len(sampled_df)}")
    print(f"Families: {sorted(sampled_df['family'].unique())}")
    print()

    # Initialize sincfold model
    print("Loading sincFold model...")
    device = 'cuda' if tr.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")

    config = {
        "device": device,
        "batch_size": 1,
        "verbose": False
    }

    model = sincfold(pretrained=True, **config)
    print("Model loaded successfully!")
    print()

    # Prepare results storage
    results = []

    # Vienna parameters (matching defaults in parser)
    vienna_weight = 1.0
    vienna_temp = 37.0
    vienna_linear = False

    print("=" * 80)
    print("Starting benchmark...")
    print("=" * 80)
    print()

    # Process each sequence
    for idx, row in tqdm(sampled_df.iterrows(), total=len(sampled_df), desc="Processing sequences"):
        seq_id = row['id']
        sequence = row['sequence']
        ref_structure = row['secondary_structure']
        family = row['family']

        # Parse reference structure
        ref_bp = dot2bp(ref_structure)
        if ref_bp is False:
            print(f"Warning: Could not parse reference structure for {seq_id}, skipping...")
            continue

        try:
            # Method 1: sincFold only
            pred_bp_sincfold = fold_sincfold_only(model, sequence, seq_id, config)
            recall_sf, prec_sf, f1_sf = calculate_f1(pred_bp_sincfold, ref_bp)

            # Method 2: sincFold + Vienna
            pred_bp_sincfold_vienna = fold_sincfold_vienna(
                model, sequence, seq_id, config,
                vienna_weight=vienna_weight,
                vienna_temp=vienna_temp,
                vienna_linear=vienna_linear
            )
            recall_sfv, prec_sfv, f1_sfv = calculate_f1(pred_bp_sincfold_vienna, ref_bp)

            # Method 3: Vienna only
            pred_bp_vienna = fold_vienna_only(sequence)
            recall_v, prec_v, f1_v = calculate_f1(pred_bp_vienna, ref_bp)

            # Store results
            results.append({
                'id': seq_id,
                'family': family,
                'length': len(sequence),
                'ref_pairs': len(ref_bp),
                # sincFold
                'sincfold_f1': f1_sf,
                'sincfold_recall': recall_sf,
                'sincfold_precision': prec_sf,
                'sincfold_pred_pairs': len(pred_bp_sincfold),
                # sincFold + Vienna
                'sincfold_vienna_f1': f1_sfv,
                'sincfold_vienna_recall': recall_sfv,
                'sincfold_vienna_precision': prec_sfv,
                'sincfold_vienna_pred_pairs': len(pred_bp_sincfold_vienna),
                # Vienna only
                'vienna_f1': f1_v,
                'vienna_recall': recall_v,
                'vienna_precision': prec_v,
                'vienna_pred_pairs': len(pred_bp_vienna),
            })

        except Exception as e:
            print(f"\nError processing {seq_id}: {e}")
            continue

    print()
    print("=" * 80)
    print("Benchmark completed!")
    print("=" * 80)
    print()

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Save detailed results
    output_file = 'benchmark_vienna_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"Detailed results saved to: {output_file}")
    print()

    # Print summary statistics
    print("=" * 80)
    print("OVERALL RESULTS")
    print("=" * 80)
    print()

    methods = [
        ('sincFold', 'sincfold_f1'),
        ('sincFold + Vienna', 'sincfold_vienna_f1'),
        ('Vienna only', 'vienna_f1')
    ]

    print(f"{'Method':<25} {'Mean F1':<12} {'Std F1':<12} {'Median F1':<12}")
    print("-" * 80)
    for method_name, col_name in methods:
        mean_f1 = results_df[col_name].mean()
        std_f1 = results_df[col_name].std()
        median_f1 = results_df[col_name].median()
        print(f"{method_name:<25} {mean_f1:<12.4f} {std_f1:<12.4f} {median_f1:<12.4f}")

    print()
    print("=" * 80)
    print("RESULTS BY FAMILY")
    print("=" * 80)
    print()

    # Group by family and calculate mean F1
    family_results = results_df.groupby('family').agg({
        'sincfold_f1': 'mean',
        'sincfold_vienna_f1': 'mean',
        'vienna_f1': 'mean',
        'length': 'mean'
    }).round(4)

    family_results = family_results.rename(columns={
        'sincfold_f1': 'sincFold',
        'sincfold_vienna_f1': 'sincFold+Vienna',
        'vienna_f1': 'Vienna',
        'length': 'Avg Length'
    })

    print(family_results.to_string())
    print()

    # Statistical comparison
    print("=" * 80)
    print("STATISTICAL COMPARISON")
    print("=" * 80)
    print()

    # Wilcoxon signed-rank test (paired comparison)
    from scipy.stats import wilcoxon

    # sincFold vs sincFold+Vienna
    try:
        stat1, p1 = wilcoxon(results_df['sincfold_f1'], results_df['sincfold_vienna_f1'])
        print(f"sincFold vs sincFold+Vienna:")
        print(f"  Wilcoxon statistic: {stat1:.4f}, p-value: {p1:.4f}")
        if p1 < 0.05:
            winner = 'sincFold+Vienna' if results_df['sincfold_vienna_f1'].mean() > results_df['sincfold_f1'].mean() else 'sincFold'
            print(f"  Result: Statistically significant difference (p < 0.05), {winner} performs better")
        else:
            print(f"  Result: No statistically significant difference (p >= 0.05)")
        print()
    except Exception as e:
        print(f"Could not perform sincFold vs sincFold+Vienna comparison: {e}")
        print()

    # sincFold vs Vienna
    try:
        stat2, p2 = wilcoxon(results_df['sincfold_f1'], results_df['vienna_f1'])
        print(f"sincFold vs Vienna:")
        print(f"  Wilcoxon statistic: {stat2:.4f}, p-value: {p2:.4f}")
        if p2 < 0.05:
            winner = 'Vienna' if results_df['vienna_f1'].mean() > results_df['sincfold_f1'].mean() else 'sincFold'
            print(f"  Result: Statistically significant difference (p < 0.05), {winner} performs better")
        else:
            print(f"  Result: No statistically significant difference (p >= 0.05)")
        print()
    except Exception as e:
        print(f"Could not perform sincFold vs Vienna comparison: {e}")
        print()

    # sincFold+Vienna vs Vienna
    try:
        stat3, p3 = wilcoxon(results_df['sincfold_vienna_f1'], results_df['vienna_f1'])
        print(f"sincFold+Vienna vs Vienna:")
        print(f"  Wilcoxon statistic: {stat3:.4f}, p-value: {p3:.4f}")
        if p3 < 0.05:
            winner = 'Vienna' if results_df['vienna_f1'].mean() > results_df['sincfold_vienna_f1'].mean() else 'sincFold+Vienna'
            print(f"  Result: Statistically significant difference (p < 0.05), {winner} performs better")
        else:
            print(f"  Result: No statistically significant difference (p >= 0.05)")
        print()
    except Exception as e:
        print(f"Could not perform sincFold+Vienna vs Vienna comparison: {e}")
        print()

    print("=" * 80)
    print("Benchmark complete! Check", output_file, "for detailed results.")
    print("=" * 80)

if __name__ == "__main__":
    main()
