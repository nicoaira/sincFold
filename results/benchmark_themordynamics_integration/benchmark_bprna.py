import pandas as pd
import numpy as np
import torch as tr
from tqdm import tqdm
import warnings
import tempfile
import os
import subprocess
from datasets import load_dataset
warnings.filterwarnings('ignore')

# Import sincfold modules
from sincfold.model import sincfold
from sincfold.metrics import f1_strict
from sincfold.utils import dot2bp, fold_with_vienna, VIENNA_AVAILABLE
from sincfold.dataset import SeqDataset, pad_batch
from torch.utils.data import DataLoader

if VIENNA_AVAILABLE:
    import RNA
else:
    print("Warning: ViennaRNA not available, will only test sincFold")

def run_cdhit_filtering(sequences, ids, identity_threshold=0.9, output_prefix="cdhit_temp"):
    """
    Filter sequences using CD-HIT at specified identity threshold.

    Args:
        sequences (list): List of sequences
        ids (list): List of sequence IDs
        identity_threshold (float): Identity threshold (0.0-1.0)
        output_prefix (str): Prefix for temporary files

    Returns:
        list: Indices of representative sequences
    """
    # Write sequences to FASTA file
    fasta_file = f"{output_prefix}.fasta"
    with open(fasta_file, 'w') as f:
        for i, (seq_id, seq) in enumerate(zip(ids, sequences)):
            f.write(f">{seq_id}\n{seq}\n")

    # Run CD-HIT
    output_file = f"{output_prefix}_filtered.fasta"
    cluster_file = f"{output_file}.clstr"

    cmd = [
        "cd-hit-est",
        "-i", fasta_file,
        "-o", output_file,
        "-c", str(identity_threshold),  # Identity threshold
        "-n", "5",  # Word length (5 for ~90% identity)
        "-M", "0",  # Unlimited memory
        "-T", "0",  # Use all available threads
        "-d", "0"   # Full sequence name in output
    ]

    print(f"Running CD-HIT with {identity_threshold*100:.0f}% identity threshold...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"CD-HIT error: {result.stderr}")
        raise RuntimeError("CD-HIT failed")

    # Parse output to get representative sequence IDs
    representative_ids = []
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip()
                representative_ids.append(seq_id)

    # Get indices of representative sequences
    id_to_idx = {seq_id: i for i, seq_id in enumerate(ids)}
    representative_indices = [id_to_idx[rep_id] for rep_id in representative_ids if rep_id in id_to_idx]

    # Clean up temporary files
    for f in [fasta_file, output_file, cluster_file]:
        if os.path.exists(f):
            os.remove(f)

    print(f"CD-HIT filtering: {len(sequences)} → {len(representative_indices)} sequences")

    return representative_indices

def fold_vienna_only(sequence):
    """Fold with Vienna RNAfold alone."""
    if not VIENNA_AVAILABLE:
        return []

    sequence_rna = sequence.upper().replace('T', 'U')
    fc = RNA.fold_compound(sequence_rna)
    structure, mfe = fc.mfe()
    base_pairs = dot2bp(structure)

    if base_pairs is False:
        base_pairs = []

    return base_pairs

def fold_sincfold_only(model, sequence, seq_id, config):
    """Fold with sincFold alone."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
        tmp.write('id,sequence\n')
        tmp.write(f'{seq_id},{sequence}\n')
        tmp_file = tmp.name

    try:
        pred_loader = DataLoader(
            SeqDataset(tmp_file, for_prediction=True, **config),
            batch_size=1,
            shuffle=False,
            num_workers=0,
            collate_fn=pad_batch,
        )

        predictions, _ = model.pred(pred_loader, logits=False, use_vienna=False)
        base_pairs = predictions.iloc[0]['base_pairs']
    finally:
        os.remove(tmp_file)

    return base_pairs

def fold_sincfold_vienna(model, sequence, seq_id, config, vienna_weight=0.7,
                         vienna_temp=37.0, vienna_linear=False):
    """Fold with sincFold + Vienna integration (optimized parameters)."""
    if not VIENNA_AVAILABLE:
        return fold_sincfold_only(model, sequence, seq_id, config)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
        tmp.write('id,sequence\n')
        tmp.write(f'{seq_id},{sequence}\n')
        tmp_file = tmp.name

    try:
        pred_loader = DataLoader(
            SeqDataset(tmp_file, for_prediction=True, **config),
            batch_size=1,
            shuffle=False,
            num_workers=0,
            collate_fn=pad_batch,
        )

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
    """Calculate strict F1 score."""
    return f1_strict(ref_bp, pred_bp)

def main():
    print("=" * 80)
    print("BPRNA-NEW BENCHMARK")
    print("=" * 80)
    print()

    # Set random seed
    RANDOM_SEED = 42
    np.random.seed(RANDOM_SEED)
    tr.manual_seed(RANDOM_SEED)
    print(f"Random seed: {RANDOM_SEED}")
    print()

    # Load dataset
    print("Loading bprna-new test split...")
    ds = load_dataset("multimolecule/bprna-new")
    test = ds['test']
    print(f"Total sequences in test split: {len(test)}")
    print()

    # Convert to lists for processing
    sequences = test['sequence']
    structures = test['secondary_structure']
    ids = test['id']

    # Filter with CD-HIT at 90% identity
    print("Filtering sequences with CD-HIT (90% identity threshold)...")
    representative_indices = run_cdhit_filtering(
        sequences,
        ids,
        identity_threshold=0.9,
        output_prefix="/tmp/bprna_cdhit"
    )
    print()

    # Get filtered sequences
    filtered_sequences = [sequences[i] for i in representative_indices]
    filtered_structures = [structures[i] for i in representative_indices]
    filtered_ids = [ids[i] for i in representative_indices]

    # Randomly sample 100 sequences
    n_sample = min(100, len(filtered_sequences))
    print(f"Randomly sampling {n_sample} sequences from {len(filtered_sequences)} filtered sequences...")

    sample_indices = np.random.choice(len(filtered_sequences), n_sample, replace=False)
    sample_sequences = [filtered_sequences[i] for i in sample_indices]
    sample_structures = [filtered_structures[i] for i in sample_indices]
    sample_ids = [filtered_ids[i] for i in sample_indices]

    print(f"Sampled {len(sample_sequences)} sequences")
    print(f"Average length: {np.mean([len(s) for s in sample_sequences]):.1f}")
    print()

    # Initialize sincFold model
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

    # Vienna parameters (optimized from previous experiments)
    vienna_weight = 0.7
    vienna_temp = 37.0
    vienna_linear = False

    print("=" * 80)
    print("Starting benchmark...")
    print("=" * 80)
    print()

    # Process each sequence
    for seq_id, sequence, ref_structure in tqdm(
        zip(sample_ids, sample_sequences, sample_structures),
        total=len(sample_sequences),
        desc="Processing sequences"
    ):
        # Parse reference structure
        ref_bp = dot2bp(ref_structure)
        if ref_bp is False:
            print(f"Warning: Could not parse reference structure for {seq_id}, skipping...")
            continue

        try:
            # Method 1: sincFold only
            pred_bp_sincfold = fold_sincfold_only(model, sequence, seq_id, config)
            recall_sf, prec_sf, f1_sf = calculate_f1(pred_bp_sincfold, ref_bp)

            # Method 2: sincFold + Vienna (if available)
            if VIENNA_AVAILABLE:
                pred_bp_sincfold_vienna = fold_sincfold_vienna(
                    model, sequence, seq_id, config,
                    vienna_weight=vienna_weight,
                    vienna_temp=vienna_temp,
                    vienna_linear=vienna_linear
                )
                recall_sfv, prec_sfv, f1_sfv = calculate_f1(pred_bp_sincfold_vienna, ref_bp)
            else:
                pred_bp_sincfold_vienna = []
                recall_sfv, prec_sfv, f1_sfv = 0.0, 0.0, 0.0

            # Method 3: Vienna only (if available)
            if VIENNA_AVAILABLE:
                pred_bp_vienna = fold_vienna_only(sequence)
                recall_v, prec_v, f1_v = calculate_f1(pred_bp_vienna, ref_bp)
            else:
                pred_bp_vienna = []
                recall_v, prec_v, f1_v = 0.0, 0.0, 0.0

            # Store results
            result = {
                'id': seq_id,
                'length': len(sequence),
                'ref_pairs': len(ref_bp),
                # sincFold
                'sincfold_f1': f1_sf,
                'sincfold_recall': recall_sf,
                'sincfold_precision': prec_sf,
                'sincfold_pred_pairs': len(pred_bp_sincfold),
            }

            if VIENNA_AVAILABLE:
                result.update({
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

            results.append(result)

        except Exception as e:
            print(f"\nError processing {seq_id}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print()
    print("=" * 80)
    print("Benchmark completed!")
    print("=" * 80)
    print()

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Save detailed results
    output_file = 'benchmark_bprna_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"Detailed results saved to: {output_file}")
    print()

    # Print summary statistics
    print("=" * 80)
    print("OVERALL RESULTS")
    print("=" * 80)
    print()

    if VIENNA_AVAILABLE:
        methods = [
            ('sincFold', 'sincfold_f1'),
            ('sincFold + Vienna', 'sincfold_vienna_f1'),
            ('Vienna only', 'vienna_f1')
        ]
    else:
        methods = [
            ('sincFold', 'sincfold_f1'),
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
    print("SEQUENCE LENGTH STATISTICS")
    print("=" * 80)
    print()
    print(f"Min length: {results_df['length'].min()}")
    print(f"Max length: {results_df['length'].max()}")
    print(f"Mean length: {results_df['length'].mean():.1f}")
    print(f"Median length: {results_df['length'].median():.1f}")
    print()

    # Statistical comparison
    if VIENNA_AVAILABLE and len(results_df) > 1:
        print("=" * 80)
        print("STATISTICAL COMPARISON")
        print("=" * 80)
        print()

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
    print(f"Benchmark complete! Check {output_file} for detailed results.")
    print("=" * 80)

if __name__ == "__main__":
    main()
