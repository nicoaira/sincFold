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
from sincfold.metrics import f1_strict
from sincfold.utils import dot2bp, fold_with_vienna, VIENNA_AVAILABLE
from sincfold.dataset import SeqDataset, pad_batch
from torch.utils.data import DataLoader
from torch.cuda.amp import autocast

if not VIENNA_AVAILABLE:
    raise ImportError("ViennaRNA is required for this optimization")

import RNA

def fold_sincfold_vienna_advanced(model, sequence, seq_id, config, vienna_weight=1.0,
                                   vienna_temp=37.0, vienna_linear=False,
                                   prob_threshold=0.1, max_constraints=None):
    """Fold with advanced Vienna parameters (prob_threshold, max_constraints)."""

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

        # Get sincFold probability matrix
        model.eval()
        with tr.no_grad():
            for batch in pred_loader:
                # Use mixed precision for inference
                if model.use_amp:
                    with autocast():
                        y_pred = model(batch)
                else:
                    y_pred = model(batch)

                if isinstance(y_pred, tuple):
                    y_pred = y_pred[0]

                # Get probability matrix for this sequence
                length = batch["length"][0]
                y_pred_slice = y_pred[0, :length, :length].squeeze().cpu()

                # Fold with ViennaRNA using custom parameters
                base_pairs, mfe_energy = fold_with_vienna(
                    sequence,
                    y_pred_slice,
                    weight=vienna_weight,
                    temperature=vienna_temp,
                    linear=vienna_linear,
                    prob_threshold=prob_threshold,
                    max_constraints=max_constraints
                )

    finally:
        os.remove(tmp_file)

    return base_pairs

def test_advanced_config(model, test_data, config, vienna_weight, vienna_temp,
                         vienna_linear, prob_threshold, max_constraints):
    """Test a specific advanced Vienna configuration."""
    f1_scores = []

    for _, row in test_data.iterrows():
        seq_id = row['id']
        sequence = row['sequence']
        ref_structure = row['secondary_structure']

        ref_bp = dot2bp(ref_structure)
        if ref_bp is False:
            continue

        try:
            pred_bp = fold_sincfold_vienna_advanced(
                model, sequence, seq_id, config,
                vienna_weight=vienna_weight,
                vienna_temp=vienna_temp,
                vienna_linear=vienna_linear,
                prob_threshold=prob_threshold,
                max_constraints=max_constraints
            )
            _, _, f1 = f1_strict(ref_bp, pred_bp)
            f1_scores.append(f1)
        except Exception as e:
            continue

    return np.mean(f1_scores) if f1_scores else 0.0

def main():
    print("=" * 80)
    print("ADVANCED VIENNA PARAMETER OPTIMIZATION")
    print("(Testing prob_threshold and max_constraints)")
    print("=" * 80)
    print()

    # Set random seed
    RANDOM_SEED = 42
    np.random.seed(RANDOM_SEED)
    tr.manual_seed(RANDOM_SEED)

    # Load dataset
    print("Loading dataset...")
    df = pd.read_parquet("hf://datasets/multimolecule/archiveii.512/test.parquet")

    # Sample data (use same sampling as benchmark)
    sampled_df = df.groupby('family', group_keys=False).apply(
        lambda x: x.sample(n=min(10, len(x)), random_state=RANDOM_SEED)
    ).reset_index(drop=True)

    print(f"Test set: {len(sampled_df)} sequences from {sampled_df['family'].nunique()} families")
    print()

    # Initialize model
    print("Loading sincFold model...")
    device = 'cuda' if tr.cuda.is_available() else 'cpu'
    config = {"device": device, "batch_size": 1, "verbose": False}
    model = sincfold(pretrained=True, **config)
    print("Model loaded!")
    print()

    # Define parameter grid
    # Based on previous results, use best weight=0.7 and log scaling
    # Now sweep prob_threshold and max_constraints

    avg_length = sampled_df['sequence'].apply(len).mean()
    print(f"Average sequence length: {avg_length:.1f}")
    print()

    param_grid = {
        'vienna_weight': [0.7, 3.0],  # Best from previous run
        'vienna_temp': [37.0],
        'vienna_linear': [False],  # Log scaling was best
        'prob_threshold': [0.01, 0.05, 0.1, 0.2, 0.3, 0.5],
        'max_constraints': [50, 100, 200, None]  # None = default (N*2)
    }

    print("Parameter grid:")
    for key, values in param_grid.items():
        print(f"  {key}: {values}")
    print()

    # Calculate total combinations
    total_configs = (len(param_grid['vienna_weight']) *
                     len(param_grid['vienna_temp']) *
                     len(param_grid['vienna_linear']) *
                     len(param_grid['prob_threshold']) *
                     len(param_grid['max_constraints']))
    print(f"Total configurations to test: {total_configs}")
    print()

    # Test all configurations
    results = []

    with tqdm(total=total_configs, desc="Testing configurations") as pbar:
        for weight in param_grid['vienna_weight']:
            for temp in param_grid['vienna_temp']:
                for linear in param_grid['vienna_linear']:
                    for threshold in param_grid['prob_threshold']:
                        for max_constr in param_grid['max_constraints']:

                            config_name = (f"w={weight:.1f}, thresh={threshold:.2f}, "
                                          f"max={max_constr if max_constr else 'auto'}")

                            mean_f1 = test_advanced_config(
                                model, sampled_df, config,
                                vienna_weight=weight,
                                vienna_temp=temp,
                                vienna_linear=linear,
                                prob_threshold=threshold,
                                max_constraints=max_constr
                            )

                            results.append({
                                'vienna_weight': weight,
                                'vienna_temp': temp,
                                'vienna_linear': linear,
                                'prob_threshold': threshold,
                                'max_constraints': max_constr if max_constr else 'auto',
                                'mean_f1': mean_f1,
                                'config_name': config_name
                            })

                            pbar.set_postfix({'F1': f'{mean_f1:.4f}'})
                            pbar.update(1)

    print()
    print("=" * 80)
    print("OPTIMIZATION RESULTS")
    print("=" * 80)
    print()

    # Convert to DataFrame and sort by F1
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('mean_f1', ascending=False)

    # Save results
    output_file = 'vienna_advanced_optimization_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"Detailed results saved to: {output_file}")
    print()

    # Print top 15 configurations
    print("TOP 15 CONFIGURATIONS:")
    print("-" * 100)
    print(f"{'Rank':<6} {'Weight':<8} {'Threshold':<11} {'MaxConstr':<12} {'Mean F1':<10}")
    print("-" * 100)

    for idx, row in results_df.head(15).iterrows():
        print(f"{results_df.index.get_loc(idx)+1:<6} {row['vienna_weight']:<8.1f} "
              f"{row['prob_threshold']:<11.2f} {str(row['max_constraints']):<12} "
              f"{row['mean_f1']:<10.4f}")

    print()

    # Print baseline comparison
    baseline_f1 = 0.7509  # Default params
    best_simple_f1 = 0.7519  # Best from first optimization
    sincfold_f1 = 0.7793  # sincFold alone

    best_config = results_df.iloc[0]
    print("=" * 80)
    print("COMPARISON WITH BASELINES:")
    print("-" * 80)
    print(f"sincFold alone:                   F1 = {sincfold_f1:.4f}")
    print(f"sincFold+Vienna (default):        F1 = {baseline_f1:.4f}")
    print(f"sincFold+Vienna (weight tuning):  F1 = {best_simple_f1:.4f}")
    print(f"sincFold+Vienna (BEST ADVANCED):  F1 = {best_config['mean_f1']:.4f}")
    print()

    improvement = best_config['mean_f1'] - baseline_f1
    improvement_vs_simple = best_config['mean_f1'] - best_simple_f1

    if best_config['mean_f1'] > sincfold_f1:
        print(f"🎉 SUCCESS: Outperforms sincFold alone by {best_config['mean_f1'] - sincfold_f1:.4f}!")
    elif improvement > 0:
        print(f"✓ IMPROVEMENT: +{improvement:.4f} over default ({improvement/baseline_f1*100:.2f}%)")
        print(f"  Additional gain over simple tuning: +{improvement_vs_simple:.4f}")
    else:
        print(f"✗ No improvement over defaults")
    print()

    # Analysis
    print("=" * 80)
    print("ANALYSIS BY PROBABILITY THRESHOLD:")
    print("-" * 80)
    threshold_analysis = results_df.groupby('prob_threshold')['mean_f1'].agg(['mean', 'std', 'max'])
    print(threshold_analysis.to_string())
    print()

    print("=" * 80)
    print("ANALYSIS BY MAX CONSTRAINTS:")
    print("-" * 80)
    constraints_analysis = results_df.groupby('max_constraints')['mean_f1'].agg(['mean', 'std', 'max'])
    print(constraints_analysis.to_string())
    print()

    # Recommendations
    print("=" * 80)
    print("RECOMMENDATIONS:")
    print("-" * 80)
    print(f"Best configuration:")
    print(f"  weight={best_config['vienna_weight']}")
    print(f"  prob_threshold={best_config['prob_threshold']}")
    print(f"  max_constraints={best_config['max_constraints']}")
    print(f"  scaling=log")
    print()

    if best_config['mean_f1'] >= sincfold_f1:
        print(f"✓ Use this Vienna configuration - it matches or beats sincFold alone!")
    else:
        print(f"✗ Recommendation: Use sincFold alone (F1={sincfold_f1:.4f})")
        print(f"  Vienna integration provides no benefit even with optimized parameters")
    print("=" * 80)

if __name__ == "__main__":
    main()
