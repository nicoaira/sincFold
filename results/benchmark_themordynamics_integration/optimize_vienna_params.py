import pandas as pd
import numpy as np
import torch as tr
from tqdm import tqdm
import warnings
import tempfile
import os
import json
from datetime import datetime
warnings.filterwarnings('ignore')

# Import sincfold modules
from sincfold.model import sincfold
from sincfold.metrics import f1_strict
from sincfold.utils import dot2bp, VIENNA_AVAILABLE
from sincfold.dataset import SeqDataset, pad_batch
from torch.utils.data import DataLoader

if not VIENNA_AVAILABLE:
    raise ImportError("ViennaRNA is required for this optimization")

import RNA

def fold_sincfold_vienna(model, sequence, seq_id, config, vienna_weight=1.0, vienna_temp=37.0, vienna_linear=False):
    """Fold RNA sequence using sincFold + ViennaRNA integration."""
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

def test_configuration(model, test_data, config, vienna_weight, vienna_temp, vienna_linear):
    """Test a specific Vienna parameter configuration."""
    f1_scores = []

    for _, row in test_data.iterrows():
        seq_id = row['id']
        sequence = row['sequence']
        ref_structure = row['secondary_structure']

        ref_bp = dot2bp(ref_structure)
        if ref_bp is False:
            continue

        try:
            pred_bp = fold_sincfold_vienna(
                model, sequence, seq_id, config,
                vienna_weight=vienna_weight,
                vienna_temp=vienna_temp,
                vienna_linear=vienna_linear
            )
            _, _, f1 = f1_strict(ref_bp, pred_bp)
            f1_scores.append(f1)
        except Exception as e:
            continue

    return np.mean(f1_scores) if f1_scores else 0.0

def main():
    print("=" * 80)
    print("VIENNA PARAMETER OPTIMIZATION")
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
    param_grid = {
        'vienna_weight': [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0],
        'vienna_temp': [37.0],  # Keep standard for now
        'vienna_linear': [True, False]
    }

    print("Parameter grid:")
    for key, values in param_grid.items():
        print(f"  {key}: {values}")
    print()

    # Calculate total combinations
    total_configs = (len(param_grid['vienna_weight']) *
                     len(param_grid['vienna_temp']) *
                     len(param_grid['vienna_linear']))
    print(f"Total configurations to test: {total_configs}")
    print()

    # Test all configurations
    results = []

    with tqdm(total=total_configs, desc="Testing configurations") as pbar:
        for weight in param_grid['vienna_weight']:
            for temp in param_grid['vienna_temp']:
                for linear in param_grid['vienna_linear']:
                    scaling = "linear" if linear else "log"
                    config_name = f"w={weight:.1f}, T={temp:.0f}, {scaling}"

                    mean_f1 = test_configuration(
                        model, sampled_df, config,
                        vienna_weight=weight,
                        vienna_temp=temp,
                        vienna_linear=linear
                    )

                    results.append({
                        'vienna_weight': weight,
                        'vienna_temp': temp,
                        'vienna_linear': linear,
                        'scaling': scaling,
                        'mean_f1': mean_f1,
                        'config_name': config_name
                    })

                    pbar.set_postfix({'config': config_name, 'F1': f'{mean_f1:.4f}'})
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
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f'vienna_optimization_results_{timestamp}.csv'
    results_df.to_csv(output_file, index=False)
    print(f"Detailed results saved to: {output_file}")
    print()

    # Print top 10 configurations
    print("TOP 10 CONFIGURATIONS:")
    print("-" * 80)
    print(f"{'Rank':<6} {'Weight':<10} {'Temp':<8} {'Scaling':<10} {'Mean F1':<10}")
    print("-" * 80)

    for idx, row in results_df.head(10).iterrows():
        print(f"{results_df.index.get_loc(idx)+1:<6} {row['vienna_weight']:<10.1f} "
              f"{row['vienna_temp']:<8.0f} {row['scaling']:<10} {row['mean_f1']:<10.4f}")

    print()

    # Print baseline comparison
    baseline_f1 = 0.7509  # From previous benchmark
    sincfold_f1 = 0.7793  # sincFold alone

    best_config = results_df.iloc[0]
    print("=" * 80)
    print("COMPARISON WITH BASELINES:")
    print("-" * 80)
    print(f"sincFold alone:                F1 = {sincfold_f1:.4f}")
    print(f"sincFold+Vienna (default):     F1 = {baseline_f1:.4f} (weight=1.0, log)")
    print(f"sincFold+Vienna (BEST):        F1 = {best_config['mean_f1']:.4f} "
          f"(weight={best_config['vienna_weight']}, {best_config['scaling']})")
    print()

    improvement = best_config['mean_f1'] - baseline_f1
    if improvement > 0:
        print(f"✓ IMPROVEMENT: +{improvement:.4f} F1 score ({improvement/baseline_f1*100:.2f}%)")
    else:
        print(f"✗ No improvement over default parameters")
    print()

    # Analyze scaling method
    print("=" * 80)
    print("ANALYSIS BY SCALING METHOD:")
    print("-" * 80)
    for scaling in ['linear', 'log']:
        subset = results_df[results_df['scaling'] == scaling]
        print(f"{scaling.upper():>8}: Mean F1 = {subset['mean_f1'].mean():.4f}, "
              f"Best F1 = {subset['mean_f1'].max():.4f}, "
              f"Worst F1 = {subset['mean_f1'].min():.4f}")
    print()

    # Analyze weight sensitivity
    print("=" * 80)
    print("ANALYSIS BY WEIGHT (averaging over scaling methods):")
    print("-" * 80)
    weight_analysis = results_df.groupby('vienna_weight')['mean_f1'].agg(['mean', 'std', 'min', 'max'])
    print(weight_analysis.to_string())
    print()

    # Recommendations
    print("=" * 80)
    print("RECOMMENDATIONS:")
    print("-" * 80)
    print(f"Best configuration: weight={best_config['vienna_weight']}, "
          f"temp={best_config['vienna_temp']:.0f}, scaling={best_config['scaling']}")

    if best_config['mean_f1'] > sincfold_f1:
        print(f"✓ This configuration OUTPERFORMS sincFold alone!")
    elif best_config['mean_f1'] > baseline_f1:
        print(f"✓ This configuration improves over default Vienna parameters")
    else:
        print(f"✗ Vienna integration does not improve performance")
        print(f"  Recommendation: Use sincFold alone (F1={sincfold_f1:.4f})")
    print("=" * 80)

if __name__ == "__main__":
    main()
