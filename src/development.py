"""Recompute saved development predictions and bootstrap intervals without fitting."""
import argparse
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, matthews_corrcoef, balanced_accuracy_score, confusion_matrix

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/training'
REPLICATES = 2000
SEED = 42
METRICS = ['ROC_AUC', 'MCC', 'Balanced_Accuracy']


def estimates(y, score, threshold):
    pred = score >= threshold
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {'N': len(y), 'positive': int(y.sum()), 'negative': int((y == 0).sum()),
            'ROC_AUC': float(roc_auc_score(y, score)), 'MCC': float(matthews_corrcoef(y, pred)),
            'Balanced_Accuracy': float(balanced_accuracy_score(y, pred)),
            'TP': int(tp), 'TN': int(tn), 'FP': int(fp), 'FN': int(fn), 'threshold': threshold}


def compute(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # Read the original search export, not the reserialized candidate copy.
    oof = pd.read_csv(BASE / 'development_publication_OOF.csv', float_precision='round_trip')
    members = pd.read_csv(BASE / 'training_members.csv')
    structures = pd.read_csv(BASE / 'training_structures.csv', float_precision='round_trip')
    groups = pd.read_csv(BASE / 'publication_groups.csv')
    result = json.loads((ROOT / 'logs/history/XGBoost_F4_result.json').read_text(encoding='utf-8'))
    threshold = result['threshold']
    for frame in [members, structures, groups]:
        if not np.array_equal(oof.compound_id, frame.compound_id) or not np.array_equal(oof.label, frame.label):
            raise RuntimeError('Development OOF membership/order/label mismatch')
    if len(oof) != 400 or oof.label.sum() != 326 or not np.isfinite(oof.score).all():
        raise RuntimeError('Unexpected development class support or invalid score')
    y, score = oof.label.to_numpy(int), oof.score.to_numpy(float)
    point = estimates(y, score, threshold)
    for metric in METRICS + ['TP', 'TN', 'FP', 'FN']:
        if abs(point[metric] - result['publication_OOF'][metric]) > 1e-12:
            raise RuntimeError('Original OOF metric mismatch: ' + metric)
    negative, positive = np.flatnonzero(y == 0), np.flatnonzero(y == 1)
    unit_names = sorted(groups.group.unique())
    units = [np.flatnonzero(groups.group.to_numpy() == name) for name in unit_names]
    draws = []
    for kind in ['stratified_compound', 'publication_connected_group']:
        rng = np.random.default_rng(SEED)
        for replicate in range(REPLICATES):
            if kind == 'stratified_compound':
                indices = np.r_[rng.choice(negative, len(negative), replace=True),
                                rng.choice(positive, len(positive), replace=True)]
            else:
                sampled = rng.integers(0, len(units), len(units))
                indices = np.concatenate([units[index] for index in sampled])
            label, probability = y[indices], score[indices]
            valid = len(np.unique(label)) == 2
            row = {'method': kind, 'replicate': replicate, 'valid_two_classes': valid,
                   'N': len(indices), 'positive': int(label.sum()), 'negative': int((label == 0).sum())}
            value = estimates(label, probability, threshold) if valid else {m: np.nan for m in METRICS}
            row.update({m: value[m] for m in METRICS})
            draws.append(row)
    boot = pd.DataFrame(draws)
    summaries = {}
    for method, frame in boot.groupby('method', sort=False):
        valid = frame[frame.valid_two_classes]
        if len(valid) < 1000:
            raise RuntimeError('Too few valid two-class publication-group bootstrap draws')
        summaries[method] = {'seed': SEED, 'attempted': REPLICATES, 'valid': len(valid),
                             'excluded_single_class_draws': REPLICATES - len(valid),
                             'percentile_95_CI': {m: [float(v) for v in np.percentile(valid[m], [2.5, 97.5])] for m in METRICS}}
    legacy = pd.read_csv(ROOT / 'logs/history/legacy_candidate_OOF.csv', float_precision='round_trip')
    if not np.array_equal(legacy.compound_id, oof.compound_id) or not np.array_equal(legacy.label, oof.label):
        raise RuntimeError('Legacy OOF membership mismatch')
    changed = np.flatnonzero((legacy.score.to_numpy() >= threshold) != (score >= threshold))
    precision = {'policy': 'Original CSV bytes + pandas float_precision=round_trip + float64 comparison >= threshold',
                 'threshold_changed': False, 'changed_classifications_vs_legacy_candidate_export': len(changed),
                 'legacy_metrics': estimates(y, legacy.score.to_numpy(), threshold),
                 'boundary_rows': [{'compound_id': oof.iloc[i].compound_id, 'original_score': float(score[i]),
                                    'legacy_export_score': float(legacy.iloc[i].score), 'threshold': threshold,
                                    'class_from_original': int(score[i] >= threshold),
                                    'class_from_legacy': int(legacy.iloc[i].score >= threshold)} for i in changed]}
    summary = {'evaluation_type': 'development selection estimate; post-selection bootstrap conditional on saved predictions',
               'independent_acceptance': False, 'model_refitted': False, 'threshold_refitted': False,
               'point_estimates': point, 'bootstrap': summaries,
               'publication_grouping': {'units': len(units), 'definition': 'Saved connected components of all associated publication IDs and connectivity keys',
                                        'source': 'publication_groups.csv; no splitting of connected publication units'},
               'limitations': 'Intervals describe resampling uncertainty of saved development OOF scores. They do not account for hyperparameter/model/threshold selection optimism and are not fresh independent validation.',
               'precision_audit': precision}
    (output / 'development_metrics.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    boot.to_csv(output / 'bootstrap_replicates.csv', index=False, encoding='utf-8-sig', float_format='%.17g')
    table = structures[['compound_id', 'isomeric_smiles', 'full_inchikey', 'connectivity_key', 'publication_group', 'label']].copy()
    table['OOF_group'] = groups.group
    table['ranking_score'] = score
    table['fixed_threshold'] = threshold
    table['fixed_threshold_class'] = (score >= threshold).astype(int)
    table['prediction_type'] = 'publication_isolated_development_OOF'
    table.to_csv(output / 'development_predictions.csv', index=False, encoding='utf-8-sig', float_format='%.17g')
    print(json.dumps({'point': point, 'bootstrap': summaries, 'groups': len(units), 'precision_boundary_rows': len(changed)}, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', default=str(BASE))
    compute(parser.parse_args().output_dir)
