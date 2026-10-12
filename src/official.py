"""Portable submission adapter for unchanged official weights and chemistry.

This submission lock is separate from the original release lock, which remains
in the local archive. Scores and chemistry follow the original prediction API.
"""
import importlib.metadata
import platform
import sys
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from src.util import ROOT, OFF, read_json, sha

sys.dont_write_bytecode = True
VERSION = 'GPR84-ML-v3-official-1'
WEIGHT_SHA256 = 'df5be03e49030ef885b24d1ed2854c2ae9fbc45858aadd298dd92af6e6328b1f'
THRESHOLD = 0.046076174825429916

def verify_lock(base=OFF):
    """Check immutable inference files and exact trained-runtime versions."""
    base=Path(base).resolve()
    lock=read_json(ROOT/'models/SUBMISSION_MODEL_LOCK.json')
    if lock['version']!=VERSION or lock['threshold']!=THRESHOLD:
        raise RuntimeError('Submission model version/threshold mismatch')
    if lock['source_model_sha256']!=WEIGHT_SHA256 or sha(base/'model.joblib')!=WEIGHT_SHA256:
        raise RuntimeError('Official model hash mismatch: model.joblib')
    for relative, expected in lock['critical_sha256'].items():
        p=(base/relative).resolve()
        if not p.is_relative_to(base) or not p.is_file() or sha(p)!=expected:
            raise RuntimeError('Submission inference hash mismatch: '+relative)
    for relative, expected in lock['training_sources_sha256'].items():
        p=(ROOT/'data/training'/relative).resolve()
        if not p.is_relative_to(ROOT/'data/training') or not p.is_file() or sha(p)!=expected:
            raise RuntimeError('Frozen training source hash mismatch: '+relative)
    env=lock['environment']
    if platform.python_version_tuple()[:2]!=tuple(env['python'].split('.')[:2]):
        raise RuntimeError('Python minor version differs from locked runtime')
    for package, version in env['packages'].items():
        if importlib.metadata.version(package)!=version:
            raise RuntimeError('Locked dependency version mismatch: '+package)
    return lock

def load_model():
    lock=verify_lock()
    bundle=joblib.load(OFF/'model.joblib')
    if bundle['threshold']!=lock['threshold'] or bundle['version']!=lock['source_model_version']:
        raise RuntimeError('Bundle metadata differs from submission lock')
    return lock,bundle

def classify(score,threshold):
    return int(float(score)>=float(threshold))

def predict(smiles):
    lock,bundle=load_model()
    sys.path.insert(0,str(OFF))
    from runtime.common import standardize
    from runtime.features import features,tanimoto
    rows,valid,positions=[],[],[]
    for position,value in enumerate(smiles):
        text=value if isinstance(value,str) else ''
        row=dict(row_id=position,SMILES=text,model_version=VERSION,
                 official_release=True,independent_acceptance=False,
                 status='invalid_SMILES',standardized_smiles=None,
                 full_inchikey=None,connectivity_key=None,ranking_score=None,
                 fixed_threshold=lock['threshold'],fixed_threshold_class=None,
                 max_training_similarity=None,descriptors_in_training_range=None,
                 applicability_domain=None)
        try:structure=standardize(text) if text.strip() else None
        except Exception:structure=None
        if structure is not None and structure['heavy_atom_count']>0:
            row.update(status='OK',standardized_smiles=structure['isomeric_smiles'],
                       full_inchikey=structure['full_inchikey'],connectivity_key=structure['connectivity_key'])
            valid.append(structure);positions.append(position)
        rows.append(row)
    if valid:
        arrays=features(pd.DataFrame(valid))
        scores=bundle['model'].predict_proba(arrays['F4'])[:,1]
        similarities=tanimoto(arrays['ECFP4'],bundle['reference_bits']).max(axis=1)
        desc=arrays['F1'][:,-10:]
        within=np.all((desc>=bundle['desc_min'])&(desc<=bundle['desc_max']),axis=1)
        for i,position in enumerate(positions):
            p,similarity=float(scores[i]),float(similarities[i])
            rows[position].update(ranking_score=p,fixed_threshold_class=classify(p,lock['threshold']),
                max_training_similarity=similarity,descriptors_in_training_range=bool(within[i]),
                applicability_domain=bool(similarity>=.4 and within[i]))
    columns=['row_id','SMILES','model_version','official_release','independent_acceptance','status',
             'standardized_smiles','full_inchikey','connectivity_key','ranking_score','fixed_threshold',
             'fixed_threshold_class','max_training_similarity','descriptors_in_training_range','applicability_domain']
    return pd.DataFrame(rows,columns=columns)
