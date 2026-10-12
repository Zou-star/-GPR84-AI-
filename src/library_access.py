"""Read frozen library identities from local data or published screening results."""
from pathlib import Path
import hashlib
import pandas as pd
from src.util import ROOT, read_json, sha

def validate(frame):
    manifest = read_json(ROOT / 'data/library/library_manifest.json')
    col = next((c for c in ['isomeric_smiles', 'SMILES', 'smiles', 'canonical_smiles'] if c in frame), None)
    if col is None or 'library_id' not in frame:
        raise ValueError('Library input requires library_id and a SMILES column')
    if len(frame) != manifest['rows'] or frame.library_id.isna().any() or frame[col].isna().any():
        raise ValueError('Library row count or required values differ from the frozen structure selection')
    ids = hashlib.sha256('\n'.join(frame.library_id.astype(str)).encode()).hexdigest()
    smiles = hashlib.sha256('\n'.join(frame[col].astype(str)).encode()).hexdigest()
    if ids != manifest['library_id_order_sha256'] or smiles != manifest['input_smiles_order_sha256']:
        raise ValueError('Library identity/order differs from the frozen structure selection')
    return col

def read_library():
    original = ROOT / 'data/library/public_purchasable_100k.csv'
    prepared = ROOT / 'data/library/reproduction_input.csv'
    if original.is_file():
        path, mode = original, 'locally_imported_original_snapshot'
        if sha(path) != read_json(ROOT / 'data/library/library_manifest.json')['original_snapshot_sha256']:
            raise ValueError('Imported original snapshot hash differs from the frozen snapshot')
        data = pd.read_csv(path, float_precision='round_trip')
    elif prepared.is_file():
        path, mode = prepared, 'locally_prepared_screening_result_structures'
        data = pd.read_csv(path, float_precision='round_trip')
    else:
        path, mode = ROOT / 'results/official_100k_scores.parquet', 'published_screening_result_structures'
        if not path.is_file():
            raise FileNotFoundError('Restore results/official_100k_scores.parquet or run prepare_library.py --input <local snapshot.csv>')
        data = pd.read_parquet(path, columns=['library_id', 'SMILES'])
        data = data.rename(columns={'SMILES': 'isomeric_smiles'})
    validate(data)
    return data, {'input_mode': mode, 'input_file': path.relative_to(ROOT).as_posix(),
                  'input_file_sha256': sha(path), 'input_columns_used': ['library_id', 'SMILES'],
                  'saved_scores_used_for_inference': False}
