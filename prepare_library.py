"""Prepare local structure inputs or import an independently obtained frozen snapshot."""
import argparse, shutil, sys
from pathlib import Path
sys.dont_write_bytecode = True
import pandas as pd
from src.util import ROOT, read_json, sha, verify_inputs
from src.library_access import validate, read_library

def main():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--from-results', action='store_true')
    g.add_argument('--input', type=Path, help='Locally obtained original frozen CSV snapshot')
    a = p.parse_args()
    verify_inputs()
    if a.input:
        source = a.input.resolve()
        if not source.is_file(): p.error('Input file does not exist')
        if sha(source) != read_json(ROOT/'data/library/library_manifest.json')['original_snapshot_sha256']:
            p.error('Original snapshot SHA-256 differs from the frozen selection')
        validate(pd.read_csv(source, float_precision='round_trip'))
        target = ROOT/'data/library/public_purchasable_100k.csv'
        if target.exists(): p.error('Target already exists; existing input is preserved')
        shutil.copyfile(source, target)
    else:
        source = ROOT/'results/official_100k_scores.parquet'
        frame = pd.read_parquet(source, columns=['library_id','SMILES']).rename(columns={'SMILES':'isomeric_smiles'})
        validate(frame)
        target = ROOT/'data/library/reproduction_input.csv'
        if target.exists(): p.error('Target already exists; existing input is preserved')
        frame.to_csv(target, index=False, encoding='utf-8-sig')
    print('VERIFIED', target.relative_to(ROOT).as_posix())
if __name__ == '__main__': main()
