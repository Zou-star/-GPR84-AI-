"""Portable paths, immutable-input verification and auditable JSON output."""
from pathlib import Path
import hashlib, json, platform, time
ROOT = Path(__file__).resolve().parents[1]
OFF = ROOT / 'models/official_model'
def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''): h.update(chunk)
    return h.hexdigest()
def read_json(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def dump(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def verify_inputs():
    manifest=read_json(ROOT/'IMMUTABLE_SHA256.json')
    for name, digest in manifest.items():
        p=(ROOT/name).resolve()
        if not p.is_relative_to(ROOT) or not p.is_file() or sha(p)!=digest:
            raise RuntimeError('Immutable submission file hash mismatch: '+name)
    return len(manifest)
def _check_output(path):
    p=Path(path).resolve()
    if p.is_relative_to(ROOT):
        allowed=[ROOT/'results',ROOT/'logs/validation']
        if not any(p.is_relative_to(s) for s in allowed):
            raise ValueError('Outputs inside project must use results/ or logs/validation/')
    return p

def output_dir(path):
    p=_check_output(path)
    if p.exists() and not p.is_dir():raise ValueError('Output directory is a file')
    p.mkdir(parents=True,exist_ok=True);return p

def output_file(path):
    p=_check_output(path)
    if p.exists() and p.is_dir():raise ValueError('Output file is a directory')
    output_dir(p.parent);return p

def scoring_signature(frame):
    """Semantic hashes of all library rows; independent of Parquet metadata."""
    import numpy as np
    columns=['ranking_score','max_training_similarity','applicability_domain',
             'fixed_threshold_class','eligible_for_ranking']
    return {'N':len(frame),
        'library_id_order_sha256':hashlib.sha256('\n'.join(frame.library_id).encode()).hexdigest(),
        'columns':{c:{'dtype':str(frame[c].to_numpy().dtype),
            'sha256':hashlib.sha256(np.ascontiguousarray(frame[c].to_numpy()).tobytes()).hexdigest()}
            for c in columns}}
def write_csv(frame,path):
    frame.to_csv(path,index=False,encoding='utf-8-sig',float_format='%.17g')
