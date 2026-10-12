"""Refresh public package manifests; scientific model locks remain separate."""
from pathlib import Path
import hashlib, json
ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''): h.update(chunk)
    return h.hexdigest()

def main():
    excluded = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache'}
    files = sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(v in p.relative_to(ROOT).parts for v in excluded))
    manifest = {}
    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        if rel in ['IMMUTABLE_SHA256.json','checksums.sha256','data/library/public_purchasable_100k.csv','data/library/reproduction_input.csv']: continue
        if rel.startswith(('results/','logs/validation/','logs/reproduction/')): continue
        manifest[rel] = digest(path)
    (ROOT/'IMMUTABLE_SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n','utf-8')
    files = sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(v in p.relative_to(ROOT).parts for v in excluded))
    (ROOT/'checksums.sha256').write_text('\n'.join(f'{digest(p)}  {p.relative_to(ROOT).as_posix()}' for p in files if p.name!='checksums.sha256')+'\n','utf-8')
    print('UPDATED', len(manifest), 'immutable inputs')
if __name__ == '__main__': main()
