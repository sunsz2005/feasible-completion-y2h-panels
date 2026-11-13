"""Retrieve the original FlyBi clone table only when needed for the input audit."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1] / 'reproducibility/flybi_validation'


def main():
    audit = json.loads((ROOT / 'data_audit.json').read_text(encoding='utf-8'))
    rel = 'raw/clone_list_20230307.tsv'
    path = ROOT / rel
    expected = audit['raw_sha256'][rel.replace('/', '\\')]
    data = path.read_bytes() if path.exists() else urlopen(audit['sources'][rel], timeout=60).read()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError('Clone table differs from the archived input; no file was written')
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_bytes(data)
    print('Verified FlyBi clone-table SHA-256:', expected)


if __name__ == '__main__':
    main()
