"""Fetch excluded yeast inputs from their original public source and verify them."""

import hashlib
from pathlib import Path
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parent.parent
SOURCE = "https://raw.githubusercontent.com/paopaopig/The-construction-of-the-CM-PIN./main/datasets/yeast"
INPUTS = {
    "YBIOGRID": "64897f90c458dd2b7ceadf922d3313c688fef91fc1a8ae588c8b8a3338be4370",
    "YDIP": "ded6cc2c13a63c08c91bead1e1d213dd689b3faf3f3bf81832c8624812700219",
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def restore(name, expected):
    archive = ROOT / "reproducibility/primary_archive/ppi_data_audit" / name / "interactions.txt"
    browsable = ROOT / "data/inputs" / name / "interactions.txt"
    for path in (archive, browsable):
        if path.exists() and digest(path.read_bytes()) != expected:
            raise ValueError(f"Existing file has an unexpected SHA-256: {path}")

    if archive.exists():
        data = archive.read_bytes()
    elif browsable.exists():
        data = browsable.read_bytes()
    else:
        url = f"{SOURCE}/{name}/interactions.txt"
        with urlopen(url, timeout=60) as response:
            downloaded = response.read()
        normalized = downloaded.replace(b"\r\n", b"\n")
        candidates = (downloaded, normalized, normalized.replace(b"\n", b"\r\n"))
        data = next((item for item in candidates if digest(item) == expected), None)
        if data is None:
            raise ValueError(f"Upstream {name} differs from the archived experiment; no file was written")

    for path in (archive, browsable):
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            path.write_bytes(data)
    print(f"OK: {name} SHA-256 {expected}")


if __name__ == "__main__":
    for source_name, source_hash in INPUTS.items():
        restore(source_name, source_hash)
