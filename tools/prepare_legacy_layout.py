"""Materialize an ignored compatibility copy for the original locked scripts."""

import hashlib
from pathlib import Path
from shutil import copy2


ROOT = Path(__file__).resolve().parents[1]
REPRODUCIBILITY = ROOT / "reproducibility"
SOURCE = REPRODUCIBILITY / "primary_archive"
COMPATIBILITY = REPRODUCIBILITY / "tmp"


def digest(path):
    return hashlib.sha256(path.read_bytes()).digest()


def prepare():
    if not SOURCE.is_dir():
        raise FileNotFoundError(SOURCE)
    copied = 0
    verified = 0
    for original in sorted(SOURCE.rglob("*")):
        if not original.is_file():
            continue
        if original.is_symlink():
            raise ValueError(f"Unexpected symlink in archive: {original}")
        destination = COMPATIBILITY / original.relative_to(SOURCE)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if not destination.is_file() or digest(destination) != digest(original):
                raise ValueError(f"Compatibility file differs: {destination}")
            verified += 1
        else:
            copy2(original, destination)
            copied += 1
    print(f"Compatibility copy ready: {copied} copied, {verified} verified")


if __name__ == "__main__":
    prepare()
