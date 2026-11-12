"""Run the unchanged primary score audit from its historical path layout."""

import subprocess
import sys

from prepare_legacy_layout import COMPATIBILITY, ROOT, prepare


for name in ("YDIP", "YBIOGRID"):
    if not (ROOT / "data" / "inputs" / name / "interactions.txt").is_file():
        raise SystemExit("Restore the excluded yeast inputs first: python tools/restore_yeast_inputs.py")

prepare()
subprocess.run(
    [sys.executable, str(COMPATIBILITY / "audit_primary_comparison.py")],
    cwd=ROOT,
    check=True,
)
