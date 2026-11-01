Primary comparison sealed-output reproduction archive

Requires Python 3.12.14, numpy==2.3.5, openpyxl==3.1.5.
Install those exact dependencies in a suitable Python environment.
From the extracted archive root run:
  python tmp/audit_primary_comparison.py

The audit verifies frozen hashes, all 360 panels, input/split identities,
paired statistics and timing arithmetic. It replays sealed-panel scores and
independently counts evaluation edges. Six mock-clock tests check timing rules.
It does not regenerate wall-clock searches or overwrite original records.
Only output/research/primary_comparison_audit_2026-10-25.json is refreshed.

Original experiment scripts and locks retain their original absolute keys.
The audit rebases those keys for verification; frozen files are unchanged.
Use manifest_sha256.json to verify all packaged originals. Hardware-dependent
search trajectories are distinct from exact sealed-output score replay.
