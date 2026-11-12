Primary comparison sealed-output reproduction archive

The original scripts, inputs, protocols, locks, panels, and results are in
primary_archive/. File contents and the historical manifest hashes are
unchanged. The original manifest and locks retain the old tmp/ path labels.

To replay the sealed primary scores, use Python 3.12.14 with numpy==2.3.5
and openpyxl==3.1.5. From the repository root run:
  python tools/restore_yeast_inputs.py
  python tools/replay_primary.py

The replay helper creates an ignored local compatibility copy at the old
path, then runs the unchanged archived audit. The audit verifies 360 panels,
input/split identities, paired statistics, timing arithmetic, and six
mock-clock boundary tests. It does not regenerate wall-clock searches or
overwrite original records. The audit rebases original absolute lock keys
for verification.
