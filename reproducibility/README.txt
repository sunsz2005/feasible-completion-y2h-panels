Primary comparison sealed-output reproduction archive

The original scripts, inputs, protocols, locks, panels, and results are in
primary_archive/. Retained file contents and the historical manifest hashes are
unchanged. Two historical lock files are omitted; see
results/primary/replay_metadata.json from the repository root. The original manifest and locks retain the old tmp/ path labels.

To replay the sealed primary scores, use Python 3.12.14 with numpy==2.3.5
and openpyxl==3.1.5. From the repository root run:
  python tools/restore_yeast_inputs.py
  python tools/replay_primary.py

The replay helper creates an ignored local compatibility copy at the old
path, then replays the original scoring routine using retained dependency hashes and portable metadata. The audit verifies 360 panels,
input/split identities, paired statistics, timing arithmetic, and six
mock-clock boundary tests. It does not regenerate wall-clock searches or
overwrite original records. The audit rebases original absolute lock keys
for verification.
