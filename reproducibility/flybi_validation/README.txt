FlyBi A/B reproducibility extract

This directory contains the FlyBi interaction-table input, derived candidate identifiers, 20 paired train/evaluation splits, 20 calibration records, 40 A/B panel records, source modules for the common solver, an A/B-only score summary, and a direct score replay.

Environment: Python 3.11.7, NumPy 1.26.4 and OpenPyXL 3.0.10. Original generation used Windows 10 build 19045 on an Intel Core i7-6700 at 3.40 GHz (four physical cores, eight logical processors). In a Python 3.11.7 environment, install dependencies from this directory with:
python -m pip install -r requirements.txt

The fresh-rerun entry point sets OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS to 1 before importing NumPy or solver modules, overriding inherited values to match the original generation limits. A and B execute serially.

Verify the archived results: python replay_ab.py
Run a fresh timed A/B pair on one archived split: python rerun_ab.py --seed 20261113

The replay verifies the original pre-evaluation hashes copied into pre_evaluation_seal_ab.json and independently recomputes capacities, exclusivity, coverage, training and held-out cross-role counts. It does not rerun timed optimization. The fresh rerun writes to a separate reruns/ directory and does not alter the archived panels. Its wall-clock outputs need not match archived outputs exactly. The src modules contain the A/B constructor and common solver; the primary analysis and protocol are also documented in the cited GitHub repository.

The pre-evaluation seal subset was extracted after scoring; it records original hashes and creation time without claiming that this subset file itself existed before scoring. Full raw provenance remains in the local research record. The shared clone-universe role eligibility is a computational assumption.

The original clone-table TSV is excluded from this Git staging tree. Its source URL and SHA-256 are in data_audit.json; python tools/restore_flybi_inputs.py retrieves and verifies it from the repository root if needed for the input audit. Replay and fresh optimization use the archived names and partitions and do not require that TSV.
