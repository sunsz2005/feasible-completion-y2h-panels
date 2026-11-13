# Feasible Completion Y2H Panels

Code and recorded results for **Feasible completion-set construction for coverage-constrained yeast two-hybrid panel optimization** by Feng Sun and Haiye Liang.

The main comparison is completion-set initialization (B) versus coverage-first initialization (A), followed by the same feasible joint solver under the same total time allowance. The effect includes each initial panel's use as a restart anchor. The repository also contains the manuscript's matched controls, role-order ablations, coverage sensitivity and FlyBi validation. All reported A/B runs are retained, including negative differences and timing fallbacks.

## Contents

| Directory | Manuscript content |
| --- | --- |
| `src/` | Seven original algorithm modules; A/B sensitivity workflow and shared input helpers. |
| `results/primary/` | Primary score and sealed-panel indexes: 360 paths and 180 pairs. |
| `results/coverage_sensitivity/` | All 360 A/B paths, paired summaries, coverage-neighborhood measurements and provenance. The 120 reused paths at tau=0.98 are in `reference_batch/`. |
| `data/` | Primary CSV exports and reusable input snapshots; source attribution in `DATA_SOURCES.md`. |
| `reproducibility/primary_archive/` | Original protocols, code, locks, inputs and run records supporting the primary experiment and supplementary controls. Retained file bytes and original hash labels are preserved; two historical locks are omitted as described below. |
| `reproducibility/flybi_validation/` | A/B-only protocol extract, 20 splits, 40 panel records, calibration records, original hash subset, score replay and fresh-rerun entry point. |
| `tools/` | Input restoration, integrity verification, result replay, summary and neighborhood-probe commands. |
| `tests/` | Original score-audit and timing-boundary checks. |

Original primary audit reports are preserved under `reproducibility/primary_audits/`.

`reproducibility/output/` and `reruns/` are generated locally and ignored by Git. They are not part of the distributed results. The primary replay helper creates an ignored compatibility directory for paths embedded in the original locks; the published archive itself is named `primary_archive/`.

## Check the download

From the repository root, run:

```bash
python tools/verify_release.py
```

This standard-library check verifies the current distribution manifest, retained primary hashes, original sensitivity panel hashes, FlyBi hash subset, and primary CSV exports. It does not run optimization. `distribution_manifest.json` identifies this packaging version; it is not a pre-experiment seal.

## Primary results and supplementary controls

Use **Python 3.12.14** with the original dependencies:

```bash
python -m pip install -r requirements.txt
python tools/restore_yeast_inputs.py
python tools/replay_primary.py
```

The two excluded yeast snapshots are retrieved from the original source and checked against their recorded hashes. The audit verifies 360 panels, input/split identities, 180 paired statistics, nine summary groups, timing arithmetic and six controlled-clock boundary tests. The original computation used an Intel Core i7-6700 at 3.40 GHz, Windows 10 Home China 22H2, NumPy 2.3.5 and OpenPyXL 3.1.5. No explicit native-library thread limit was set for these primary runs.

Two historical `lock.json` files (the primary initialization-time and human external stages) are omitted from this distribution. Their original bytes remain in the separate local research archive. `results/primary/replay_metadata.json` preserves their original digests, the primary configuration and software versions, and retained dependency hashes using relative paths. This is post-scoring metadata, not a replacement pre-run lock. Use `tools/replay_primary.py` for the current package; direct execution of the archived historical audit requires the original locks. The current verifier reports the two omissions explicitly. For a portability check in a different software environment, add `--allow-version-difference`; its report states whether versions match.

The primary archive also retains the development, coverage-timing, role-order, fixed-time and EP-Prox-inspired controls described in Supplementary Sections S1-S7. Their locks form part of the primary audit's dependency chain. They are retained even where a control did not favor the proposed constructor.

## Coverage sensitivity

Use a separate **Python 3.11.7** environment:

```bash
python -m pip install -r requirements_followup.txt
python tools/restore_yeast_inputs.py
python src/coverage_sensitivity_benchmark.py score
python tools/summarize_coverage_sensitivity.py --scores reproducibility/output/coverage_sensitivity/scores.json
python tools/probe_output_neighborhood.py
```

The wrappers set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` to 1. Original generation used NumPy 1.26.4 and OpenPyXL 3.0.10 on Windows 10 build 19045 with the same Intel Core i7-6700 (four physical cores, eight logical processors). Replay outputs go to `reproducibility/output/coverage_sensitivity/`; archived scores and panels are not overwritten. The probe is post hoc and descriptive, and does not establish a causal effect of neighborhood size.

The 0.90/0.95 paths and the reused 0.98 paths came from two separately sealed batches. `provenance.json` preserves all 360 original panel hashes, recorded environment, original seal/lock digests and retained input/solver hashes. It and `protocol_extract.md` were assembled after scoring. They do not claim to be the original pre-run documents; the complete historical records remain in a separate local archive. The extracted shared input helper removes dependencies on experiments outside this manuscript without changing the original algorithm modules or recorded results.

For a new timed **0.90/0.95** batch, choose a new directory:

```bash
python src/coverage_sensitivity_benchmark.py lock --output reruns/sensitivity
python src/coverage_sensitivity_benchmark.py generate --output reruns/sensitivity
python src/coverage_sensitivity_benchmark.py score --output reruns/sensitivity
python tools/summarize_coverage_sensitivity.py --scores reruns/sensitivity/replay/scores.json
```

As in the original design, this command reuses the archived 0.98 A/B batch. It does not present all three thresholds as newly generated. Time allowances use the archived per-instance calibration; wall-clock trajectories can change with hardware and load. New-run locks record the actual rerun environment. No evaluation edges are used to select optimization states.

## FlyBi validation

Use the same Python 3.11.7 follow-up environment:

```bash
python reproducibility/flybi_validation/replay_ab.py
python reproducibility/flybi_validation/rerun_ab.py --seed 20261113
```

The first command verifies and rescores the 40 archived panels without timed optimization. The second generates a fresh timed A/B pair with the three native-library thread limits set to 1 before NumPy is imported. It writes to that directory's ignored `reruns/` folder. Timing outputs need not match the archive byte for byte.

The protocol extract and seal subset were packaged after scoring; their hash entries refer to the original pre-evaluation seal. Shared role eligibility is a computational assumption. The raw clone-list TSV is excluded from this Git tree; its URL and exact hash are in `data_audit.json`. If needed for the input audit, retrieve it with `python tools/restore_flybi_inputs.py`. The fixed candidate identifiers and splits suffice for replay.

## Citation and licensing

See `CITATION.cff` for the manuscript authors and title. Original project code is under the MIT license. Third-party data retain their source terms; the MIT license does not apply to those datasets. Cite the original publications and data providers in `DATA_SOURCES.md` when reusing inputs. The excluded YDIP/YBIOGRID snapshots must not be included in a public upload without redistribution permission.
