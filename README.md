# Feasible completion sets for dual-role screening panels

Code, input snapshots, sealed panels, and results for **“Feasible completion-set construction for coverage-constrained yeast two-hybrid panel optimization”** by Feng Sun and Haiye Liang. This directory accompanies a submission to *Algorithms for Molecular Biology*.

## Contents

| Path | Contents |
| --- | --- |
| `src/` | Seven original research scripts plus the later independent SA benchmark script. |
| `tests/` | Original score-replay audit and timing-boundary tests. |
| `data/inputs/` | Six included input and derived-mapping files. Two yeast network snapshots are obtained separately. |
| `data/primary_summary.csv` | Nine source-by-budget summary rows. |
| `data/primary_paired.csv` | 180 paired evaluation comparisons. |
| `data/primary_runs.csv` | All 360 individual optimization paths. |
| `results/primary/` | Original score and sealed-panel indexes. |
| `results/independent_sa/` | Follow-up simulated-annealing protocol, sealed panels, scores, and paired summaries. |
| `reproducibility/` | Original archive layout, including protocols, calibrations, all panel files, locks, and audit output; two yeast network snapshots are excluded from Git. |
| `reproducibility/manifest_sha256.json` | SHA-256 digests for all 479 original research files, including the two excluded inputs. |
| `tools/verify_release.py` | Standard-library check of the included archive, browsable copies, and CSV exports. |
| `tools/restore_yeast_inputs.py` | Fetch and hash-check the two excluded yeast network snapshots from their upstream repository. |
| `DATA_SOURCES.md` | Sources, attribution, and data-use notes. |

The files in `src/`, `tests/`, `data/inputs/`, and `results/primary/` are browsable copies of the original execution snapshot. The historical layout lives in `reproducibility/`; restore the two excluded yeast inputs before running the full score audit. Included archived files are unchanged from the audit archive. Original lock files contain historical Windows paths that the audit script resolves relative to this checkout. The verifier checks the included files against the manifest and their browsable copies. The CSV files are convenient exports of the original `scores.json`; the JSON is authoritative.

The later independent-SA benchmark is in `src/independent_sa_benchmark.py` and `results/independent_sa/`. It adds a separate experiment without changing the seven original scripts or the archived primary results. Its A and B paths are rerun beside the independent SA path in one software environment; they should not be substituted for the original primary summaries. The benchmark's protocol, code and input hashes, panel files, seal, scores, and paired summaries permit a separate audit.

The authors' original software code is offered under the [MIT License](LICENSE). Third-party input datasets are excluded from that license and retain their source terms; see [DATA_SOURCES.md](DATA_SOURCES.md).

## Verify the package

From this directory, run:

```bash
python tools/verify_release.py
```

This requires only the Python standard library and checks 477 included archived files, their browsable copies, CSV exports, panel count, and primary result counts. The two excluded yeast inputs remain listed in the manifest and are checked if present. This step does not rerun the time-limited optimization.

## Restore the two excluded yeast inputs

The YDIP and YBIOGRID `interactions.txt` snapshots came from the [CM-PIN upstream repository](https://github.com/paopaopig/The-construction-of-the-CM-PIN./tree/main/datasets/yeast). That repository does not state redistribution terms for these files, so they are not committed here. For a full score replay, run:

```bash
python tools/restore_yeast_inputs.py
python tools/verify_release.py
```

The restoration script downloads each file from that upstream repository, checks its SHA-256 against the historical manifest, and places verified copies in both `data/inputs/` and `reproducibility/tmp/ppi_data_audit/`. It accepts the upstream file's original or normalized line endings only when the historical hash matches. If the upstream files change, the script stops without substituting different data. The original YDIP and YBIOGRID SHA-256 values are `ded6cc2c13a63c08c91bead1e1d213dd689b3faf3f3bf81832c8624812700219` and `64897f90c458dd2b7ceadf922d3313c688fef91fc1a8ae588c8b8a3338be4370`, respectively.

## Replay sealed-panel scoring

Use **Python 3.12.14** with the pinned dependencies:

```bash
python -m pip install -r requirements.txt
python reproducibility/tmp/audit_primary_comparison.py
```

After the yeast inputs are restored, the audit checks input and implementation hashes, all 360 panel constraints and evaluation scores, 180 paired statistics, nine summary groups, and six clock-controlled boundary tests. It writes `reproducibility/output/research/primary_comparison_audit_2026-10-25.json`. Run it in a disposable copy if you need the original archive tree to remain byte-for-byte unchanged. The audit replays **sealed outputs**; it does not regenerate timed search trajectories.

## Independent simulated-annealing comparison

The follow-up comparison uses the same 20 yeast splits and 20 human search seeds at the archived primary `5*T0` budget. It compares complete A and B pipelines with a separate constrained simulated-annealing joint optimizer, rerunning all three serially on the same computer and software environment. Initialization time is included for each path. Evaluation interactions are parsed and scored only after all 180 panel outputs are sealed. See `results/independent_sa/protocol.md` for the fixed move, cooling, coverage, and timing rules. The two yeast input snapshots must be restored locally before regeneration; they remain excluded from Git.

With Python 3.11.7 and the versions in `requirements_independent_sa.txt` installed, the steps are:

```bash
python src/independent_sa_benchmark.py lock
python src/independent_sa_benchmark.py generate
python src/independent_sa_benchmark.py score
python tools/summarize_independent_sa.py
```

The lock command is for a fresh checkout only; it refuses to overwrite an existing lock. Generation likewise refuses to overwrite a completed seal. To replay scores without rerunning the timed paths, use the existing lock and seal and run only `score` and `summarize_independent_sa.py` after restoring the yeast inputs.

The lock records a SHA-256 digest of the evaluation file before generation to identify the input version. The generation algorithm does not parse evaluation-edge records or use their scores; parsing and scoring occur after the panel seal.

For the descriptive training-versus-evaluation comparison, install `requirements_train_eval_plot.txt` and run `python tools/plot_train_eval_deltas.py`. It reads the sealed `paired.csv` and writes `train_eval_correlation.json` and `train_eval_deltas.pdf` in the same results directory. The 20 human points reuse one graph pair, so the reported Spearman coefficients describe search-run variation rather than independent biological samples.

The separate coverage-threshold follow-up is in `src/coverage_sensitivity_benchmark.py` and `results/coverage_sensitivity/`. It fixes `q=0.8`, `k=200`, and the same `5*T0` total budget while testing `tau=0.90, 0.95, 0.98`. The 0.98 A/B paths come from the sealed independent-SA rerun. After restoring the excluded yeast snapshots and installing `requirements_independent_sa.txt`, its commands are `lock`, `generate`, and `score`, followed by `python tools/summarize_coverage_sensitivity.py`. The protocol explains the fixed-panel feasible-replacement probe and why its results are descriptive.

For a more informative post hoc check of the final panels' local neighborhoods, run `python tools/probe_output_neighborhood.py` after scoring. It writes `output_neighborhood.json`; `results/coverage_sensitivity/output_neighborhood_protocol.md` defines the sampled move and its limitations.

## Reproduction scope

The original comparison used 20 yeast split/search seeds for each of YDIP and YBIOGRID, and 20 search seeds on one historical human graph pair. Budgets were calibrated per training instance at `T0`, `5T0`, and `20T0`. Evaluation edges were counted only after outputs had been selected by training scores and sealed. The run environment was an Intel Core i7-6700 CPU at 3.40 GHz with eight logical processors, Windows 10 Home China 22H2, Python 3.12.14, NumPy 2.3.5, and openpyxl 3.1.5. Optimization paths ran sequentially. No OpenMP, OpenBLAS, or MKL thread limit was set; NumPy's effective native thread count was not logged. Wall-clock trajectories can vary on other machines.

## Before public release

The public repository is <https://github.com/sunsz2005/feasible-completion-y2h-panels>. The two yeast snapshots remain excluded pending redistribution permission. A release DOI has not yet been assigned; update the manuscript's availability statement with the repository URL and, when available, the release DOI.

