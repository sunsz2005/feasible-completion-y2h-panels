# Release checks

1. Run `python tools/verify_release.py` in the download copy.
2. In the documented environments, replay primary results, coverage sensitivity and FlyBi A/B results. Keep locally restored inputs, compatibility copies and rerun outputs out of Git.
3. Preserve all reported A/B results, timing fallbacks and supplementary controls. Original panel hashes must not be regenerated to mask a change.
4. The sensitivity and FlyBi protocol extracts and hash subsets were packaged after scoring. Preserve this distinction when describing them.
5. At the final public update, record the new commit/release identifier and align the manuscript's repository citation with that version. The current local restructuring does not change the existing public commit.
6. Check source attribution and exclusions in `DATA_SOURCES.md`. MIT covers the project code, not the upstream datasets.

Two historical lock files are omitted. The primary replay reports this explicitly and uses the relative-path metadata in `results/primary/replay_metadata.json`. Retained hashes and panel records must still pass validation.
