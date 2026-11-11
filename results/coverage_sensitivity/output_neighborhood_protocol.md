# Post hoc output-neighborhood probe

After both the new `tau=0.90/0.95` panels and the reused `tau=0.98` panels were
sealed, `tools/probe_output_neighborhood.py` sampled 1,000 selected-to-unselected
protein replacements from each of the 360 final A/B panels. It used the same
training network and that panel's own coverage threshold. The proposal index
stream is fixed by source and seed and reused across methods and thresholds;
actual protein pairs differ because the output panels differ. The exact
coverage function is recomputed for each hypothetical replacement. Three
vectorized results per panel are checked against the ordinary `State.swap`
implementation. The probe does not modify, select, or rerun any timed output.

The fraction estimates only the local feasible replacement neighborhood under
this uniform proposal rule. It does not count all exchanges, role swaps, or
states reachable through longer paths. Because the panels themselves differ
across methods and thresholds, changes in the fraction cannot be assigned
causally to the threshold or used to explain an evaluation-score gain.

The reference-panel probe in `scores.json` used the same 1,000-proposal rule
before any joint search. Nearly all replacements were feasible at those high
coverage starting panels, so it was not used as a performance explanation.
