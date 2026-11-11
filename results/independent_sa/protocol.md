# Independent simulated-annealing benchmark protocol

Frozen before scoring, 2026-11-11. This is a new follow-up comparison; the
original 360-panel A/B experiment and its archived results remain untouched.

Use the same 20 training/evaluation edge splits for each yeast source and the
same fixed human training/evaluation graph pair with 20 search seeds. Compare
three methods at the archived instance-specific primary total budget `5*T0`:
A (coverage-first plus the frozen common multistart solver), B (completion-set
start plus that solver), and a separately implemented simulated-annealing (SA)
joint optimizer. Rerun all three in the same Python/NumPy environment on the
same machine, sequentially and with rotating method order. Graph loading and
matrix construction are common offline costs, excluded for all three. Every
method's coverage-reference construction and role initialization count against
its own total budget. The original calibration fixes only `T0` and the coverage
floor; no archived A/B panel is reused.

SA starts from the coverage-reference panel after the same fixed-panel role
partition used by A. It never invokes the common multistart search or the
completion-set constructor. At each proposal it chooses, with equal
probability, a DB-AD role swap or a same-role selected-to-unselected protein
replacement. For a replacement, choose the selected protein uniformly from
one of the two roles; choose the outside protein from the top 1,000
training-degree candidates with probability 0.5 and uniformly from all
candidates otherwise. Reject any replacement violating the exact training
coverage floor. Accept a feasible move with probability
`min(1, exp(training_delta / T))`. Use four geometric cooling cycles between
`T=4.0` and `T=0.2` over the available total wall-clock budget. Return the
best feasible training-score panel encountered. No evaluation edge or score is
read during generation. The same 0.01-second return reserve and 5% or 0.02-s
timing tolerance apply to all methods.

The comparison reports all 20 paired differences per source, separately for
B versus SA and A versus SA, as well as training scores, coverage, runtime,
invalid timing paths, and the full ranges. Medians of paired differences are
not inferred from marginal medians. Because the networks and model choices
were used during method development, these are computational comparisons on
fixed backgrounds, not an independent biological validation or evidence of
universal superiority over other optimizers. The human graph assumes
dual-role eligibility. No algorithm is selected or tuned using evaluation
edges. All 180 panels must be hashed and sealed before scoring begins.

The run requires local YDIP and YBIOGRID snapshots, which are excluded from
the public repository pending redistribution permission. Their hashes are
checked by the experiment lock. The public release contains the restoration
script and source provenance, not the snapshots themselves.
