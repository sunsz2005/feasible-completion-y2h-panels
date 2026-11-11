# Coverage-threshold sensitivity protocol

This follow-up varies only the training coverage floor at `tau = 0.90, 0.95,
0.98`. The training fraction is `q = 0.8`, total panel capacity is `k = 200`,
and every method has the same instance-specific `5*T0` wall-clock budget as
the primary experiment. All initialization time counts. The two new thresholds
rerun pipelines A and B on the original 20 seeds in each of YDIP, YBIOGRID,
and Human Space-I. The `0.98` paths are taken from the sealed independent-SA
rerun because it used the same A/B implementation, seeds, computer, software
environment, and `5*T0` budget. Those paths are reverified and rescored.

The method order alternates by source, seed, and threshold. All paths run
serially. The input graph and evaluation split are fixed within a matched
source/seed pair. No evaluation interaction is used to generate or choose a
panel. New panel files are SHA-256 sealed before evaluation interactions are
read. The same training coverage function is recomputed for every panel during
scoring. The source snapshots are excluded from the public repository; the
existing retrieval instructions apply.

To probe the local feasible replacement space without changing the timed
budget, the scorer draws 1,000 selected-to-unselected protein pairs from the
same coverage-reference panel for every source/seed. It evaluates each proposed
replacement under all three thresholds, reporting the feasible fraction. This
is a fixed-panel sample, not an exact count of all possible exchanges or a
measure of the states visited by the optimizer. The probe is run after the
timed panels have been sealed and is not charged to any method.

The human 20 runs share one historical training/evaluation graph pair and vary
only by search seed. All sensitivity results are descriptive; the thresholds
and backgrounds had already informed method development.
