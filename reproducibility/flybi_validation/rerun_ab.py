"""Fresh timed A/B rerun on one archived FlyBi split; writes to reruns/."""
from pathlib import Path
import argparse, json, os, sys, time

# Set native-library limits before importing NumPy or any solver module.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
import balanced_assay_pilot as pilot
import balanced_assay_initialization_time as joint
from balanced_assay_tradeoff import feasible_role_greedy

parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=int, required=True, choices=range(20261113,20261133))
args = parser.parse_args()
train = np.load(ROOT / "splits" / f"{args.seed}_train.npy")
names = json.loads((ROOT / "names.json").read_text(encoding="utf-8"))
graph = pilot.matrix(len(names),train)
started = time.perf_counter()
selected = pilot.coverage_order(graph,200)
reference = pilot.partition(graph,selected)
reference_coverage = reference.coverage()
floor = .98*reference_coverage
feasible_role_greedy(graph,selected,floor)
t0 = time.perf_counter()-started
budget = 5*t0
out = ROOT / "reruns" / str(args.seed)
out.mkdir(parents=True,exist_ok=True)
for method in ("A","B"):
    panel, timing = joint.pipeline(graph,method,floor,budget,args.seed)
    db = np.flatnonzero(panel.roles==1).tolist()
    ad = np.flatnonzero(panel.roles==2).tolist()
    assert len(db)==len(ad)==100 and not set(db)&set(ad)
    assert panel.coverage()+1e-7 >= floor
    record = {"seed":args.seed,"method":method,"t0_seconds":t0,"reference_coverage":reference_coverage,"floor":floor,"db":db,"ad":ad,"train_cross":int(panel.score),"coverage":panel.coverage(),**timing}
    (out/f"panel_{method}.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(method, "train", panel.score, "timing_valid", timing["timing_valid"])
