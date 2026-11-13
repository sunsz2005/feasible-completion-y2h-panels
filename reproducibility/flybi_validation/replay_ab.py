"""Verify FlyBi A/B panels and independently recount training/evaluation edges."""
from pathlib import Path
import hashlib, json, sys
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
import balanced_assay_pilot as pilot

def load(path): return json.loads(path.read_text(encoding="utf-8"))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def stats(values):
    a = np.asarray(values, dtype=float)
    return (float(np.median(a)), float(np.quantile(a,.25)), float(np.quantile(a,.75)), int((a>0).sum()), int((a==0).sum()), int((a<0).sum()))

seal = load(ROOT / "pre_evaluation_seal_ab.json")
assert len(seal["files"]) == 60
for rel, digest in seal["files"].items(): assert sha(ROOT / rel) == digest, rel
split_manifest = load(ROOT / "split_manifest.json")
for rel, digest in split_manifest.items(): assert sha(ROOT / rel.replace("\\", "/")) == digest, rel
names = load(ROOT / "names.json")
published = load(ROOT / "scores_ab.json")
replayed = []
for seed in range(20261113,20261133):
    train = np.load(ROOT / "splits" / f"{seed}_train.npy")
    evaluation = np.load(ROOT / "splits" / f"{seed}_evaluation.npy")
    assert not (set(map(tuple,train)) & set(map(tuple,evaluation)))
    graph = pilot.matrix(len(names),train)
    cal = load(ROOT / "results" / f"calibration_{seed}.json")
    for method in ("A","B"):
        panel = load(ROOT / "results" / f"panel_{seed}_{method}.json")
        db, ad = panel["db"], panel["ad"]
        assert len(db)==len(ad)==100 and len(set(db))==len(set(ad))==100 and not set(db)&set(ad)
        chosen = db+ad
        counts = graph[:,chosen].sum(axis=1)
        terms = 1-np.power(.2,counts)
        terms[chosen] = 1
        coverage = float(terms.sum())
        assert abs(coverage-panel["coverage"]) < 1e-6
        assert coverage+1e-7 >= cal["floor"]
        train_cross = int(graph[np.ix_(db,ad)].sum())
        assert train_cross == panel["train_cross"]
        roles = np.zeros(len(names),dtype=np.int8)
        roles[db], roles[ad] = 1, 2
        eval_cross = int(( ((roles[evaluation[:,0]]==1)&(roles[evaluation[:,1]]==2)) | ((roles[evaluation[:,0]]==2)&(roles[evaluation[:,1]]==1)) ).sum())
        expected = next(r for r in published["rows"] if r["seed"]==seed and r["method"]==method)
        assert expected["train_cross"] == train_cross and expected["evaluation_cross"] == eval_cross
        assert panel["timing_valid"] and expected["timing_valid"]
        replayed.append((seed, method, train_cross, eval_cross))
delta = [next(r[3] for r in replayed if r[0]==seed and r[1]=="B")-next(r[3] for r in replayed if r[0]==seed and r[1]=="A") for seed in range(20261113,20261133)]
assert stats(delta) == (2.5,-.25,9.25,14,1,5)
print("Verified 40 A/B panels, timing flags, scores, and 20 paired evaluation differences.")
