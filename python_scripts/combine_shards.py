"""Combine independent BoxCorr ensemble shards of the same (x,T) point into one
high-statistics estimate.  Usage: combine_shards.py SHARD_DIR...
Each shard's box_corr.txt has one row: x  MC  stderr  (for one T)."""
import json, os, sys
import numpy as np

dirs = sys.argv[1:]
means, ses, ns = [], [], []
x = T = None
for d in dirs:
    raw = np.atleast_2d(np.loadtxt(os.path.join(d, "box_corr.txt")))
    cfg = json.load(open(os.path.join(d, "configuration.json")))
    means.append(raw[0, 1]); ses.append(raw[0, 2]); ns.append(cfg["nsims"])
    x = raw[0, 0]; T = cfg["T_values"][0]
means, ses, ns = np.array(means), np.array(ses), np.array(ns)

# equal nsims => simple average; combined variance = (1/n^2) sum se_i^2
n = len(dirs)
m = means.mean()
se_prop = np.sqrt((ses ** 2).sum()) / n            # propagated from per-shard se
se_emp = means.std(ddof=1) / np.sqrt(n)            # empirical from shard spread (cross-check)

print(f"point x={x:g}  T={T:g}  shards={n}  total nsims={int(ns.sum())}")
print(f"per-shard means : {np.array2string(means, precision=6)}")
print(f"per-shard stderr: {np.array2string(ses, precision=6)}")
print(f"COMBINED mean   = {m:.8f}")
print(f"COMBINED stderr = {se_prop:.8f}   (propagated)")
print(f"  cross-check   = {se_emp:.8f}   (empirical shard spread; should agree)")
