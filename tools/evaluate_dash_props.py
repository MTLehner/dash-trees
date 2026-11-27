"""Compare populated DASH trees with the reference values of the test molecules.

  python evaluate_dash_props.py <data folder> out.json name=tree.dash:default name=tree.dash:legacy

The test molecules are the 10 % that build_dash_props.py --split train leaves
out. The script reads each tree with the tie-breaking rule after its colon,
'default' or 'legacy'. For each tree and property it reports the RMSE, MAE and
R2 against the reference value of every conformer, and the RMSE against the
conformer mean per atom. 'conformerFloor' is the scatter of the conformers
around that mean, which no model that reads only the molecular graph can beat.
'charges' are the MBIS charges after the std-weighted normalization of
GetPartialCharges().

For the first two trees, it also reports a paired bootstrap of the RMSE of the
first tree minus the RMSE of the second (95 % interval, resampled over
molecules). It also counts the atoms that get a different MBIS charge.
'symmetry' counts the classes of symmetry-equivalent atoms whose raw MBIS
values ('result') differ by more than 1e-4 e.
"""

import json
import os
import sys

import numpy as np

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDASHTree

from build_dash_props import PROPERTIES, Dataset

THREADS = 0  # every hardware thread
USAGE = __doc__.split("\n\n")[1]


def parseSpecs(specs):
  """[(name, path, rule)] for the name=path:rule arguments."""
  parsed = []
  for spec in specs:
    name, _, rest = spec.partition("=")
    path, _, rule = rest.rpartition(":")
    if not name or name == "comparison" or rule not in ("default", "legacy") or \
        not os.path.isfile(path):
      sys.exit(f"bad tree argument {spec!r}, use name=tree.dash:default or :legacy\n{USAGE}")
    parsed.append((name, path, rule))
  names = [name for name, _, _ in parsed]
  if not parsed or len(set(names)) != len(names):
    sys.exit(f"give one or more trees, each with its own name\n{USAGE}")
  return parsed


def readParams(rule):
  params = rdDASHTree.DASHParams()
  params.legacyTieBreaking = rule == "legacy"
  return params


def errors(preds, refs):
  """The metrics, and the squared-error sum and count per molecule."""
  sse, cnt = np.zeros(len(refs)), np.zeros(len(refs))
  sae = sseMean = cntMean = floor = 0.0
  missing = 0
  ys, rs = [], []
  for i, (y, ref) in enumerate(zip(preds, refs)):
    y = np.asarray(y, dtype=np.float64)
    e = ref - y
    has = ~np.isnan(e)
    missing += int(np.sum(~np.isnan(ref) & np.isnan(y)))
    sse[i], cnt[i], sae = np.sum(e[has]**2), has.sum(), sae + np.abs(e[has]).sum()
    ys.append(np.broadcast_to(y, ref.shape)[has])
    rs.append(ref[has])
    if not has.any():
      continue
    mean = np.nanmean(ref[:, has.any(0)], axis=0)
    sseMean += np.sum((mean - y[has.any(0)])**2)
    cntMean += mean.size
    floor += np.nansum((ref[:, has.any(0)] - mean)**2)
  n = cnt.sum()
  if n == 0:
    return dict(n=0, noPrediction=missing), sse, cnt
  ys, rs = np.concatenate(ys), np.concatenate(rs)
  metrics = dict(n=int(n), noPrediction=missing, rmse=float(np.sqrt(sse.sum() / n)),
                 mae=float(sae / n),
                 r2=float(1 - np.sum((rs - ys)**2) / np.sum((rs - rs.mean())**2)),
                 rmseVsConformerMean=float(np.sqrt(sseMean / cntMean)),
                 conformerFloor=float(np.sqrt(floor / n)))
  return metrics, sse, cnt


def bootstrap(a, b, reps=2000):
  """The RMSE of a minus the RMSE of b, and its 95 % interval over molecules."""
  (sseA, cntA), (sseB, cntB) = a, b
  pick = np.random.default_rng(0).integers(0, len(sseA), size=(reps, len(sseA)))
  rmseA = np.sqrt(sseA[pick].sum(1) / cntA[pick].sum(1))
  rmseB = np.sqrt(sseB[pick].sum(1) / cntB[pick].sum(1))
  return dict(delta=float(np.sqrt(sseA.sum() / cntA.sum()) - np.sqrt(sseB.sum() / cntB.sum())),
              ci95=[float(np.percentile(rmseA - rmseB, q)) for q in (2.5, 97.5)])


def symmetry(mols, charges, tol=1e-4):
  classes = broken = 0
  for mol, q in zip(mols, charges):
    q = np.asarray(q)
    ranks = np.array(Chem.CanonicalRankAtoms(mol, breakTies=False))
    for r in np.unique(ranks):
      members = q[ranks == r]
      if members.size > 1:
        classes += 1
        broken += int(members.max() - members.min() > tol)
  return dict(classes=classes, broken=broken)


def main(folder, out, specs):
  trees = parseSpecs(specs)
  RDLogger.DisableLog("rdApp.*")
  mols, refs = [], {p: [] for p in PROPERTIES}
  for mol, values in Dataset(folder).records("test"):
    mols.append(mol)
    for p in PROPERTIES:
      refs[p].append(values[p])
  print(f"{len(mols)} test molecules", flush=True)

  results, perMolecule, mbis = {}, {}, {}
  for name, path, rule in trees:
    tree = rdDASHTree.DASHTree(path)
    metrics = {}
    for p in PROPERTIES:
      if p in tree.GetAvailablePropertyNames():
        preds = tree.GetMolPropertyBatch(mols, p, readParams(rule), THREADS)
        metrics[p], *perMolecule[(name, p)] = errors(preds, refs[p])
        if p == "result":
          metrics["symmetry"] = symmetry(mols, preds)
    options = rdDASHTree.NormalizationOptions()
    options.params = readParams(rule)
    mbis[name] = tree.GetPartialChargesBatch(mols, options, THREADS)
    metrics["charges"], *perMolecule[(name, "charges")] = errors(mbis[name], refs["result"])
    results[name] = dict(path=path, rule=rule, nodes=tree.GetNumNodes(), metrics=metrics)
    print(name, {p: round(m["rmse"], 5) for p, m in metrics.items() if "rmse" in m}, flush=True)
    del tree

  if len(trees) > 1:
    a, b = trees[0][0], trees[1][0]
    shared = [
      p for p in PROPERTIES + ["charges"] if (a, p) in perMolecule and (b, p) in perMolecule
    ]
    dq = np.concatenate([np.subtract(x, y) for x, y in zip(mbis[a], mbis[b])])
    results["comparison"] = dict(
      trees=[a,
             b], bootstrap={p: bootstrap(perMolecule[(a, p)], perMolecule[(b, p)])
                            for p in shared},
      charges=dict(atoms=int(dq.size), differ=int((abs(dq) > 1e-6).sum()), over0_01=int(
        (abs(dq) > 0.01).sum()), maxAbs=float(abs(dq).max()), rms=float(np.sqrt(np.mean(dq**2)))))
  with open(out, "w") as f:
    json.dump(results, f, indent=1)
  print("wrote", out)


if __name__ == "__main__":
  if len(sys.argv) < 4:
    sys.exit(USAGE)
  main(sys.argv[1], sys.argv[2], sys.argv[3:])
