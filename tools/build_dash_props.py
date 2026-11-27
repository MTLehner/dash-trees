"""Populate a DASH tree with the DASH-properties data set.

The data set is "Calculated Atomic and Molecular Properties for the DASH
Dataset", ETH Research Collection, doi:10.3929/ethz-b-000670546 (CC BY-SA
4.0). Download these files into one folder. Unpack the zip files there:

  sdf_qmugs500_mbis_collect.sdf.zip   871,906 conformers of 293,822 molecules
  props.h5                            per-atom wave-function properties
  am1bcc.zip                          per-atom AM1-BCC charges, which unpack
                                      to am1bcc/charges/*.csv

The tree comes from the full MBIS charge tree of the DASH-tree package,
converted with python -m rdkit.Chem.DASHTree.convert. The script keeps the
nodes of the tree and replaces its columns. For each property P, each node
gets the mean (P), the sample standard deviation (P_std) and the number of
values (P_count) of the atoms that reach it. Each conformer adds one value per
atom. AM1-BCC is the exception: each molecule adds one value per atom. The MBIS
charge is P = 'result', with 'std' and 'size' for the other two columns, so
GetPartialCharges() reads it by default.

  python build_dash_props.py <data folder> <full.dash> <out.dash>

Use --legacy-tie-breaking for a tree that you read with
DASHParams.legacyTieBreaking. --split train leaves out about one molecule in
ten (chosen by a hash of its canonical SMILES), so that evaluate_dash_props.py can test on the
molecules that the tree did not see. The script needs an RDKit with
rdkit.Chem.DASHTree, and numpy, pandas and h5py.
"""

import argparse
import glob
import hashlib
import os

import h5py
import numpy as np
import pandas as pd

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDASHTree
from rdkit.Chem.DASHTree import populate

SDF = "sdf_qmugs500_mbis_collect.sdf"
# columns of props.h5 that are per-atom properties, with their names in the tree
H5_COLUMNS = {
  "mulliken": "mulliken",
  "resp1": "resp1",
  "resp2": "resp2",
  "dual": "dual",
  "mbis_dipole_strength": "mbis_dipole_strength"
}
PROPERTIES = ["result", *H5_COLUMNS.values(), "am1bcc"]
RENAME = {"result_std": "std", "result_count": "size"}


class Dataset:
  """The SDF records, their atoms, and every property per atom, in SDF order."""

  def __init__(self, folder):
    self.folder = folder
    self._indexSdf()
    self._groupConformers()
    self.values = {"result": self.mbis}
    self._readH5()
    self._readAm1bcc()

  def _indexSdf(self):
    """One pass over the SDF: where each record starts, its atom count and
    DASH_IDX, and per atom its element and MBIS charge."""
    offsets, nAtoms, dash, elements, mbis = [0], [], [], [], []
    pos, lineNo, want, n = 0, 0, False, 0
    with open(os.path.join(self.folder, SDF), "rb") as f:
      for line in f:
        pos += len(line)
        if lineNo == 3:
          n = int(line[:3])
          nAtoms.append(n)
        elif 3 < lineNo <= 3 + n:
          elements.append(line[31:34].strip())
        lineNo += 1
        if want:
          field, want = want, False
          if field == b"MBIScharge":
            mbis.append(np.array(line.split(b"|"), dtype=np.float64))
          else:
            dash.append(line.strip().decode())
        elif line.startswith(b">  <MBIScharge>") or line.startswith(b">  <DASH_IDX>"):
          want = line[4:line.index(b">", 4)]
        elif line.startswith(b"$$$$"):
          offsets.append(pos)
          lineNo = 0
    self.offsets = np.array(offsets, dtype=np.int64)
    self.nAtoms = np.array(nAtoms, dtype=np.int64)
    self.start = np.concatenate([[0], np.cumsum(self.nAtoms)])
    self.dash = np.array(dash)
    self.elements = np.array(elements, dtype="S2")
    self.mbis = np.concatenate(mbis)
    if not len(dash) == len(mbis) == self.nAtoms.size or self.mbis.size != self.start[-1]:
      raise ValueError("every record needs a DASH_IDX and one MBIS charge per atom")
    print(f"{self.nAtoms.size} records, {self.start[-1]} atoms", flush=True)

  def _atomsOf(self, r):
    return self.elements[self.start[r]:self.start[r + 1]].tobytes()

  def _groupConformers(self):
    """A molecule is a run of records with the same DASH_IDX and the same atoms."""
    first = [0]
    for r in range(1, self.nAtoms.size):
      if self.dash[r] != self.dash[r - 1] or self._atomsOf(r) != self._atomsOf(r - 1):
        first.append(r)
    self.first = np.array(first)
    self.nConf = np.diff(np.append(self.first, self.nAtoms.size))
    print(f"{self.first.size} molecules", flush=True)

  def _readH5(self):
    """props.h5 is a pandas frame with one row per atom, in SDF order. h5py
    reads its numeric blocks only, so no pickled object column is loaded."""
    frame = h5py.File(os.path.join(self.folder, "props.h5"), "r")["df"]
    names = [x.decode() for x in frame["block0_items"][:]]
    index = [x.decode() for x in frame["block1_items"][:]]
    atomIdx = frame["block1_values"][:, index.index("atom_idx")]
    inRecord = np.arange(self.start[-1]) - np.repeat(self.start[:-1], self.nAtoms)
    charges = frame["block0_values"][:, names.index("MBIScharge")]
    if not (np.array_equal(atomIdx, inRecord) and np.array_equal(charges, self.mbis)):
      raise ValueError("props.h5 is not in the atom order of the SDF")
    for column, name in H5_COLUMNS.items():
      self.values[name] = frame["block0_values"][:, names.index(column)]

  def _readAm1bcc(self):
    """am1bcc/charges/*.csv has one conformer-averaged value per atom and
    molecule. It goes on the molecule's first conformer, so it counts once."""
    folder = os.path.join(self.folder, "am1bcc", "charges")
    paths = glob.glob(os.path.join(folder, "*.csv"))
    if not paths:
      raise FileNotFoundError(f"no AM1-BCC tables in {folder}")
    firstOf = {}
    for d, r in zip(self.dash[self.first], self.first):
      firstOf.setdefault(d, r)
    values = np.full(self.start[-1], np.nan)
    placed = rejected = 0
    for path in paths:
      table = pd.read_csv(
        path, usecols=["mol_index", "dash_index", "atom_index", "element", "am1Bcc_averaged"])
      for _, rows in table.groupby("mol_index", sort=False):
        r = firstOf.get(rows.dash_index.iloc[0])
        atoms = rows.atom_index.to_numpy()
        if r is None or len(rows) != self.nAtoms[r] or atoms.max() >= self.nAtoms[r] or \
            not np.array_equal(self.elements[self.start[r] + atoms],
                               rows.element.to_numpy().astype("S2")):
          rejected += 1
          continue
        values[self.start[r] + atoms] = rows.am1Bcc_averaged.to_numpy()
        placed += 1
    if not placed:
      raise ValueError(f"no AM1-BCC table in {folder} matches a molecule of the SDF")
    print(f"AM1-BCC: {placed} conformer tables placed, {rejected} without a matching molecule",
          flush=True)
    self.values["am1bcc"] = values

  def molecules(self):
    """(g, the first conformer of molecule g) for every molecule, parsed from
    its mol block alone."""
    with open(os.path.join(self.folder, SDF), "rb") as f:
      for g, r in enumerate(self.first):
        f.seek(self.offsets[r])
        block = f.read(self.offsets[r + 1] - self.offsets[r]).split(b"M  END", 1)[0]
        yield g, Chem.MolFromMolBlock(block.decode() + "M  END\n", removeHs=False)

  def records(self, split="all"):
    """(mol, {property: one row per conformer}) for each molecule of the split."""
    for g, mol in self.molecules():
      if mol is None or (split != "all" and isTest(mol) != (split == "test")):
        continue
      s, e = self.start[self.first[g]], self.start[self.first[g] + self.nConf[g]]
      yield mol, {p: v[s:e].reshape(self.nConf[g], -1) for p, v in self.values.items()}


def isTest(mol):
  """About one molecule in ten, chosen by a hash of its canonical SMILES."""
  key = Chem.MolToSmiles(Chem.RemoveHs(mol, sanitize=False), isomericSmiles=False)
  return int.from_bytes(hashlib.blake2b(key.encode(), digest_size=8).digest(), "little") % 10 == 0


def main():
  ap = argparse.ArgumentParser(description=__doc__,
                               formatter_class=argparse.RawDescriptionHelpFormatter)
  ap.add_argument("data", help="the folder with the unpacked data set")
  ap.add_argument("tree", help="the full MBIS tree as a .dash container")
  ap.add_argument("out", help="the .dash file to write")
  ap.add_argument("--split", choices=["all", "train"], default="all",
                  help="train: leave out the test molecules of evaluate_dash_props.py")
  ap.add_argument("--legacy-tie-breaking", action="store_true",
                  help="match atoms as the DASH-tree python package does")
  args = ap.parse_args()
  if os.path.exists(args.out) and os.path.samefile(args.tree, args.out):
    ap.error("the output would overwrite the input")
  RDLogger.DisableLog("rdApp.*")
  data = Dataset(args.data)
  params = rdDASHTree.DASHParams()
  params.legacyTieBreaking = args.legacy_tie_breaking
  columns, failed = populate.collectStatistics(args.tree, data.records(args.split), PROPERTIES,
                                               params)
  columns = [(RENAME.get(name, name), tag, arr) for name, tag, arr in columns]
  populate.writeTree(args.tree, args.out, columns, keepExisting=False)
  print(f"{failed} molecules with an atom outside the tree's classes; wrote {args.out}")


if __name__ == "__main__":
  main()
