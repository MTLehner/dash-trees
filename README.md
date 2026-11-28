# DASH trees for the RDKit

This repository holds data files for the RDKit's DASH module
(`Code/GraphMol/DASHTree`, `rdkit.Chem.DASHTree`). The RDKit repository does
not contain them. Each file is a `.dash` container, and
`rdkit.Chem.DASHTree.GetDASHTree(name)` downloads it by name into a local
cache. The RDKit build also downloads `initial_test` into `Data/DASHTree` for
its tests, and `GetDASHTree()` opens that one when you name no other tree.

The RDKit downloads each file in git from a tag, so such a file never changes
under its checksum. `dash_props` is too large for a file in git, so it is an
asset of the release `v4`. Compare its md5 with the table below.

## Available trees

The MBIS charge trees:

| name | size | nodes | charge error vs the full tree (RMS / max) | tag | md5 |
|---|---|---|---|---|---|
| [`mbis_pruned_0.01`](#mbis_pruned_001) | 63.7 MB | 3,981,565 | 0.0036 e / 0.036 e | `v2` | `b5f8ce1d7b9bd71c80c73983b9b7a15b` |
| [`initial_test`](#initial_test) | 38.2 MB | 2,385,193 | 0.0068 e / 0.044 e | `v1` | `3dc6a5c2f221e7a3be2edea1b38518b3` |
| [`mbis_pruned_0.1`](#mbis_pruned_01) | 4.2 MB | 259,275 | 0.024 e / 0.10 e | `v2` | `e5264d130cb2c4534d5ceb90914a91dc` |

The trees with seven properties:

| name | size | nodes | MBIS charge RMS error on held-out molecules ([Accuracy](#accuracy)) | tag | md5 |
|---|---|---|---|---|---|
| [`dash_props`](#dash_props) | 774.0 MB | 9,439,485 | 0.0250 e | `v4` (release asset) | `22444b59bbcd84e46b7a1172481d9024` |
| [`dash_props_pruned_0.02`](#dash_props_pruned_002) | 93.9 MB | 1,738,553 | 0.0254 e | `v4` | `3fc18d68f56ae34aacb9443079dc5f2e` |

## The MBIS charge trees

The errors in the first table compare std-weighted partial charges with those
of the full, unpruned tree (9,439,485 nodes) on the 100 reference molecules in
the RDKit's `Code/GraphMol/DASHTree/test_data` (5019 atoms).

All three are pruned from the MBIS partial-charge tree of the
[DASH-tree](https://github.com/rinikerlab/DASH-tree) package
([MIT licence](https://github.com/rinikerlab/DASH-tree/blob/main/LICENSE)),
version 2.0.1 (commit `6cf1b23`, `serenityff/charge/data/default_dash_tree`).
They carry that tree's `result` (charge) and `std` columns. They differ only in
the pruning tolerance, so the same two commands make all of them:

```
python -m rdkit.Chem.DASHTree.convert <default_dash_tree> full.dash
python -m rdkit.Chem.DASHTree.prune full.dash <name>.dash --tolerance <EPS> --props result,std
```

Both steps are deterministic, so each file is reproducible byte for byte.
Pruning drops a subtree when no value in it differs from its parent's by EPS
or more. An atom whose descent reaches a cut reads the deepest node that is
still there. So every tree assigns values to every molecule the full tree
accepts.

### mbis_pruned_0.01

This is the most accurate of the three, pruned at `--tolerance 0.01`. Use it
when the charges should stay as close as possible to the full tree's. It is the
largest (sha256
`16445710f5a594964ada60a8087c2f152f57d534e3210543b2c0b41fe7525eda`).

### initial_test

This is the default, pruned at `--tolerance 0.02`. The RDKit's tests run on it,
and `GetDASHTree()` opens it when you name no other tree. It balances size and
accuracy (sha256
`8c9e35f69bf9d81cf0fb19922105cddd045bf5b4def7e18461fec7bb845a3e1b`).

### mbis_pruned_0.1

This is the smallest, pruned at `--tolerance 0.1`. It is for quick estimates,
and for environments where a 4 MB download matters more than the last
hundredths of an electron (sha256
`140463ca3603256f38e7eabe1d3f9e6ae40d50b1eddc14b89e542d22350f06da`).

## The trees with seven properties

`dash_props` has the nodes of the full MBIS tree. `tools/build_dash_props.py`
filled them with the RDKit's populator (`rdkit.Chem.DASHTree.populate`) and the
DASH-properties data set:
"Calculated Atomic and Molecular Properties for the DASH Dataset", ETH Research
Collection, [doi:10.3929/ethz-b-000670546](https://doi.org/10.3929/ethz-b-000670546)
(871,906 conformers of 293,822 molecules). The populator matches each atom and
gives its value to every node on the atom's path. Each node then stores the
mean, the sample standard deviation and the number of values.

| column | property |
|---|---|
| `result`, `std`, `size` | MBIS partial charge (e) |
| `mulliken` | Mulliken charge (e) |
| `resp1`, `resp2` | the two RESP charges of the data set (e) |
| `am1bcc` | AM1-BCC charge (e) |
| `dual` | dual descriptor |
| `mbis_dipole_strength` | magnitude of the MBIS atomic dipole |

Each property except the MBIS charge has a `_std` and a `_count` column. Each conformer
adds one value per atom. AM1-BCC is the exception: the data set gives one
value per atom and molecule. The MBIS charge uses the column names that
`GetPartialCharges()` reads by default.

The data set also has DFT-D4 C6 values. Their molecule indices refer to an
earlier data set (doi:10.3929/ethz-b-000613415). That data set has a different
licence (CC BY-NC-SA 4.0), so these trees do not have the C6 values.

### dash_props

This is the full tree, with all 21 columns (sha256 `f96aa474e6984d85989bf049f3147bdda8ce850d5c508d6dd72368180f718c17`). Use it for
the most accurate values.

### dash_props_pruned_0.02

This is `dash_props` pruned at `--tolerance 0.02` on `result`. It has the means
and standard deviations, but no `_count` columns and no `size` column (sha256
`8c23769b7db87f720cb3646fb3d681e151d20c95f40ed437dd82d2dfe2292cba`). It is small enough for a file in git. Its MBIS charges are 1.5 %
less accurate. Its AM1-BCC charges are 2.3 % less accurate.

### How the trees were made

`tools/build_dash_props.py` reads the data set and populates the full MBIS tree
(`full.dash`, the output of the convert command above):

```
python tools/build_dash_props.py <data folder> full.dash dash_props.dash
python -m rdkit.Chem.DASHTree.prune dash_props.dash dash_props_pruned_0.02.dash --tolerance 0.02 --props result,std,mulliken,mulliken_std,resp1,resp1_std,resp2,resp2_std,dual,dual_std,mbis_dipole_strength,mbis_dipole_strength_std,am1bcc,am1bcc_std
```

Both steps are deterministic. The docstring of `build_dash_props.py` lists the
files of the data set that it reads. The trees of `v4` come from the populator
of [rdkit/rdkit#9795](https://github.com/rdkit/rdkit/pull/9795) at commit
`add5de390`, which computes the standard deviations without cancellation. The
tag `v3` has an earlier build of `dash_props_pruned_0.02`. Two of its standard
deviations (`mulliken_std` and `resp1_std`, one node each) differ by 2⁻²⁴, and
all its other values are the same.

### Accuracy

To measure the accuracy, `build_dash_props.py --split train` populated the tree
from 90 % of the molecules. `tools/evaluate_dash_props.py` then compared its
values with every conformer of the other 10 % (29,469 molecules, 1.23 million
atoms). The table gives the RMS errors:

| column | full tree | pruned at 0.02 e | conformer scatter |
|---|---|---|---|
| `result` (MBIS charge) | 0.0250 e | 0.0254 e | 0.0119 e |
| `result` after `GetPartialCharges()` | 0.0243 e | 0.0247 e | 0.0119 e |
| `mulliken` | 0.0374 e | 0.0377 e | 0.0241 e |
| `resp1` | 0.0676 e | 0.0673 e | 0.0420 e |
| `resp2` | 0.0699 e | 0.0683 e | 0.0423 e |
| `am1bcc` | 0.0278 e | 0.0285 e | 0 |
| `dual` | 0.0477 | 0.0479 | 0.0164 |
| `mbis_dipole_strength` | 0.0106 | 0.0107 | 0.0070 |

The conformer scatter is the RMS distance of the conformers' values from their
mean. No tree can be more accurate than that, because a tree reads the
molecular graph only.

The populator used the RDKit's default tie-breaking rule. So read these trees
with the default `DASHParams`, not with `legacyTieBreaking`. A tree populated
and read with the legacy rule has a 0.3 % larger MBIS charge error (0.02510 e
against 0.02504 e). A paired bootstrap over the molecules shows that the
difference is significant. For the other properties, the default rule is as
accurate or more accurate. With the default rule, symmetry-equivalent atoms
get the same values. With the legacy rule, 9 % of the classes of
symmetry-equivalent atoms get different values. If you populate a tree with
one rule and read it with the other, its MBIS charge error is 8 % to 42 %
larger.

## Licences

The MBIS charge trees come from the DASH-tree package under the MIT licence.
The tools are under the MIT licence of this repository. The `dash_props` trees
come from the DASH-properties data set, so they carry its licence,
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). If you use
them, cite the data set and the J. Chem. Phys. 2024 paper below.

## References

- M. Lehner et al., J. Chem. Inf. Model. 2023, 63, 6014
- M. Lehner et al., J. Chem. Phys. 2024, 161, 074103
- Calculated Atomic and Molecular Properties for the DASH Dataset, ETH Research
  Collection, 2024, doi:10.3929/ethz-b-000670546
