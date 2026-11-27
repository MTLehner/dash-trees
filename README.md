# DASH trees for the RDKit

This repository holds data files for the RDKit's DASH module
(`Code/GraphMol/DASHTree`, `rdkit.Chem.DASHTree`). The RDKit repository does
not contain them. Each file is a `.dash` container, and
`rdkit.Chem.DASHTree.GetDASHTree(name)` downloads it by name into a local
cache. The RDKit build also downloads `initial_test` into `Data/DASHTree` for
its tests, and `GetDASHTree()` opens that one when you name no other tree.

The RDKit downloads every file from a tag, so a published file never changes
under its checksum.

## Available trees

| name | size | nodes | charge error vs the full tree (RMS / max) | tag | md5 |
|---|---|---|---|---|---|
| [`mbis_pruned_0.01`](#mbis_pruned_001) | 63.7 MB | 3,981,565 | 0.0036 e / 0.036 e | `v2` | `b5f8ce1d7b9bd71c80c73983b9b7a15b` |
| [`initial_test`](#initial_test) | 38.2 MB | 2,385,193 | 0.0068 e / 0.044 e | `v1` | `3dc6a5c2f221e7a3be2edea1b38518b3` |
| [`mbis_pruned_0.1`](#mbis_pruned_01) | 4.2 MB | 259,275 | 0.024 e / 0.10 e | `v2` | `e5264d130cb2c4534d5ceb90914a91dc` |

The errors compare std-weighted partial charges with those of the full,
unpruned tree (9,439,485 nodes) on the 100 reference molecules in the RDKit's
`Code/GraphMol/DASHTree/test_data` (5019 atoms).

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

## mbis_pruned_0.01

This is the most accurate of the three, pruned at `--tolerance 0.01`. Use it
when the charges should stay as close as possible to the full tree's. It is the
largest (sha256
`16445710f5a594964ada60a8087c2f152f57d534e3210543b2c0b41fe7525eda`).

## initial_test

This is the default, pruned at `--tolerance 0.02`. The RDKit's tests run on it,
and `GetDASHTree()` opens it when you name no other tree. It balances size and
accuracy (sha256
`8c9e35f69bf9d81cf0fb19922105cddd045bf5b4def7e18461fec7bb845a3e1b`).

## mbis_pruned_0.1

This is the smallest, pruned at `--tolerance 0.1`. It is for quick estimates,
and for environments where a 4 MB download matters more than the last
hundredths of an electron (sha256
`140463ca3603256f38e7eabe1d3f9e6ae40d50b1eddc14b89e542d22350f06da`).

## References

- M. Lehner et al., J. Chem. Inf. Model. 2023, 63, 6014
- M. Lehner et al., J. Chem. Phys. 2024, 161, 074103
