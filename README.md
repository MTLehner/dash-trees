# DASH trees for the RDKit

Data files for the RDKit's DASH module (`Code/GraphMol/DASHTree`,
`rdkit.Chem.DASHTree`), kept outside the RDKit repository. Each file is a
`.dash` container; `rdkit.Chem.DASHTree.GetDASHTree(name)` downloads it by name.

| name | file | md5 |
|---|---|---|
| `initial_test` | `initial_test.dash` (38 MB) | `3dc6a5c2f221e7a3be2edea1b38518b3` |

The RDKit build downloads `initial_test` into `Data/DASHTree` for its tests, and
`GetDASHTree()` uses it when no other tree is given. Files are fetched from a
tag (`v1`), so a published file never changes under its checksum.

## initial_test

The MBIS partial-charge tree of the
[DASH-tree](https://github.com/rinikerlab/DASH-tree) package
([MIT licence](https://github.com/rinikerlab/DASH-tree/blob/main/LICENSE)),
version 2.0.1 (commit `6cf1b23`, `serenityff/charge/data/default_dash_tree`),
converted and pruned with the RDKit's tools:

```
python -m rdkit.Chem.DASHTree.convert <default_dash_tree> full.dash
python -m rdkit.Chem.DASHTree.prune full.dash initial_test.dash --tolerance 0.02 --props result,std
```

Both steps are deterministic, so the file is reproducible byte for byte
(sha256 `8c9e35f69bf9d81cf0fb19922105cddd045bf5b4def7e18461fec7bb845a3e1b`). It
carries the `result` and `std` columns; on the 100 reference molecules in the
RDKit's test data its charges stay within 0.044 e of the full tree's
(0.0068 e RMS).

References:
- M. Lehner et al., J. Chem. Inf. Model. 2023, 63, 6014
- M. Lehner et al., J. Chem. Phys. 2024, 161, 074103
