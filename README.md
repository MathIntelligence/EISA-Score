# EISA-Score

There are three folders (a) `src`, (b) `examples`, and (c) `utils`. The `src` folder contains the source code that can be used to generate EISA features with local and global surface area methods described in the paper. The surface area calculations require the evaluation of the atomic density function at each grid point making it computationally expensive. Therefore, it is recommended to use parallel computation to generate EISA features for each molecular complex in a dataset. The code `eisa_score_local.py` and `eisa_score_global.py` returns EISA features for local and global surface, respectively, for a given protein-ligand complex PDBID. Examples code for generating features for a dataset are provided in the `examples` folder.

## Steps to generate EISA features

One can use the script `generate_eisa_features.py` from `examples` folder to generate EISA features for a given data set.

```shell
python generate_eisa_features.py --surface_type='global' --csv_file_path='../utils/PDBbindv2007_RefinedSet.csv' --pdb_index=0
```
Note: to generate features for all complexes in the dataset, run the script for `--pdb_index` from 0 to the length of the dataset.
