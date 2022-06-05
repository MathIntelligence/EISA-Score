# EISA-Score

There are three folders (a) `src`, (b) `examples`, and (c) `utils`. The `src` folder contains the source code that can be used to generate EISA features with local and global surface area methods described in the paper. The surface area calculations require the evaluation of the atomic density function at each grid point making it computationally expensive. Therefore, it is recommended to use parallel computation to generate EISA features for each molecular complex in a dataset. The code `eisa_score_local.py` and `eisa_score_global.py` returns EISA features for local and global surface, respectively, for a given protein-ligand complex PDBID. Examples code for generating features for a dataset are provided in the `examples` folder. The `utils` folder contains the csv files of the dataset containing 'PDBID' and 'pK' values.

## Steps to generate EISA features

One can use the scripts `get_eisa_features_global.py` and `get_eisa_features_local.py`from `examples` folder to generate EISA features for a given data set. 

```shell

# Generate EISA features using global surface area method

python get_eisa_features_global.py --dataset_csv_file '../utils/PDBbindv2016_RefinedSet.csv' --data_folder '../PDBbindDataset/v2016_refined_set' --out_dir '../features' --kernel_type 'exponential' --kernel_tau 1.0 --kernel_power 3.0 --cutoff 12.0 --pdbid_index 0

# Generate EISA features using local surface area method

python get_eisa_features_local.py --dataset_csv_file '../utils/PDBbindv2016_RefinedSet.csv' --data_folder '../PDBbindDataset/v2016_refined_set' --out_dir '../features' --kernel_type 'exponential' --kernel_tau 0.5 --kernel_power 15.0 --cutoff 6.5 --isovalue 0.15 --pdbid_index 0

```
Note: to generate features for all complexes in the dataset, run the script for `--pdb_index` from 0 to the length of the dataset.
