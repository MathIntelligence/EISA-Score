# EISA-Score

There are two folder (a) *src* and (b) *examples*. The *src* folder contains the source code that can be used to generate EISA features with local and global surface area methods described in the paper. The surface area calculations require the evaluation of the atomic density function at each grid point making it computationally expensive. Therefore, it is recommended to use parallel computation to generate EISA features for each molecule in a dataset. The code *eisa_score_local.py* and *eisa_score_global.py* returns EISA features for local and global surface, respectively, for a given protein-ligand complex PDB id. Examples code for generating features for a dataset are provided in the *examples* folder.

## Steps to generate EISA features for PDBbind v2007 refined set

