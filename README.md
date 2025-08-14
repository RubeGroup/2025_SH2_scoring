# 2025_SH2_Models
This repository contains a script for scoring peptide sequences using the free energy models for SH2 domains introduced in [Gagoski & Rube et al. (2024)](https://www.biorxiv.org/content/10.1101/2024.12.23.630085v2.full). For example, a peptide sequence (such as `LGHQRYHNITPGA`) can be scored using a SH2 binding model of interest (such as `Blk`) by executing
```
./scoreProteinSequences.py -d Blk -s LGHQRYHNITPGA
```
By default, the script reports the sum of the predicted relative affinities at all available binding windows, but prediction at each window can be retrieved using the argument `--profile`. 

Many sequences can be scored using the arguments `-t file.txt` (where `file.txt` contains one sequence per line) or `-f file.fa` (where `file.fa` is a FASTA file). 

By default, the script uses the binding models for  `c-Src`, `Fyn`, `Grb2`, `Lyn`, `Yes`, and `Blk` shown in figures 4 and 5, but the other models (identified using the `fit_id` values listed in the supplemental table and on the dryad repository (https://datadryad.org/dataset/doi:10.5061/dryad.msbcc2g7w)) can be selected using the argumnet `-i fit_id`. In addition, a custom scoring model JSON object (stored in the a file `model.json`) can be selecte using the argument `-j model.json`.
