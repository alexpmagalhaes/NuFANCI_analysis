#!/usr/bin/env bash
# Step 00: remove contaminants and single-peptide proteins from the FragPipe protein table.
# Run in the folder with combined_protein.tsv (from PRIDE PXD069855). Writes combined_protein_f.tsv.
# REPO = path to the NuFANCI_MS_Proteomics folder of this repository.

REPO=${REPO:-../NuFANCI_MS_Proteomics}

python ${REPO}/scripts/00_common/filter_proteins.py \
  -i combined_protein.tsv \
  -c ${REPO}/reference_lists/Contaminants.tsv
