#!/usr/bin/env bash
# Step 01: differential abundance with FragPipe-Analyst, run from the command line
# (no Shiny app) via scripts/00_common/run_FragPipe-Analyst_headless.R.
#
# Settings of the original GUI run (no log was saved; inferred from the exported
# matrices and confirmed by reproducing DE_results.csv exactly):
#   LFQ, MaxLFQ intensities, contaminants removed, no missing-value filter,
#   VSN normalisation, Perseus-type ("man") imputation, limma all pairwise contrasts,
#   BH, significant if p.adj <= 0.05 and |log2FC| >= 1.
#
# Run in the folder with combined_protein_f.tsv (step 00) and the annotation files
# from metadata/03_speckle_ActD/.
#   REPO = path to the NuFANCI_MS_Proteomics folder of this repository
#   FPA  = path to FragPipe-Analyst (git clone, commit 3e1af92)

REPO=${REPO:-../NuFANCI_MS_Proteomics}
FPA=${FPA:-../FragPipe-Analyst}
RUN="Rscript ${REPO}/scripts/00_common/run_FragPipe-Analyst_headless.R --fpa_dir=${FPA} --input=combined_protein_f.tsv"

## DMSO / ActD, NuFANCI / input, all pairwise contrasts -> DE_results.csv
$RUN --annotation=experiment_annotation.tsv --imputation=man
