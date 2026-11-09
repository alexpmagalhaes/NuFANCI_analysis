#!/usr/bin/env bash
# Step 01: differential abundance with FragPipe-Analyst, run from the command line
# (no Shiny app) via scripts/00_common/run_FragPipe-Analyst_headless.R.
#
# Settings = those of the original GUI run (see results/02_nucleolus_vs_speckle/
# differential_abundance/FragPipe-Analyst_report.pdf):
#   LFQ, MaxLFQ intensities, contaminants removed, no missing-value filter,
#   VSN normalisation, Perseus-type ("man") imputation, limma, BH,
#   significant if p.adj <= 0.05 and |log2FC| >= 1.
#
# Run in the folder with combined_protein_f.tsv (step 00) and the annotation files
# from metadata/02_nucleolus_vs_speckle/.
#   REPO = path to the NuFANCI_MS_Proteomics folder of this repository
#   FPA  = path to FragPipe-Analyst (git clone, commit 3e1af92)

REPO=${REPO:-../NuFANCI_MS_Proteomics}
FPA=${FPA:-../FragPipe-Analyst}
RUN="Rscript ${REPO}/scripts/00_common/run_FragPipe-Analyst_headless.R --fpa_dir=${FPA} --input=combined_protein_f.tsv"

## NPKL_SRRM2_Ab_FANCI vs NUCL_NPM1_AB_FANCI -> DE_results.csv
$RUN --annotation=experiment_annotation.tsv --imputation=man
