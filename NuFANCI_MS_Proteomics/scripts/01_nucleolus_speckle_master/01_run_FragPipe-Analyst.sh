#!/usr/bin/env bash
# Step 01: differential abundance with FragPipe-Analyst, run from the command line
# (no Shiny app) via scripts/00_common/run_FragPipe-Analyst_headless.R.
#
# Settings = those of the original GUI run (see 01_FragPipe-Analyst_log_20260520.txt):
#   LFQ, MaxLFQ intensities, contaminants removed, no missing-value filter,
#   VSN normalisation, min imputation, limma all pairwise contrasts, BH,
#   significant if p.adj <= 0.05 and |log2FC| >= 1.
# (The saved GUI log is from an earlier run with 47 samples; the published results
#  use the 38 samples in experiment_annotation.tsv. Run A reproduces
#  DE_results_minImputation.csv exactly.)
#
# Run in the folder with combined_protein_f.tsv (step 00) and the annotation files
# from metadata/01_nucleolus_speckle_master/.
#   REPO = path to the NuFANCI_MS_Proteomics folder of this repository
#   FPA  = path to FragPipe-Analyst (git clone, commit 3e1af92)

REPO=${REPO:-../NuFANCI_MS_Proteomics}
FPA=${FPA:-../FragPipe-Analyst}
RUN="Rscript ${REPO}/scripts/00_common/run_FragPipe-Analyst_headless.R --fpa_dir=${FPA} --input=combined_protein_f.tsv"

## A) all samples (NuFANCI and input) -> DE_results_minImputation.csv, dual comparisons
$RUN --annotation=experiment_annotation.tsv --suffix=_minImputation \
  --dual=NUCL_HCT116_NPM1_FANCI_vs_NUCL_HCT116_NPM1_Input,NUCL_HCT116_SRRM2_FANCI_vs_NUCL_HCT116_SRRM2_Input \
  --dual=NUCL_HCT116_NPM1_FANCI_vs_NUCL_HCT116_NPM1_Input,NUCL_HSLiver_NPM1_FANCI_vs_NUCL_HSLiver_NPM1_Input \
  --dual=NUCL_HCT116_NPM1_FANCI_vs_NUCL_HCT116_NPM1_Input,NUCL_HepG2_NPM1_FANCI_vs_NUCL_HepG2_NPM1_Input \
  --dual=NUCL_HepG2_NPM1_FANCI_vs_NUCL_HepG2_NPM1_Input,NUCL_HSLiver_NPM1_FANCI_vs_NUCL_HSLiver_NPM1_Input

## B) NuFANCI samples only -> DE_results_onlyFANCI_minImputation.csv
$RUN --annotation=experiment_annotation_FANCI.tsv --suffix=_onlyFANCI_minImputation
