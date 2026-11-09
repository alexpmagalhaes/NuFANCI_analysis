# NuFANCI mass-spectrometry proteomics

> Part of the manuscript **"Dual and latent condensate localization signatures encoded in intrinsically disordered protein regions"**.
> Author and contact: **Alexandre P. Magalhaes**, Email: magalhae@molgen.mpg.de

Label-free quantitative proteomics of NuFANCI-sorted nuclear bodies compared with matched input.

Raw files and protein/peptide tables (FragPipe `combined_protein.tsv`, filtered, normalised and imputed matrices) are deposited in PRIDE ([PXD069855](https://www.ebi.ac.uk/pride/archive/projects/PXD069855)) and are not stored here.

## Experiments

| # | Folder | Content |
|---|--------|---------|
| 01 | `01_nucleolus_speckle_master` | Nucleolar NuFANCI (NPM1 sort) in HCT116, HepG2 and human liver; speckle NuFANCI (SRRM2 sort) in HCT116; NUP98-DDX10 NuFANCI in HCT116; each with matched input |
| 02 | `02_nucleolus_vs_speckle` | Nucleolar (NPM1 antibody sort) vs nuclear speckle (SRRM2 antibody sort) NuFANCI in HCT116 |
| 03 | `03_speckle_ActD` | Nuclear speckle NuFANCI and input from HCT116 treated with DMSO or Actinomycin D (ActD) |

## Pre-processing

Raw data (Bruker `.d`, DDA) were searched with **FragPipe** (not included) and quantified with MaxLFQ (`combined_protein.tsv`).

**Step 00, filtering** (`scripts/00_common/filter_proteins.py`, run by `00_filter_proteins.sh` in each experiment): removes from `combined_protein.tsv`
- proteins on the contaminant list (`reference_lists/Contaminants.tsv`; Frankenfield et al., *J Proteome Res* 2022, doi:10.1021/acs.jproteome.2c00145);
- contaminant entries from the search FASTA (`contam_sp…` / `Cont_`);
- proteins with ≤ 1 total peptide.

The output is `combined_protein_f.tsv`.

**Step 01, differential abundance** with **FragPipe-Analyst** ([MonashProteomics/FragPipe-Analyst](https://github.com/MonashProteomics/FragPipe-Analyst), commit `3e1af92`), LFQ mode on MaxLFQ intensities:
- VSN normalisation;
- imputation: minimum value (experiment 01) or Perseus-type (experiments 02, 03);
- limma, Benjamini–Hochberg correction;
- significant if adjusted p ≤ 0.05 and |log2 FC| ≥ 1.

The analyses were run in a local FragPipe-Analyst installation (Shiny app in RStudio). To run them without the app, use `scripts/00_common/run_FragPipe-Analyst_headless.R`. It calls FragPipe-Analyst's own R functions and writes the same files as the app's downloads: DE results, original/filtered/normalised/imputed matrices, dual comparisons and a run log. It uses the local installation (default `~/FragPipe-Analyst`), including its renv package library, so the package versions are those of the app. Each experiment's `01_run_FragPipe-Analyst.sh` holds the exact settings.

```bash
Rscript NuFANCI_MS_Proteomics/scripts/00_common/run_FragPipe-Analyst_headless.R --help     # all options

Rscript NuFANCI_MS_Proteomics/scripts/00_common/run_FragPipe-Analyst_headless.R \
  --input=combined_protein_f.tsv --annotation=experiment_annotation.tsv \
  --lfq_type=MaxLFQ --normalization=vsn --imputation=min --outdir=out --plots=TRUE
```

Main options (defaults are the FragPipe-Analyst app defaults):

| Option | Values (default first) |
|--------|------------------------|
| `--data_type` | `LFQ` (combined_protein.tsv), `LFQ-peptide` (combined_peptide.tsv), `DIA` (FragPipe abundance_protein / DIA-NN pg_matrix), `TMT` (TMT-Integrator report) |
| `--lfq_type` | `Intensity`, `MaxLFQ`, `Spectral Count` |
| `--remove_contaminants` | `TRUE`, `FALSE` |
| `--min_global`, `--min_condition` | % of non-missing values required (`0` = no filter) |
| `--normalization` | `none`, `vsn`, `MD` (median-centred) |
| `--imputation` | `man` (Perseus-type; `none` for TMT), `none`, `knn`, `MLE`, `min`, `zero`, `MinProb`, `MinDet`, `QRILC`, `bpca`, `RF`; `--man_shift`, `--man_scale` |
| `--fdr` | `BH` (limma), `fdrtool` |
| `--de_type` | `all` pairs, `control` (with `--control=A;B`), `others` |
| `--p`, `--lfc` | `0.05`, `1` |
| `--dual` | `CONTRAST1,CONTRAST2` (repeatable); `--dual_alpha`, `--dual_lfc`, `--dual_p_adj` |
| `--plots` | `FALSE`, `TRUE`: PCA, correlation, heatmap, volcano per contrast, missing values, density, protein numbers, coverage, CVs |
| `--matrices`, `--outdir`, `--suffix` | output control |
| `--fpa_dir`, `--use_renv` | FragPipe-Analyst installation (default `~/FragPipe-Analyst`; its renv library is used) |

The script checks options and file names before it runs. Some imputation methods need extra R packages (`knn`: impute, `bpca`: pcaMethods, `RF`: missForest, `MLE`: norm, `MinProb`/`MinDet`/`QRILC`: imputeLCMD), and `fdrtool` needs the fdrtool package; the script stops with a clear message if one is missing. If FragPipe-Analyst is not installed: `git clone https://github.com/MonashProteomics/FragPipe-Analyst.git ~/FragPipe-Analyst && cd ~/FragPipe-Analyst && git checkout 3e1af92`, then restore its renv library.

**Verification:** run through each experiment's `01_run_FragPipe-Analyst.sh`, the script reproduced the app's `DE_results` for all three experiments, including the NuFANCI-only run of experiment 01. The proteins, significance calls and imputed matrices are identical, and numeric differences are ≤ 1e-13.

## 01: nucleolus / speckle master experiment

Sample annotation: `metadata/01_nucleolus_speckle_master/experiment_annotation.tsv` (all samples) and `experiment_annotation_FANCI.tsv` (NuFANCI samples only).

| Step | Script | Output (in `results/01_nucleolus_speckle_master/`) |
|------|--------|---------------------------------------------------|
| 00 | `00_filter_proteins.sh` → `scripts/00_common/filter_proteins.py`: contaminant and single-peptide filtering | `combined_protein_f.tsv` (PRIDE [PXD069855](https://www.ebi.ac.uk/pride/archive/projects/PXD069855)) |
| 01 | `01_run_FragPipe-Analyst.sh` (VSN, min imputation, all pairwise contrasts; A: all samples, B: NuFANCI only). The GUI log is in `01_FragPipe-Analyst_log_20260520.txt` | `differential_abundance/DE_results_minImputation.csv`, `DE_results_onlyFANCI_minImputation.csv`; `co_enrichment/dual_comparison_*.tsv` |
| 02 | `02_clustering_heatmap_FANCI_and_Input.ipynb`: clustering of all samples, annotated with published marker lists | `clustering/Heatmap_clustered_Input_FANCI_centered.*`, `den_Heatmap_*FANCI_Input*`, `FANCI_Input_ordered_row_df_heatmap_wClustering.csv` |
| 03 | `03_clustering_heatmap_FANCI_only.ipynb`: same, NuFANCI samples only | `clustering/Heatmap_clustered_FANCI_centered.*`, `den_Heatmap_*FANCI_log*`, `FANCI_ordered_row_df_heatmap_wClustering.csv` |
| 04 | `04_co_enrichment_NPM1_vs_SRRM2.Rmd`: nucleolar vs speckle log2 FC scatter | `co_enrichment/protein_srrm2_vs_nucleoli.pdf`, `dual_comparison_scatter_*.pdf` |
| 05 | `05_VolcanoPlot.Rmd`: volcano plots | – |
| 06 | `06_GO_enrichment.Rmd`: GO enrichment (Enrichr, via FragPipe-Analyst functions) | – |
| 07 | `07_Gene_plots.Rmd`: intensities of selected proteins | – |

Steps 05–07 come from the earlier analysis of the nucleolar samples (2025-10-19).

Run the notebooks from the folder that holds the FragPipe / FragPipe-Analyst tables from PRIDE ([PXD069855](https://www.ebi.ac.uk/pride/archive/projects/PXD069855)), with the module's `reference_lists/` folder copied or linked into it (`ln -s <repo>/NuFANCI_MS_Proteomics/reference_lists .`). In the R notebooks, set the data folder with the `setwd()` line at the top.

## 02: nucleolus vs speckle

Sample annotation: `metadata/02_nucleolus_vs_speckle/experiment_annotation.tsv` (NPM1 and SRRM2 antibody NuFANCI sorts, 3 replicates each).

| Step | Script | Output (in `results/02_nucleolus_vs_speckle/`) |
|------|--------|-----------------------------------------------|
| 00 | `00_filter_proteins.sh` | `combined_protein_f.tsv` (PRIDE [PXD069855](https://www.ebi.ac.uk/pride/archive/projects/PXD069855)) |
| 01 | `01_run_FragPipe-Analyst.sh` (VSN, Perseus-type imputation) | `differential_abundance/DE_results.csv`, `FragPipe-Analyst_report.pdf` (original app report) |
| 02 | `02_analysis_nucleolus_vs_speckle.ipynb`: detection, Euler diagram, correlation, PCA, clustered heatmap with marker lists, intensity ranking, volcano plot | `figures/`, `tables/NPKL_SRRM2_Ab_FANCI.csv`, `tables/NUCL_NPM1_AB_FANCI.csv`, `differential_abundance/Volcano_NPKL_vs_NUCL_Data.csv` |
| 03 | `03_select_proteins_for_GO.ipynb`: protein selection from Supplementary Table S3 | `GO_enrichment/filtered_results_GreyDots_plus5values.csv` |
| 04 | `04_GO_enrichment.Rmd`: GO enrichment of nucleolus-, speckle-enriched and selected proteins | `GO_enrichment/GO*.pdf`, `GO_enrichment_GO_Biological_Process_2023.csv` |

## 03: nuclear speckle NuFANCI after Actinomycin D

Sample annotation: `metadata/03_speckle_ActD/experiment_annotation.tsv` (DMSO / ActD × NuFANCI / input).

| Step | Script | Output (in `results/03_speckle_ActD/`) |
|------|--------|---------------------------------------|
| 00 | `00_filter_proteins.sh` | `combined_protein_f.tsv` (PRIDE [PXD069855](https://www.ebi.ac.uk/pride/archive/projects/PXD069855)) |
| 01 | `01_run_FragPipe-Analyst.sh` (VSN, Perseus-type imputation, all pairwise contrasts) | `differential_abundance/DE_results.csv` |
| 02 | `02_analysis_speckle_ActD.ipynb`: detection, correlation, PCA, clustered heatmap with marker lists | `figures/`, `tables/NSPKL_*_FANCI.csv` (detected proteins) |
| 03 | `03_Gene_plots.Rmd`: intensities of selected speckle proteins per condition | `figures/Marker_gene_ActD.pdf` |

The ActD notebook was adapted from experiment 02. Its later cells (intensity ranking, volcano plot) still use experiment 02's sample names.

## Reference lists (`reference_lists/`)

| File | Content / source |
|------|------------------|
| `Contaminants.tsv` | UniProt IDs of common contaminants removed in step 00 (from the universal contaminant libraries of Frankenfield et al., *J Proteome Res* 2022; 21(9): 2104–2113, [doi:10.1021/acs.jproteome.2c00145](https://doi.org/10.1021/acs.jproteome.2c00145)) |
| `Hsapiens_Uniprot_annotation.tsv` | Human UniProt annotation (Protein, Entry Name, Gene, Description) |
| `NucleolusList/Andersen2005*.{txt,tsv}` | Nucleolar proteome, Andersen et al. 2005 (with UniProt ID mapping) |
| `NucleolusList/SubCell_markerProteins*.csv` | SubCellBarcode marker proteins (Orre et al. 2019) |
| `NucleolusList/msb209469-sup-0002-datasetev1.xlsx` | Stenström et al. 2020, Mol Syst Biol, Dataset EV1 |
| `NucleolusList/Zang2025.csv` | Nucleolar proteins, Zhang et al. 2025 |
| `Speackle_lists_media-1.xlsx` | Nuclear speckle and paraspeckle proteins (Dyakov et al. 2024) |
| `Dopie_etall_2020_jcb_201910207_tables1.xlsx` | Dopie et al. 2020, J Cell Biol, Table S1 (top 400 speckle proteins) |

## Software

See `envs/software_versions.md`.
