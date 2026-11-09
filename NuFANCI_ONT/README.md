# NuFANCI ONT: CpG methylation of nucleolar DNA

> Part of the manuscript **"Dual and latent condensate localization signatures encoded in intrinsically disordered protein regions"**.
>
> Authors:
> - **Alexandre P. Magalhaes**, Email: magalhae@molgen.mpg.de (contact)
> - **Henri Niskanen**, Email: niskanen@molgen.mpg.de (R plotting scripts)

Oxford Nanopore (ONT) sequencing of nucleolar NuFANCI-sorted DNA and matched input from HCT116 cells, used to measure CpG methylation at the rDNA genes.

Raw (POD5), aligned (BAM) and methylation (bedMethyl / bigWig) data are on GEO (GSE_TBD).

## Pre-processing (step 00, not included)

1. **Basecalling:** POD5 data from ONT R10.4.1 flow cells were base-called with **dorado v1.0.2** (`sup` mode) using the 5mC_5hmC v5.0.0 DNA model.
2. **Trimming:** reads were converted to FASTQ (**samtools fastq**), and barcodes and adapters were removed with **porechop v0.2.4**.
3. **Back to BAM:** trimmed reads were converted back to BAM (**samtools import**), and the MM/ML modification tags were repaired with **modkit v0.5.0**.
4. **Alignment:** reads were mapped to **T2T-CHM13 v2.0** with **minimap2**.
5. **Methylation calling:** **modkit** produced bedMethyl files, keeping CpGs with at least 5× coverage.

For plotting, the CpG (+ strand) methylation values were converted to bigWig.

## Scripts

Scripts written by **Henri Niskanen** (niskanen@molgen.mpg.de).

| Script | Samples | Output |
|--------|---------|--------|
| `scripts/01_rDNA_CGmeth_heatmap_HCT116_GFPNPM1sort.R` | HCT116 GFP-NPM1 sort, NuFANCI vs input (ONT run 20250618) | CpG methylation heatmaps at the 45S rDNA genes |
| `scripts/02_rDNA_CGmeth_heatmap_HCT116_NPM1sort_rep2.R` | HCT116 NPM1 sort replicate 2, NuFANCI vs input | CpG methylation heatmaps at the 45S rDNA genes |

Each script plots methylation (0–100 %) around the five annotated 45S rDNA copies (RNA45S1-1 … RNA45S5-1, from the CAT/Liftoff CHM13 annotation `catLiftOffGenesV1.bed`) in two views:
1. gene body ± 15 kb;
2. 30 kb upstream (IGS/promoter) + gene body.

## Software

See `envs/software_versions.md`.
