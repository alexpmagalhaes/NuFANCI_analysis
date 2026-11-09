# NuFANCI DNA-seq

> Part of the manuscript **"Dual and latent condensate localization signatures encoded in intrinsically disordered protein regions"**.
> Author and contact: **Alexandre P. Magalhaes**, Email: magalhae@molgen.mpg.de

Identification of nucleolus-associated domains (NADs) and speckle-associated domains (SADs) from NuFANCI-sorted DNA-seq, compared with matched input. All data are aligned to **T2T-CHM13 v2.0**.

Datasets:

| ID | Sample | Compartment | Replicates |
|----|--------|-------------|------------|
| HCT116_NUCL | HCT116, NPM1 Ab sort / GFP-NPM1 Ab sort / GFP-NPM1 GFP sort | nucleoli | 2 each |
| HCT116_NPKL | HCT116, speckle sort (`NPKL_DNA_FANCI`) | nuclear speckles | 4 (rep3 excluded) |
| HepG2_NUCL | HepG2, NPM1 Ab sort | nucleoli | 4 |

Raw reads, BAM, bedGraph, bigWig and domain BED files are on GEO (GSE_TBD).

## Pipeline

**Step 00, pre-processing (in-house scripts, not included):** reads were aligned to T2T-CHM13 v2.0 with **BWA-MEM**. Alignments were sorted, duplicate-marked and indexed with **Samtools**. Coverage tracks were generated from the BAM files with **deepTools** `bamCoverage`: binned bedGraphs (`*.w50k_b1k.bedgraph`), which are the input to step 01, and bigWig files for genome browsers. BAM, bedGraph and bigWig files are available on GEO.

| Step | Script(s) | Input | Output | Results in repo |
|------|-----------|-------|--------|-----------------|
| 00 | in-house (BWA-MEM, Samtools, deepTools) | FASTQ | BAM; bedGraph and bigWig coverage tracks | – |
| 01 | `01_domain_calling/01_run_FPCaller.sh` → [FPCaller](https://github.com/alexpmagalhaes/FPCaller) | NuFANCI and input bedGraphs | per-replicate domain BED, stats, chr14 QC plot | `results/01_domain_calling/` |
| 02a | `02_consensus_filtering/02a_consensus_replicates.sh` | replicate BEDs | consensus `*_merged.bed` | – |
| 02b | `02_consensus_filtering/02b_run_n_stats.sh` → `n_stats.py` | consensus BEDs | N5–N95 domain sizes | `results/02_consensus_filtering/` |
| 02c | `02_consensus_filtering/02c_filter_by_N95.sh` | consensus BEDs + N95 | `*.flt.bed` final domains (≥ N95, no centromere overlap) | – (GEO) |
| 03 | `03_signal_tracks/03_merge_beds_make_bigwigs.sh` | replicate BEDs/bedGraphs | merged BEDs, per-replicate and mean bigWigs | – (GEO) |
| 04 | `04_domain_enrichment/04_run_coverage_enrichment.sh` → `run_enrichment_NADs.py` | final domains + BAMs (deepTools `multiBamSummary`) | read counts per domain; fold change, Poisson p-value, BH FDR | `results/04_domain_enrichment/` |
| 05 | Intervene (see `05_domain_overlaps/README.md`) | final domains, HCT116 Hi-C B compartment | Venn/Euler/Jaccard overlaps | `results/05_domain_overlaps/` |
| 06 | `06_rDNA_rRNA_quantification/06_count_rDNA_DNAseq.sh` | DNA-seq BAMs, CHM13 rRNA gene GTF | rDNA read counts per rRNA gene copy | `results/06_rDNA_rRNA_quantification/` |
| 07a | `07_SAD_NAD_genes/07a_GFF_to_CSV_SAF.ipynb` | CHM13 v2.0 GFF3 | gene table / SAF | – |
| 07b | `07_SAD_NAD_genes/07b_intersect_genes_SADs_NADs.sh` | genes, SADs, NADs | SAD, NAD and other gene sets | – |
| 07c | `07_SAD_NAD_genes/07c_GFF_to_SAF_exon_intron.R` | gene-set GFFs | exon/intron SAF files | – |
| 07d | Salmon spliced/unspliced quantification, see [`NuFANCI_RNAseq`](../NuFANCI_RNAseq) step 03c | NPKL RNA-seq FASTQ | transcript quantification | – |
| 07e | `07_SAD_NAD_genes/07e_RankingSADGenes.Rmd` | gene counts, `Gene_cov_analysis_masterDF.csv` | TPM, gene rankings, density plots | `results/07_SAD_NAD_genes/ranking/` |
| 07f | `07_SAD_NAD_genes/07f_densityPlots.Rmd` | `density_graph_ratio.csv` | SAD vs non-SAD ratio density | `results/07_SAD_NAD_genes/density/` |

`scripts/legacy/` keeps earlier versions that are not used for the final results: `calculate_coverage_enrichment.py`, which step 04 replaced, and the first draft of the ranking notebook.

### Domain calling: FPCaller

Domains are called with **FPCaller**, released as a separate package:
**https://github.com/alexpmagalhaes/FPCaller** (v1.0.0, GPL-3.0).

```bash
pip install "git+https://github.com/alexpmagalhaes/FPCaller.git@v1.0.0"
fpcaller --fancy NuFANCI.w50k_b1k.bedgraph --input INPUT.w50k_b1k.bedgraph --sample NAME
```

Method in brief:
1. Keep autosomes only (`chr\d+`).
2. log2 CPM ratio NuFANCI/input per chromosome (pseudocount 1).
3. Baseline correction per chromosome: subtract a 3rd-degree polynomial fit, then take a rolling mean over 5 bins.
4. Z-score across the genome.
5. Domains = contiguous bins with z ≥ 0 that are longer than 2 kb.
6. Remove domains that overlap CHM13 centromeres (also in `references/chm13_centromeres.bed`).

RNA-seq quantification used in step 07 (gene counts, SAD exon/intron counts, Salmon) is in [`NuFANCI_RNAseq`](../NuFANCI_RNAseq).

## References

- `references/chm13_chrm.sizes`: chromosome sizes used for bigWigs and `bedtools complement`.
- `references/rRNA_loci.chm13.SAF.tsv`: rDNA arrays on chr13/14/15/21/22.
- `references/chm13_centromeres.bed`: CHM13 centromere intervals (same as in FPCaller), used in step 02c.
- External (not included): T2T-CHM13 v2.0 genome and annotation (`chm13.draft_v2.0.gene_annotation.gff3`); human hs1-rDNA reference v1.0 (`HS1_chrR`, `hs1-rDNA_v1.0.SAF`).

## Software

[FPCaller](https://github.com/alexpmagalhaes/FPCaller) v1.0.0; Python 3 (pandas, numpy, scipy, statsmodels, seaborn, matplotlib, pybedtools, gffpandas); R (tidyverse, Rsubread, ggprism, ggpubr, DGEobj.utils, introdataviz); BWA-MEM, Samtools, bedtools, deepTools, UCSC `bedGraphToBigWig`, wiggletools, featureCounts (Subread), Intervene. Versions: see `envs/` (TBD).

## Large tables not in the repo

These are suggested as GEO supplementary files:

- `Gene_cov_analysis_masterDF.csv` / `_filtered.csv`
- `NPKL_RNA_quant_TPM.csv`
- `chm13.draft_v2.0.gene.counts.tsv`
- exon/intron count tables
- Salmon `quant`
- `NUC_NPM1_Ab_sort_HepG2_merged.tsv`
