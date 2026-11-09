# NuFANCI RNA-seq: quantification

> Part of the manuscript **"Dual and latent condensate localization signatures encoded in intrinsically disordered protein regions"**.
> Author and contact: **Alexandre P. Magalhaes**, Email: magalhae@molgen.mpg.de

Quantification of RNA from NuFANCI-sorted nucleoli (NUCL) and nuclear speckles (NPKL) compared with matched input, in HCT116 cells. It covers two analyses:

- **rRNA quantification:** reads on the human rDNA unit (hs1-rDNA reference)
- **gene-level quantification:** CHM13 genes, SAD vs other genes, spliced/unspliced

Raw and processed data (FASTQ, BAM, bigWig, full count tables) are on GEO ([GSE335756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335756)).

## Datasets

| Dataset | Samples | Replicates | Library |
|---------|---------|------------|---------|
| NUCL_RNA | nucleolar NuFANCI sort + input | 5 + 5 | total RNA-seq, PE, reverse-stranded |
| NPKL_RNA | speckle NuFANCI sort + input | 5 + 5 | total RNA-seq, PE, reverse-stranded |
| NPKL_SMARTseq_nbGFP | speckle NuFANCI + input, Dox / noDox | 3 each | SMART-seq, PE, UMI |
| NPKL_SMARTseq_ActD | speckle NuFANCI + input, DMSO / Actinomycin D | 4 each | SMART-seq, PE, UMI |

The full sample list is in `metadata/samples.tsv`.

## Pipeline

**Step 00, pre-processing (in-house scripts, not included):**
- reads were trimmed with Cutadapt and checked with FastQC/MultiQC;
- SMART-seq reads were UMI-extracted and deduplicated with UMI-tools;
- reads were aligned to T2T-CHM13 v2.0 with **STAR**, and to the rDNA reference with **STAR** and **bowtie2** (step 02);
- **Samtools** was used for sorting, filtering and indexing;
- **deepTools** `bamCoverage` was used to make bigWig coverage tracks.

| Step | Script | Input | Output | Results in repo |
|------|--------|-------|--------|-----------------|
| 01 | `01_rDNA_reference/01_STAR_index_hs1_rDNA.sh` | hs1-rDNA FASTA + GTF | STAR rDNA index | – |
| 02a | `02_rRNA_quantification/02a_rRNA_bowtie2_featureCounts_NUCL.sh` | NUCL trimmed FASTQ | bowtie2 rDNA BAMs; read pairs per rDNA feature (featureCounts) | `results/02_rRNA_quantification/NUCL_bowtie2/` |
| 02b | `02_rRNA_quantification/02b_rRNA_STAR_NUCL.sh` | NUCL trimmed FASTQ (reps 1, 4, 5) | STAR rRNA gene counts, splice junctions | `.../NUCL_STAR/` |
| 02c | `02_rRNA_quantification/02c_rRNA_STAR_NPKL.sh` | NPKL trimmed FASTQ (reps 2, 4, 5) | STAR rRNA gene counts | `.../NPKL_STAR/` |
| 02d | `02_rRNA_quantification/02d_rRNA_STAR_SMARTseq_nbGFP.sh` | SMART-seq deduplicated FASTQ | STAR rRNA gene counts | `.../SMARTseq_nbGFP_STAR/` |
| 02e | `02_rRNA_quantification/02e_rDNA_strand_coverage_bigwig.sh` | STAR rDNA BAMs (02b) | strand-specific BPM/CPM bigWigs | – (GEO); junctions and sashimi in `.../rDNA_splice_junctions/` |
| 02f | `02_rRNA_quantification/02f_rRNA_STAR_SMARTseq_ActD.sh` (reconstructed) | SMART-seq ActD/DMSO deduplicated FASTQ | STAR rRNA gene counts; `ActD_rRNA_reads.csv` (% reads on rDNA) | `.../ActD_DMSO_STAR/` |
| 03a | `03_gene_quantification/03a_featureCounts_genes_NPKL.sh` | NPKL STAR CHM13 BAMs | gene counts | `results/03_gene_quantification/` (summary) |
| 03b | `03_gene_quantification/03b_featureCounts_SAD_genes_exon_intron_NPKL.sh` | NPKL STAR CHM13 BAMs, SAD/other gene SAFs | exon and intron counts | `results/03_gene_quantification/` (summaries) |
| 03c | `03_gene_quantification/03c_salmon_spliced_unspliced_NPKL.sh` | NPKL trimmed FASTQ | Salmon spliced/unspliced quantification | – (GEO) |

**Reading STAR rRNA counts:** `*ReadsPerGene.out.tab` has 4 columns: gene, unstranded, forward, reverse. The libraries are reverse-stranded, so **use column 4**.

**ActD / DMSO (step 02f):** `ActD_rRNA_reads.csv` gives, per sample, the number of input reads, the reads mapping uniquely to the rDNA, and the % rRNA (= uniquely mapped / input × 100), taken from STAR `Log.final.out`. The original script for this run was not kept. 02f was rebuilt from the other SMART-seq run (02d), and its summary reproduces all 16 rows of the CSV exactly.

### Links to the DNA-seq module

- The SAF files used in 03a/03b (CHM13 genes; SAD and other genes, exon/intron) are made in [`NuFANCI_DNAseq`](../NuFANCI_DNAseq) steps 07a–07c.
- The gene counts and Salmon output are used for the SAD gene ranking in `NuFANCI_DNAseq` step 07e.
- The rDNA content of the DNA-seq samples is quantified in `NuFANCI_DNAseq` step 06.

## References

- `references/hs1-rDNA_v1.0.gtf`: one gene `rRNA` over the 45S unit, used for STAR GeneCounts.
- `references/hs1-rDNA_v1.0.SAF`: rDNA features (IGS, promoter, 5'ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3'ETS ...), used for featureCounts (02a).
- External: human hs1-rDNA genome v1.0 (`HS1_chrR.fa`), T2T-CHM13 v2.0 genome and gene annotation.

## Software

See `envs/software_versions.md`.

## Large tables not in the repo

These are on GEO:

- per-sample STAR `ReadsPerGene.out.tab` for CHM13
- `chm13.draft_v2.0.gene.counts.tsv`
- SAD/other gene exon and intron count tables
- Salmon `quant.sf`
- STAR rDNA splice-junction tables (`*SJ.out.tab`)
