#!/usr/bin/env bash
# Step 03a: gene-level read counts for nuclear-speckle (NPKL) RNA-seq.
#
# Read pairs from the STAR CHM13 alignments (in-house pipeline) are counted per gene
# with featureCounts on the CHM13 v2.0 gene SAF (genes only; made by NuFANCI_DNAseq step
# 07a). Reverse-stranded (-s 2); reads overlapping several genes are counted for each (-O).
#
# Output: chm13.draft_v2.0.gene.counts.tsv (+ .summary)
#   Used for TPM / gene ranking in NuFANCI_DNAseq step 07e. Table on GEO.
#
# Command recovered from the header of chm13.draft_v2.0.gene.counts.tsv (featureCounts v2.1.1).

mkdir -p counts && cd counts

featureCounts -F SAF -p --countReadPairs -T 40 -O -s 2 \
  -a ./../SAF/chm13.draft_v2.0.gene_annotation.onlygenes.SAF \
  -o chm13.draft_v2.0.gene.counts.tsv \
  ./../BAM/NPKL_RNA_FANCI_rep1_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_FANCI_rep2_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_FANCI_rep3_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_FANCI_rep4_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_FANCI_rep5_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_Input_rep1_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_Input_rep2_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_Input_rep3_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_Input_rep4_chm13.star.srt.bam \
  ./../BAM/NPKL_RNA_Input_rep5_chm13.star.srt.bam
