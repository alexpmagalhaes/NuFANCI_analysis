#!/usr/bin/env bash
# Step 03b: exon and intron read counts of SAD genes vs all other genes (NPKL RNA-seq).
#
# SAF files (made by NuFANCI_DNAseq steps 07b-07c):
#   NPKL_DNA_FANCI_SADs_genes.ftl.{exon,intron}.saf       genes inside SADs
#   NPKL_DNA_FANCI_All_otherGenes.ftl.{exon,intron}.saf   all other genes
# Primary alignments only, reverse-stranded (-s 2), paired-end.
#
# Output: <saf name>.tsv (+ .summary) -> summaries in results/03_gene_quantification/
#   Used in NuFANCI_DNAseq step 07 (exon/intron SAD gene plots). Tables on GEO.
#
# Original: NPKL_RNA_20250721/run_count_reads.sh. Fix: the last command was
# missing a space before "-o" ("...intron.saf-o"), so featureCounts failed on it.

mkdir -p counts && cd counts

for set in NPKL_DNA_FANCI_All_otherGenes NPKL_DNA_FANCI_SADs_genes; do
  for feature in exon intron; do
    featureCounts -F SAF --primary -s 2 -p -T 40 \
      -a ./../SAF/${set}.ftl.${feature}.saf \
      -o ${set}.ftl.${feature}.tsv \
      ./../BAM/*.bam
  done
done
