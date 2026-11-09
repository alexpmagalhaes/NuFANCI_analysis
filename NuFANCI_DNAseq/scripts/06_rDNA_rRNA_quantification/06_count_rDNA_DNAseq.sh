#!/usr/bin/env bash
# Step 06: rDNA content of NuFANCI DNA-seq vs input.
#
# Reads from the BWA CHM13 alignments (duplicate-marked) are counted on the CHM13 v2.0
# rRNA gene annotation (chm13.draft_v2.0.rRNA.fix.gtf; 45S/28S/18S/5.8S/5S copies on the
# acrocentric arms and the 5S array). Multi-mapping reads are counted (-M) because rDNA
# copies are near-identical; duplicates are ignored; only properly paired fragments (-B).
#
# Output: chm13.draft_v2.0.rrna_gene_counts.tsv (+ .summary)
#         -> results/06_rDNA_rRNA_quantification/total/
#
# BAM names: GFPcells = HCT116 GFP-NPM1, WTcells = HCT116 WT (NPM1 Ab sort), NPKLS = speckle sort.
# Command recovered from the featureCounts header of the output table (featureCounts v2.0.3).
# RNA-seq rRNA quantification is in NuFANCI_RNAseq step 02.

featureCounts -T 40 -B -M --ignoreDup -p --countReadPairs \
  -t gene -g gene_id \
  -a chm13.draft_v2.0.rRNA.fix.gtf \
  -o chm13.draft_v2.0.rrna_gene_counts.tsv \
  ./bam/GFPcells_Ab_sort_rep1_CHM13.bwa.srt.md.bam \
  ./bam/GFPcells_Ab_sort_rep1_input_CHM13.bwa.srt.md.bam \
  ./bam/GFPcells_Ab_sort_rep2_CHM13.bwa.srt.md.bam \
  ./bam/GFPcells_Ab_sort_rep2_input_CHM13.bwa.srt.md.bam \
  ./bam/GFPcells_GFP_sort_rep1_CHM13.bwa.srt.md.bam \
  ./bam/GFPcells_GFP_sort_rep1_input_CHM13.bwa.srt.md.bam \
  ./bam/GFPcells_GFP_sort_rep2_CHM13.bwa.srt.md.bam \
  ./bam/GFPcells_GFP_sort_rep2_input_CHM13.bwa.srt.md.bam \
  ./bam/NPKLS_Ab_sort_rep1_c50_CHM13.bwa.srt.md.bam \
  ./bam/NPKLS_Ab_sort_rep1_c50_input_CHM13.bwa.srt.md.bam \
  ./bam/WTcells_Ab_sort_rep1_CHM13.bwa.srt.md.bam \
  ./bam/WTcells_Ab_sort_rep1_input_CHM13.bwa.srt.md.bam \
  ./bam/WTcells_Ab_sort_rep2_CHM13.bwa.srt.md.bam \
  ./bam/WTcells_Ab_sort_rep2_input_CHM13.bwa.srt.md.bam
