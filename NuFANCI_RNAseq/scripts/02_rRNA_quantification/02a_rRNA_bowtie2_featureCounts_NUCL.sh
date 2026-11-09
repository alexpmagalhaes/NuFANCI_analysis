#!/usr/bin/env bash
# Step 02a: rRNA quantification of nucleolar (NUCL) NuFANCI and input RNA-seq with bowtie2.
#
# Trimmed reads are aligned to the hs1-rDNA reference (chrR) with bowtie2; mapped
# reads with MAPQ >= 1 are kept, sorted and indexed, then read pairs are counted per
# rDNA feature (5'ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3'ETS, IGS ...) with featureCounts.
#
# Output: NUC_RNA_FANCI_input (+ .summary) -> results/02_rRNA_quantification/NUCL_bowtie2/
#
# Fixes vs. the original script: the loops used {1..$END}, which bash does not
# expand (loop ran once with i="{1..5}"); featureCounts was missing "-F SAF".

WRKDIR=/project/hnisz_seq_data/FANCI_seq/RNA_seq/NUCL_RNA_20250721/rRNA_workflow_bowtie
SAF=${WRKDIR}/genome/Human_hs1-rDNA_genome_v1.0_annotation/hs1-rDNA_v1.0.SAF

END=5

for sample in NUC_RNA_FANCI NUC_RNA_Input; do
  for i in $(seq 1 $END); do
    prefix=${sample}_rep${i}.HS1_chrR.bwt2

    # Align (max fragment length 2 kb)
    bowtie2 -X 2000 -p 20 \
      -x ${WRKDIR}/genome/HS1_chrR \
      -1 ${sample}_rep${i}_R1.trm.fq.gz \
      -2 ${sample}_rep${i}_R2.trm.fq.gz \
      -S ${prefix}.sam

    # Keep mapped reads with MAPQ >= 1, sort, index
    samtools view -h -F 4 -q 1 -bS ${prefix}.sam > ${prefix}.bam
    samtools sort -@ 20 ${prefix}.bam -o ${prefix}.srt.bam
    samtools index -@ 20 ${prefix}.srt.bam
  done
done

## Count read pairs per rDNA feature
featureCounts -T 20 \
  -F SAF \
  -a ${SAF} \
  -o NUC_RNA_FANCI_input \
  -p --countReadPairs *.HS1_chrR.bwt2.srt.bam
