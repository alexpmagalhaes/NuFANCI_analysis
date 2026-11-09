#!/usr/bin/env bash
# Step 02c: rRNA quantification of nuclear-speckle (NPKL) RNA-seq with STAR (replicates 2, 4, 5).
#
# Same procedure as step 02b. Libraries are reverse-stranded: use column 4 of
# *ReadsPerGene.out.tab.
#
# Output: NPKL_RNA_{FANCI,Input}_rep{2,4,5}.HS1_chrR.strReadsPerGene.out.tab
#         -> results/02_rRNA_quantification/NPKL_STAR/
#
# Fix vs. the original script: samtools sort read "NUC_RNA_FANCI_rep*" files,
# but STAR writes "NPKL_RNA_*" here.

GENOME=/project/hnisz_seq_data/FANCI_seq/RNA_seq/NUCL_RNA_20250721/rRNA_workflow_bowtie/genome/STAR

for sample in NPKL_RNA_FANCI NPKL_RNA_Input; do
  for i in 2 4 5; do

    STAR --runThreadN 40 \
      --runMode alignReads \
      --genomeDir ${GENOME} \
      --sjdbGTFfile ${GENOME}/hs1-rDNA_v1.0.gtf \
      --readFilesIn ./fastq/${sample}_rep${i}_R1.trm.fq.gz ./fastq/${sample}_rep${i}_R2.trm.fq.gz \
      --outFileNamePrefix ${sample}_rep${i}.HS1_chrR.str \
      --outStd Log \
      --readFilesCommand gzip -dc \
      --quantMode GeneCounts

    samtools sort -@ 40 ${sample}_rep${i}.HS1_chrR.strAligned.out.sam \
      -o ${sample}_rep${i}.HS1_chrR.srt.bam
  done
done
