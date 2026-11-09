#!/usr/bin/env bash
# Step 02b: rRNA quantification of nucleolar (NUCL) RNA-seq with STAR (replicates 1, 4, 5).
#
# Trimmed reads are aligned to the hs1-rDNA reference (index from step 01) and reads
# per rRNA gene are counted with STAR --quantMode GeneCounts.
# Libraries are reverse-stranded: use column 4 of *ReadsPerGene.out.tab.
#
# Output: <sample>.HS1_chrR.strReadsPerGene.out.tab, <sample>.HS1_chrR.strSJ.out.tab
#         -> results/02_rRNA_quantification/NUCL_STAR/
# The sorted NuFANCI BAMs are used for strand-specific coverage tracks (step 02e).

GENOME=/project/hnisz_seq_data/FANCI_seq/RNA_seq/NUCL_RNA_20250721/rRNA_workflow_bowtie/genome/STAR

for sample in NUC_RNA_FANCI NUC_RNA_Input; do
  for i in 1 4 5; do

    STAR --runThreadN 40 \
      --runMode alignReads \
      --genomeDir ${GENOME} \
      --sjdbGTFfile ${GENOME}/hs1-rDNA_v1.0.gtf \
      --readFilesIn ./trim_fastq/${sample}_rep${i}_R1.trm.fq.gz ./trim_fastq/${sample}_rep${i}_R2.trm.fq.gz \
      --outFileNamePrefix ${sample}_rep${i}.HS1_chrR.str \
      --outStd Log \
      --readFilesCommand gzip -dc \
      --outSAMattributes XS \
      --quantMode GeneCounts

    samtools sort -@ 40 ${sample}_rep${i}.HS1_chrR.strAligned.out.sam \
      -o ${sample}_rep${i}.HS1_chrR.srt.bam
  done
done
