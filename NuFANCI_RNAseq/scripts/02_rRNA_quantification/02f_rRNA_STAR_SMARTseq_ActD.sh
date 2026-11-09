#!/usr/bin/env bash
# Step 02f: rRNA quantification of the SMART-seq Actinomycin D experiment
# (HCT116, DMSO vs ActD, speckle NuFANCI sort and input, 4 replicates each) with STAR,
# followed by a per-sample summary of the fraction of reads on the rDNA.
#
# NOTE: reconstructed. The original script for this run was not kept. It follows the
# same STAR settings as the other SMART-seq run (step 02d). The summary part
# reproduces ActD_rRNA_reads.csv exactly (checked against the STAR count tables).
#
# Input reads are the UMI-extracted, trimmed and deduplicated FASTQs from the
# in-house SMART-seq pipeline (*.extr.trm.dedup.fq).
#
# Output:
#   <sample>.HS1_chrR.strReadsPerGene.out.tab, <sample>.HS1_chrR.strLog.final.out
#   ActD_rRNA_reads.csv   Sample, Input_Reads, Aligned_Reads, % rRNA
#                         Input_Reads   = "Number of input reads"
#                         Aligned_Reads = "Uniquely mapped reads number" (on the rDNA)
#                         % rRNA        = Aligned_Reads / Input_Reads * 100
#   -> results/02_rRNA_quantification/ActD_DMSO_STAR/

GENOME=/project/hnisz_seq_data/FANCI_seq/RNA_seq/NUCL_RNA_20250721/rRNA_workflow_bowtie/genome/STAR

SAMPLES="HCT116_DMSO_Input HCT116_DMSO_FANCI HCT116_ActD_Input HCT116_ActD_FANCI"

## Align to the hs1-rDNA reference and count reads on the rRNA gene
for sample in $SAMPLES; do
  for i in 1 2 3 4; do

    STAR --runThreadN 40 \
      --runMode alignReads \
      --genomeDir ${GENOME} \
      --sjdbGTFfile ${GENOME}/hs1-rDNA_v1.0.gtf \
      --readFilesIn ./../fastq/${sample}_rep${i}_R1.extr.trm.dedup.fq ./../fastq/${sample}_rep${i}_R2.extr.trm.dedup.fq \
      --outFileNamePrefix ${sample}_rep${i}.HS1_chrR.str \
      --outStd Log \
      --quantMode GeneCounts
  done
done

## Summarise: % of input reads that map uniquely to the rDNA
echo "Sample,Input_Reads,Aligned_Reads,% rRNA" > ActD_rRNA_reads.csv
for sample in $SAMPLES; do
  for i in 1 2 3 4; do
    log=${sample}_rep${i}.HS1_chrR.strLog.final.out
    input=$(grep "Number of input reads" $log | cut -f2)
    unique=$(grep "Uniquely mapped reads number" $log | cut -f2)
    pct=$(awk -v a=$unique -v b=$input 'BEGIN{printf "%.8f", a/b*100}')
    echo "${log},${input},${unique},${pct}" >> ActD_rRNA_reads.csv
  done
done
