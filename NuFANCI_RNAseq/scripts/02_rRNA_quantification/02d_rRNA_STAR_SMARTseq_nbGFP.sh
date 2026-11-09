#!/usr/bin/env bash
# Step 02d: rRNA quantification of the SMART-seq nbGFP experiment (HCT116, Dox vs noDox,
# NuFANCI and input, 3 replicates each) with STAR.
#
# Input reads are the UMI-extracted, trimmed and deduplicated FASTQs produced by the
# in-house SMART-seq pipeline (*.extr.trm.dedup.fq).
#
# Output: <sample>.HS1_chrR.strReadsPerGene.out.tab
#         -> results/02_rRNA_quantification/SMARTseq_nbGFP_STAR/

GENOME=/project/hnisz_seq_data/FANCI_seq/RNA_seq/NUCL_RNA_20250721/rRNA_workflow_bowtie/genome/STAR

for sample in FANCI_nbGFP_2xKS_Dox_HCT116_Pb FANCI_nbGFP_2xKS_noDox_HCT116_Pb \
              Input_nbGFP_2xKS_Dox_HCT116_Pb Input_nbGFP_2xKS_noDox_HCT116_Pb; do
  for i in 1 2 3; do

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
