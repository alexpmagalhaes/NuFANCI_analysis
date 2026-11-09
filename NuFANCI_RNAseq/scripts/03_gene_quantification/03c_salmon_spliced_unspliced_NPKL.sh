#!/usr/bin/env bash
# Step 03c: transcript quantification of NPKL RNA-seq (replicates 2, 4, 5) with Salmon
# against a CHM13 index that contains both spliced and unspliced (pre-mRNA) transcripts,
# to separate mature from nascent RNA.
#
# Output: salmon_counts/<sample>_quant/quant.sf (on GEO); used in NuFANCI_DNAseq step 07.


cd /project/hnisz_seq_data/FANCI_seq/RNA_seq/NPKL_RNA_20250721

mkdir -p salmon_counts

salmon quant -i /project/hnisz_bioinf/Genomes/CHM13/Genes/Homo_sapiens_T2T_splice_unspliced -l A \
-1 ./rnaseq_NPKL_RNA_FANCI_rep2_chm13/NPKL_RNA_FANCI_rep2_R1.trm.fq.gz \
-2 ./rnaseq_NPKL_RNA_FANCI_rep2_chm13/NPKL_RNA_FANCI_rep2_R2.trm.fq.gz \
--validateMappings -o ./salmon_counts/NPKL_RNA_FANCI_rep2_quant

salmon quant -i /project/hnisz_bioinf/Genomes/CHM13/Genes/Homo_sapiens_T2T_splice_unspliced -l A \
-1 ./rnaseq_NPKL_RNA_FANCI_rep4_chm13/NPKL_RNA_FANCI_rep4_R1.trm.fq.gz \
-2 ./rnaseq_NPKL_RNA_FANCI_rep4_chm13/NPKL_RNA_FANCI_rep4_R2.trm.fq.gz \
--validateMappings -o ./salmon_counts/NPKL_RNA_FANCI_rep4_quant

salmon quant -i /project/hnisz_bioinf/Genomes/CHM13/Genes/Homo_sapiens_T2T_splice_unspliced -l A \
-1 ./rnaseq_NPKL_RNA_FANCI_rep5_chm13/NPKL_RNA_FANCI_rep5_R1.trm.fq.gz \
-2 ./rnaseq_NPKL_RNA_FANCI_rep5_chm13/NPKL_RNA_FANCI_rep5_R2.trm.fq.gz \
--validateMappings -o ./salmon_counts/NPKL_RNA_FANCI_rep5_quant

salmon quant -i /project/hnisz_bioinf/Genomes/CHM13/Genes/Homo_sapiens_T2T_splice_unspliced -l A \
-1 ./rnaseq_NPKL_RNA_Input_rep2_chm13/NPKL_RNA_Input_rep2_R1.trm.fq.gz \
-2 ./rnaseq_NPKL_RNA_Input_rep2_chm13/NPKL_RNA_Input_rep2_R2.trm.fq.gz \
--validateMappings -o ./salmon_counts/NPKL_RNA_Input_rep2_quant

salmon quant -i /project/hnisz_bioinf/Genomes/CHM13/Genes/Homo_sapiens_T2T_splice_unspliced -l A \
-1 ./rnaseq_NPKL_RNA_Input_rep4_chm13/NPKL_RNA_Input_rep4_R1.trm.fq.gz \
-2 ./rnaseq_NPKL_RNA_Input_rep4_chm13/NPKL_RNA_Input_rep4_R2.trm.fq.gz \
--validateMappings -o ./salmon_counts/NPKL_RNA_Input_rep4_quant

salmon quant -i /project/hnisz_bioinf/Genomes/CHM13/Genes/Homo_sapiens_T2T_splice_unspliced -l A \
-1 ./rnaseq_NPKL_RNA_Input_rep5_chm13/NPKL_RNA_Input_rep5_R1.trm.fq.gz \
-2 ./rnaseq_NPKL_RNA_Input_rep5_chm13/NPKL_RNA_Input_rep5_R2.trm.fq.gz \
--validateMappings -o ./salmon_counts/NPKL_RNA_Input_rep5_quant
