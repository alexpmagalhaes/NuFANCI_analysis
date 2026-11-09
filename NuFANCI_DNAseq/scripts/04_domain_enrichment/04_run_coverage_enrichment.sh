
cd /project/hnisz_seq_data/FANCI_seq/DNA_seq/HCT116_NUCL_DNA_20250831/bed_coverage
# megadepth=/project/hnisz_bioinf/bin/megadepth

# BED_dir=/project/hnisz_seq_data/FANCI_seq/DNA_seq/HCT116_NUCL_DNA_20250831/BED/HCT116_NUCL
# BAM_dir=/project/hnisz_seq_data/FANCI_seq/DNA_seq/HCT116_NUCL_DNA_20250831/bam
#
# # bedtools coverage -a ${BED_dir}/NUC_HCT116_NPM1_Ab_sort_merged.flt.bed -b file.bam > file.cov.txt
# # bedtools coverage -a ${BED_dir}/NUC_HCT116_GFPNPM1_Ab_sort_merged.flt.bed -b file.bam > file.cov.txt
# # bedtools coverage -a ${BED_dir}/NUCL_HCT116_GFPNPM1_GFP_sort_merged.flt.bed -b file.bam > file.cov.txt
#
#
# multiBamSummary BED-file \
#   --BED ${BED_dir}/NUC_HCT116_NPM1_Ab_sort_merged.flt.bed \
#   --outRawCounts NUC_HCT116_NPM1_Ab_sort_merged.coverage.tsv \
#   --samFlagExclude 256 \
#   --smartLabels \
#   --numberOfProcessors 30 \
#   -o NUC_HCT116_NPM1_Ab_sort_merged.coverage.npz \
#   --bamfiles ${BAM_dir}/NUC_NPM1_Ab_sort_rep1_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_NPM1_Ab_sort_rep2_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_NPM1_Ab_sort_rep1_input_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_NPM1_Ab_sort_rep2_input_CHM13.bwa.srt.md.bam
#
# multiBamSummary BED-file \
#   --BED ${BED_dir}/NUC_HCT116_GFPNPM1_Ab_sort_merged.flt.bed \
#   --outRawCounts NUC_HCT116_GFPNPM1_Ab_sort_merged.coverage.tsv \
#   --samFlagExclude 256 \
#   --smartLabels \
#   --numberOfProcessors 12 \
#   -o NUC_HCT116_GFPNPM1_Ab_sort_merged.coverage.npz \
#   --bamfiles ${BAM_dir}/NUC_GFPNPM1_GFP_sort_rep1_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_GFPNPM1_GFP_sort_rep2_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_GFPNPM1_GFP_sort_rep1_input_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_GFPNPM1_GFP_sort_rep2_input_CHM13.bwa.srt.md.bam
#
# multiBamSummary BED-file \
#   --BED ${BED_dir}/NUCL_HCT116_GFPNPM1_GFP_sort_merged.flt.bed \
#   --outRawCounts NUC_HCT116_GFPNPM1_GFP_sort_merged.coverage.tsv\
#   --samFlagExclude 256 \
#   --smartLabels \
#   --numberOfProcessors 12 \
#   -o NUC_HCT116_GFPNPM1_GFP_sort_merged.coverage.npz \
#   --bamfiles ${BAM_dir}/NUC_GFPNPM1_Ab_sort_rep1_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_GFPNPM1_Ab_sort_rep2_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_GFPNPM1_Ab_sort_rep1_input_CHM13.bwa.srt.md.bam ${BAM_dir}/NUC_GFPNPM1_Ab_sort_rep2_input_CHM13.bwa.srt.md.bam
#


python ./run_enrichment_NADs.py -i NUC_HCT116_GFPNPM1_GFP_sort_merged.coverage.tsv \
-o NUC_HCT116_GFPNPM1_GFP_sort_merged.coverage.pval.csv --coords chr start end \
--inputs NUC_GFPNPM1_Ab_sort_rep1_input NUC_GFPNPM1_Ab_sort_rep2_input \
--ips NUC_GFPNPM1_Ab_sort_rep1 NUC_GFPNPM1_Ab_sort_rep2


python ./run_enrichment_NADs.py -i NUC_HCT116_GFPNPM1_Ab_sort_merged.coverage.tsv \
-o NUC_HCT116_GFPNPM1_Ab_sort_merged.coverage.pval.csv --coords chr start end \
--inputs NUC_GFPNPM1_GFP_sort_rep1_input NUC_GFPNPM1_GFP_sort_rep2_input \
--ips NUC_GFPNPM1_GFP_sort_rep1 NUC_GFPNPM1_GFP_sort_rep2


python ./run_enrichment_NADs.py -i NUC_HCT116_NPM1_Ab_sort_merged.coverage.tsv \
-o NUC_HCT116_NPM1_Ab_sort_merged.coverage.pval.csv --coords chr start end \
--inputs NUC_NPM1_Ab_sort_rep1_input NUC_NPM1_Ab_sort_rep2_input \
--ips NUC_NPM1_Ab_sort_rep1 NUC_NPM1_Ab_sort_rep2
