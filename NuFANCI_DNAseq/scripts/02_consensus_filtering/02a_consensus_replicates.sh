#!/usr/bin/env bash
# Step 02a: consensus domains = intersection of replicate domain BEDs (from step 01).
# Split from the original run.sh.

bedtools intersect -a NUC_GFPNPM1_Ab_sort_rep1.bed  -b NUC_GFPNPM1_Ab_sort_rep2.bed  > NUC_HCT116_GFPNPM1_Ab_sort_merged.bed
bedtools intersect -a NUC_GFPNPM1_GFP_sort_rep1.bed -b NUC_GFPNPM1_GFP_sort_rep2.bed > NUCL_HCT116_GFPNPM1_GFP_sort_merged.bed
bedtools intersect -a NUC_NPM1_Ab_sort_rep1.bed     -b NUC_NPM1_Ab_sort_rep2.bed     > NUC_HCT116_NPM1_Ab_sort_merged.bed

bedtools intersect -a NPKL_DNA_FANCI_rep1.w50k_b1k.bed \
  -b NPKL_DNA_FANCI_rep2.w50k_b1k.bed NPKL_DNA_FANCI_rep4.w50k_b1k.bed > NPKL_DNA_FANCI_merged.w50k_b1k.bed

bedtools intersect -a NUC_NPM1_Ab_sort_HepG2_R1.bed \
  -b NUC_NPM1_Ab_sort_HepG2_R2.bed NUC_NPM1_Ab_sort_HepG2_R3.bed NUC_NPM1_Ab_sort_HepG2_R4.bed > NUC_NPM1_Ab_sort_HepG2_merged.bed
