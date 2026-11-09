#!/usr/bin/env bash
# Step 02c: remove domains shorter than the N95 domain size (from step 02b), then sort.
# Split from the original run.sh.

## HCT116 nucleoli (consensus)
awk '{if($3-$2 >= 61000) print}' NUCL_HCT116_GFPNPM1_GFP_sort_merged.bed > NUCL_HCT116_GFPNPM1_GFP_sort_merged.flt.bed
awk '{if($3-$2 >= 66000) print}' NUC_HCT116_NPM1_Ab_sort_merged.bed      > NUC_HCT116_NPM1_Ab_sort_merged.flt.bed
awk '{if($3-$2 >= 67000) print}' NUC_HCT116_GFPNPM1_Ab_sort_merged.bed   > NUC_HCT116_GFPNPM1_Ab_sort_merged.flt.bed

## HCT116 speckles (consensus)
awk '{if($3-$2 >= 24000) print}' NPKL_DNA_FANCI_merged.w50k_b1k.bed | sort -k1,1 -k2,2n > NPKL_DNA_FANCI_merged.w50k_b1k.flt.sort.bed

## Remove any final domain overlapping a CHM13 centromere (same intervals as FPCaller).
## FPCaller already filters per replicate; this re-applies the filter to the final consensus sets.
CEN=$(dirname "$0")/../../references/chm13_centromeres.bed
for f in NUCL_HCT116_GFPNPM1_GFP_sort_merged.flt.bed NUC_HCT116_NPM1_Ab_sort_merged.flt.bed \
         NUC_HCT116_GFPNPM1_Ab_sort_merged.flt.bed NPKL_DNA_FANCI_merged.w50k_b1k.flt.sort.bed; do
  bedtools intersect -v -a "$f" -b "$CEN" > "${f%.bed}.tmp" && mv "${f%.bed}.tmp" "$f"
done
