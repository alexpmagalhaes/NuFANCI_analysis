#!/usr/bin/env bash
# Step 02b: N-statistics (N5..N95) of consensus domain sizes.
# The N95 value sets the minimum domain size used in step 02c.
# Split from the original run.sh.

NSTATS=$(dirname "$0")/n_stats.py

python $NSTATS NUC_HCT116_NPM1_Ab_sort_merged.bed       > NUC_HCT116_NPM1_Ab_sort_merged.n_stats.txt
python $NSTATS NUC_HCT116_GFPNPM1_Ab_sort_merged.bed    > NUC_HCT116_GFPNPM1_Ab_sort_merged.n_stats.txt
python $NSTATS NUCL_HCT116_GFPNPM1_GFP_sort_merged.bed  > NUCL_HCT116_GFPNPM1_GFP_sort_merged.n_stats.txt
python $NSTATS NPKL_DNA_FANCI_merged.w50k_b1k.bed       > NPKL_DNA_FANCI_merged.w50k_b1k.stats.tsv
python $NSTATS NUC_NPM1_Ab_sort_HepG2_merged.bed        > NUC_NPM1_Ab_sort_HepG2_merged.n_stats.txt
