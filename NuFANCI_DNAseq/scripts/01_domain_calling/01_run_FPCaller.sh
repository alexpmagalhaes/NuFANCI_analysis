#!/usr/bin/env bash
# Step 01: call NuFANCI-enriched domains (NADs for nucleolar sorts, SADs for speckle/NPKL sorts)
# from 50 kb binned bedGraphs of NuFANCI vs. input (CHM13 / T2T reference).
#
# Domain caller: FPCaller v1.0.0, https://github.com/alexpmagalhaes/FPCaller
# Install:  pip install "git+https://github.com/alexpmagalhaes/FPCaller.git@v1.0.0"
#
# The CHM13 centromere filter is on by default.
#
# Per sample, fpcaller writes:
#   <sample>.bed                     domains (positive baseline-corrected z-score)
#   <sample>.stats.tsv               n domains, total bp, mean size, N10/N50/N90, mean gap
#   <sample>.chr14_qc.{svg,png,pdf}  chr14 QC plot
# (The original script named these <sample>.stats.txt and <sample>.seaborn_plot.{svg,png};
#  results/01_domain_calling/ keeps those original files.)
#
# Input bedGraphs (*.w50k_b1k.bedgraph) are available on GEO.
# Split from the original run.sh.

BG=${BG:-./bedgraph}

run () { fpcaller --genome chm13 --fancy "$1" --input "$2" --sample "$3"; }

## HCT116 nucleoli (NUCL)
for s in NUC_GFPNPM1_Ab_sort NUC_GFPNPM1_GFP_sort NUC_NPM1_Ab_sort; do
  for r in 1 2; do
    run $BG/HCT116_NUCL/reps/${s}_rep${r}.w50k_b1k.bedgraph \
        $BG/HCT116_NUCL/reps/${s}_rep${r}_input.w50k_b1k.bedgraph \
        ${s}_rep${r}
  done
done

## HCT116 nuclear speckles (NPKL). rep3 excluded (low quality).
for r in 1 2 4; do
  run $BG/HCT116_NPKL/reps/NPKL_DNA_FANCI_rep${r}.w50k_b1k.bedgraph \
      $BG/HCT116_NPKL/reps/NPKL_DNA_FANCI_rep${r}_input.w50k_b1k.bedgraph \
      NPKL_DNA_FANCI_rep${r}.w50k_b1k
done

## HepG2 nucleoli (NPM1 Ab sort)
for r in 1 2 3 4; do
  run $BG/HepG2_NUCL/NUC_NPM1_Ab_sort_HepG2_R${r}.CHM13.w50k_b1k.bedgraph \
      $BG/HepG2_NUCL/NUC_NPM1_Ab_sort_HepG2_R${r}_input.CHM13.w50k_b1k.bedgraph \
      NUC_NPM1_Ab_sort_HepG2_R${r}
done
