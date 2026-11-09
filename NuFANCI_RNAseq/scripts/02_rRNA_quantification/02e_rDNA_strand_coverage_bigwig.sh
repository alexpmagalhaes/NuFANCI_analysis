#!/usr/bin/env bash
# Step 02e: strand-specific rDNA coverage tracks (bigWig) from the STAR rDNA alignments
# of step 02b, normalised as BPM and CPM, with deepTools bamCoverage.
#
# Output: <bam>.cov.{BPM,CPM}.{for,rev}.bw  (on GEO)
# effectiveGenomeSize = length of the hs1-rDNA reference (chrR).

for bam_file in NUC_RNA_*_rep*.HS1_chrR.srt.bam; do
  samtools index ${bam_file}
  for n in BPM CPM; do
    for f in forward reverse; do
      bamCoverage \
        --bam ${bam_file} \
        --normalizeUsing ${n} \
        --filterRNAstrand ${f} \
        --outFileName ${bam_file}.cov.${n}.${f:0:3}.bw \
        --numberOfProcessors 30 \
        --effectiveGenomeSize 44838
    done
  done
done
