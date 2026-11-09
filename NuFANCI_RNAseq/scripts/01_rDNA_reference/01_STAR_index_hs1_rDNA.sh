#!/usr/bin/env bash
# Step 01: build a STAR index for the single-copy human rDNA reference (hs1-rDNA v1.0, contig chrR).
# Used by all STAR rRNA quantifications in step 02.
#
# Inputs (references/ and external):
#   HS1_chrR.fa          hs1-rDNA v1.0 rDNA unit (external, see README)
#   hs1-rDNA_v1.0.gtf    one gene "rRNA" spanning the 45S transcription unit (references/)
#
# Command recovered from genomeParameters.txt of the index (STAR 2.7.4a).
# --genomeSAindexNbases 6 because the reference is only ~45 kb.

cd genome/STAR

STAR --runMode genomeGenerate \
  --runThreadN 30 \
  --genomeDir ./ \
  --genomeFastaFiles ./HS1_chrR.fa \
  --genomeSAindexNbases 6 \
  --sjdbGTFfile ./hs1-rDNA_v1.0.gtf \
  --sjdbOverhang 149
