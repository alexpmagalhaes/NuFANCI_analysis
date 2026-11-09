#!/usr/bin/env bash

# bedtools intersect -nonamecheck -wb -a /Users/magalhae/Desktop/FANCY_DNAcaller/FANCI_DNAseq_T2T_CHM3/BED/HCT116_NPKL/NPKL_DNA_FANCI_merged.w50k_b1k.flt.sort.bed -b ./chm13.draft_v2.0.gene_annotation.intronicAdded.gff > NPKL_DNA_FANCI_SADs_genes.gff
# bedtools complement -i /Users/magalhae/Desktop/FANCY_DNAcaller/FANCI_DNAseq_T2T_CHM3/BED/HCT116_NPKL/NPKL_DNA_FANCI_merged.w50k_b1k.flt.sort.bed -g /Users/magalhae/Desktop/FANCY_DNAcaller/FANCI_DNAseq_T2T_CHM3/chm13_chrm.sizes > NPKL_DNA_FANCI_merged.w50k_b1k.complement.bed
#
# bedtools intersect -nonamecheck -wb -a NPKL_DNA_FANCI_merged.w50k_b1k.complement.bed -b ./chm13.draft_v2.0.gene_annotation.intronicAdded.gff > NPKL_DNA_FANCI_All_otherGenes.gff
# gcut --complement -f1-3 ./NPKL_DNA_FANCI_All_otherGenes.gff > NPKL_DNA_FANCI_All_otherGenes.ftl.gff
#
# gcut --complement -f1-3 ./NPKL_DNA_FANCI_SADs_genes.gff > NPKL_DNA_FANCI_SADs_genes.ftl.gff

GENES=/Users/magalhae/Desktop/FANCY_DNAcaller/FANCI_DNAseq_T2T_CHM3/NPKL_SAD_genes/References/chm13.draft_v2.0.genes.bed
SADS=/Users/magalhae/Desktop/FANCY_DNAcaller/FANCI_DNAseq_T2T_CHM3/BED/HCT116_NPKL/NPKL_DNA_FANCI_merged.w50k_b1k.flt.sort.bed
NADS=/Users/magalhae/Desktop/FANCY_DNAcaller/FANCI_DNAseq_T2T_CHM3/BED/HCT116_NUCL/NUCL_HCT116_GFPNPM1_GFP_sort_merged.flt.bed

bedtools intersect -nonamecheck -wa -a $GENES -b $SADS > NPKL_DNA_FANCI_SADs_genes.bed
bedtools intersect -nonamecheck -wa -a $GENES -b $NADS > NUCL_DNA_FANCI_NADS_genes.bed

cat NPKL_DNA_FANCI_SADs_genes.bed NUCL_DNA_FANCI_NADS_genes.bed > NAD_SAD_genes.bed
sort -k1,1 -k2,2n NAD_SAD_genes.bed > NAD_SAD_genes.srt.bed

bedtools intersect -nonamecheck -v -a $GENES -b NAD_SAD_genes.srt.bed > Other_genes.bed




# awk -F "\t" '$3 ~ /exon/' NPKL_DNA_FANCI_SADs_genes.ftl.gff |                            # extract all lines of type 'exon' (3rd column)
# grep gene= |                                                                # delete lines without "gene=" (<50 lines that allocate to mitochondrial DNA, removal allows uniform format to create SAF)
# awk '{FS="\t"; OFS="\t"; sub(/.*;gene=/,"",$9); sub(/;.*/,"",$9)}1' |       # reduce column 9 to only include the value of "gene="
# awk 'BEGIN {FS="\t"; OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $9,$1,$4,$5,$7}' > NPKL_DNA_FANCI_SADs_genes.ftl.exon.saf
#
#
# awk -F "\t" '$3 ~ /intron/' NPKL_DNA_FANCI_SADs_genes.ftl.gff |                            # extract all lines of type 'exon' (3rd column)
# grep gene= |                                                                # delete lines without "gene=" (<50 lines that allocate to mitochondrial DNA, removal allows uniform format to create SAF)
# awk '{FS="\t"; OFS="\t"; sub(/.*;gene=/,"",$9); sub(/;.*/,"",$9)}1' |       # reduce column 9 to only include the value of "gene="
# awk 'BEGIN {FS="\t"; OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $9,$1,$4,$5,$7}' > NPKL_DNA_FANCI_SADs_genes.ftl.intron.saf
