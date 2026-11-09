# =============================================================================
# NuFANCI ONT: CpG methylation at rDNA genes, HCT116 NPM1 sort replicate 2
#
# Author: Henri Niskanen
# email:  niskanen@molgen.mpg.de
#
# Plots CpG methylation of nucleolar NuFANCI-sorted DNA vs input (ONT) across the
# 45S rDNA genes.
#
# Input:
#   - CpG methylation bigWigs (from modkit bedMethyl; CG, + strand; % methylation in "score")
#   - catLiftOffGenesV1.bed: T2T-CHM13 v2.0 gene annotation (CAT/Liftoff), used to
#     select the 45S rDNA gene copies RNA45S1-1 ... RNA45S5-1
# Output:
#   - CG_meth_rDNA_genes_heatmap*.pdf: heatmaps of CpG methylation (0-100 %) around
#     the rDNA genes: (1) gene body +/- 15 kb; (2) 30 kb upstream (IGS/promoter) + gene body
#
# Requires: R packages EnrichedHeatmap, circlize, rtracklayer
# Note: human liver samples were removed from the published version of this script.
# =============================================================================

library(EnrichedHeatmap)
library(circlize)
library(rtracklayer)


FANCI_HCT_R2_bw <- import("HCT116_NPM1_sort_FANCI_rep2_CG_positive.sorted.bw")
input_HCT_R2_bw <- import("HCT116_NPM1_sort_Input_rep2_CG_positive.sorted.bw")

# Plot CpG methylation across individual rDNA gene loci

genes <- import("catLiftOffGenesV1.bed")

rRNA_genes <- genes[genes$name %in% c("RNA45S1-1","RNA45S2-1","RNA45S3-1","RNA45S4-1","RNA45S5-1")]


FANCI_HCT_R2_bw <- import("HCT116_NPM1_sort_FANCI_rep2_CG_positive.sorted.bw")
input_HCT_R2_bw <- import("HCT116_NPM1_sort_Input_rep2_CG_positive.sorted.bw")


mat_HCT_FANCI_rRNA_genes <- normalizeToMatrix(FANCI_HCT_R2_bw, rRNA_genes, value_column = "score", extend = c(15000,15000), mean_mode = "absolute", w = 500, empty_value = NA, target_ratio = 0.5)
mat_HCT_input_rRNA_genes <- normalizeToMatrix(input_HCT_R2_bw, rRNA_genes, value_column = "score", extend = c(15000,15000), mean_mode = "absolute", w = 500, empty_value = NA, target_ratio = 0.5)


col_fun <- colorRamp2(c(0, 50, 100), c("blue", "white", "red"))


pdf("CG_meth_rDNA_genes_heatmap_HCT116_NPM1sort_rep2.pdf", width=20, height=5)

  EnrichedHeatmap(mat_HCT_FANCI_rRNA_genes, col = col_fun, axis_name_rot= 90, column_title = "HCT_R2_FANCI", name="HCT_R2_FANCI", show_heatmap_legend=T,
                  top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0,100))), width=unit(3.5,"cm")) +
  EnrichedHeatmap(mat_HCT_input_rRNA_genes, col = col_fun, axis_name_rot= 90, column_title = "HCT_R2_input", name="HCT_R2_input", show_heatmap_legend=T,
                  top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0,100))), width=unit(3.5,"cm")) 

dev.off()


mat2_HCT_FANCI_rRNA_genes <- normalizeToMatrix(FANCI_HCT_R2_bw, rRNA_genes, value_column = "score", extend = c(30000,0), mean_mode = "absolute", w = 200, empty_value = NA, target_ratio = 0.3)
mat2_HCT_input_rRNA_genes <- normalizeToMatrix(input_HCT_R2_bw, rRNA_genes, value_column = "score", extend = c(30000,0), mean_mode = "absolute", w = 200, empty_value = NA, target_ratio = 0.3)

col_fun <- colorRamp2(c(0, 50, 100), c("blue", "white", "red"))


pdf("CG_meth_rDNA_genes_heatmap2_HCT116_NPM1sort_rep2.pdf", width=20, height=5)

  EnrichedHeatmap(mat2_HCT_FANCI_rRNA_genes, col = col_fun, axis_name_rot= 90, column_title = "HCT_R2_FANCI", name="HCT_R2_FANCI", show_heatmap_legend=T,
                  top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0,100))), width=unit(3.5,"cm")) +
  EnrichedHeatmap(mat2_HCT_input_rRNA_genes, col = col_fun, axis_name_rot= 90, column_title = "HCT_R2_input", name="HCT_R2_input", show_heatmap_legend=T,
                  top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0,100))), width=unit(3.5,"cm")) 

dev.off()
