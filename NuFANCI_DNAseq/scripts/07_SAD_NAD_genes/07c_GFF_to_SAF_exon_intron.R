


library(Rsubread)


setwd("/Users/magalhae/Desktop/FANCY_DNAcaller/FANCI_DNAseq_T2T_CHM3/NPKL_SAD_genes/")
SAF_exon <- flattenGTF("NPKL_DNA_FANCI_SADs_genes.ftl.gff", GTF.featureType = "exon", GTF.attrType = "Name")
SAF_intron <- flattenGTF("NPKL_DNA_FANCI_SADs_genes.ftl.gff", GTF.featureType = "intron", GTF.attrType = "Name")

write.table(SAF_exon, file="NPKL_DNA_FANCI_SADs_genes.ftl.exon.saf", quote=FALSE, sep='\t', col.names = F, row.names = F)
write.table(SAF_intron, file="NPKL_DNA_FANCI_SADs_genes.ftl.intron.saf", quote=FALSE, sep='\t', col.names = F, row.names = F)


SAF_exon_others <- flattenGTF("NPKL_DNA_FANCI_All_otherGenes.ftl.gff", GTF.featureType = "exon", GTF.attrType = "Name")
SAF_intron_others <- flattenGTF("NPKL_DNA_FANCI_All_otherGenes.ftl.gff", GTF.featureType = "intron", GTF.attrType = "Name")

write.table(SAF_exon_others, file="NPKL_DNA_FANCI_All_otherGenes.ftl.exon.saf", quote=FALSE, sep='\t', col.names = F, row.names = F)
write.table(SAF_intron_others, file="NPKL_DNA_FANCI_All_otherGenes.ftl.intron.saf", quote=FALSE, sep='\t', col.names = F, row.names = F)
