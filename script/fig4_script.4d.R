

obj=readRDS("rename_epi.dc.330.2026(input_4d).rds")
obj$celltype.final=gsub("epi_LYZ_Low","Epithelial 2",obj$celltype.final)

obj$celltype.final=gsub("epi_FABP_LYZ_High","Epithelial 1",obj$celltype.final)


pdf("fig.s4.marker.pdf",15,10)
gene <- c(
  "AGR2", "FABP1", "FABP2", "EPCAM", "MUC2",
  "CD14", "LYZ", "SIRPA", "LAMP3", "CD1A",
  "CD1C", "HLA-DRA", "CLEC10A", "ITGAX"
)

FeaturePlot_scCustom(integrated,features=gene)
#,cols = viridis::viridis(100))
dev.off()


#saveRDS(integrated,"input_fig4d.rds")







