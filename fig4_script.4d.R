dc=readRDS("/local/projects-t3/BTRAN/analysis/bma_project/snrna_bing_1106/integrate_two_DC_anchor.base/result_813/result_825/result_825/result_2026/integrated.DC.only.proportion.rds")
DefaultAssay(dc)="RNA"
dc1 <- DietSeurat(dc, assays = "RNA")
#combined <- merge(
# x = epi,
#y = dc1
# add.cell.ids = c("s1", "s2", "s3", "s4"))


dc1 <- SCTransform(dc1, vst.flavor = "v2", verbose = FALSE)
dc2 <- SCTransform(epi, vst.flavor = "v2", verbose = FALSE)

data1=RunPCA(dc1,npcs = 30, verbose = FALSE)
data2=RunPCA(dc2,npcs = 30, verbose = FALSE)


data.list <- list(data1, data2)
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)

data.list[[1]]$sample_id <- "dc"
data.list[[2]]$sample_id <- "epi"

data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features, verbose = TRUE)

# Step 4: Find integration anchors
anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", 
                                  anchor.features = features, verbose = TRUE)

# Step 5: Integrate the datasets
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = TRUE)

# Step 6: Proceed with downstream analysis on integrated data
# Set the default assay to the integrated data
DefaultAssay(integrated) <- "integrated"

# Run PCA
integrated <- RunPCA(integrated, verbose = TRUE)

# Run UMAP for visualization
integrated <- FindNeighbors(integrated, dims = 1:30)

integrated <- RunUMAP(integrated, dims = 1:30, verbose = TRUE)

integrated <- FindClusters(integrated, resolution = 1.0, graph.name = "integrated_snn")
# Find neighbors and clusters

saveRDS(integrated,"integrated.epi.dc.2026.rds")

integrated=readRDS("/local/projects-t3/BTRAN/analysis/bma_project/snrna_bing_1106/integrate_two_DC_anchor.base/result_813/result_825/result_825/result_2026/result_329/rename_epi.dc.330.2026.rds")
integrated$celltype.final=gsub("epi_LYZ_Low","Epithelial 2",integrated$celltype.final)

integrated$celltype.final=gsub("epi_FABP_LYZ_High","Epithelial 1",integrated$celltype.final)


pdf("fig.s4.marker.pdf",15,10)
gene <- c(
  "AGR2", "FABP1", "FABP2", "EPCAM", "MUC2",
  "CD14", "LYZ", "SIRPA", "LAMP3", "CD1A",
  "CD1C", "HLA-DRA", "CLEC10A", "ITGAX"
)

FeaturePlot_scCustom(integrated,features=gene)
#,cols = viridis::viridis(100))
dev.off()


saveRDS(integrated,"input_fig4d.rds")



seu=integrated

DefaultAssay(seu) <- "RNA"

data.input <- GetAssayData(seu, assay = "RNA", layer = "data")
meta <- seu@meta.data[, "celltype.final", drop = FALSE]
colnames(meta) <- "celltype.final"

cc <- createCellChat(object = data.input, meta = meta, group.by = "celltype.final")
cc@DB <- CellChatDB.human

cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
cc <- computeCommunProb(cc)
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)

cellchat=cc


df.net <- subsetCommunication(cellchat)
write.table(df.net,"cellchat.result.txt",sep="\t",quote=F)


cellchat@net$count <- rename_dimnames(cellchat@net$count)
cellchat@net$weight <- rename_dimnames(cellchat@net$weight)



pdf("cellchat.fig1.pdf",5,5)
groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,vertex.label.cex = 0.5,
  title.name = "Interaction weights/strength"
)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,vertex.label.cex = 0.5,
  title.name = "Interaction weights/strength")

dev.off()

pdf("pathway.pdf")

mat <- cellchat@net$weight
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize,vertex.label.cex = 0.5,
                   weight.scale = T, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

dev.off()

saveRDS(cc,"dc_epi.input.cellchat.2026.rds")


pathways.show <- cellchat@netP$pathways
pathways.show
length(pathways.show)


pdf("all_pathways_circle.pdf", width = 8, height = 8)

for (pw in pathways.show) {
  netVisual_aggregate(cellchat, signaling = pw, layout = "circle")
  title(main = pw)
}


dev.off()


pdf("all_pathways_chord.pdf", width = 8, height = 8)

for (pw in pathways.show) {
  circlize::circos.clear()
  netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
  #title(main = pw)
}

dev.off()


pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
#vertex.receiver = seq(1,4)

pdf("path.hirachy.829.pdf",15,5)

path=cellchat@netP$pathways
plot_list = list()
for (i in 1:length(path)) {
  pathways.show=path[i]
  vertex.receiver <- c(1,4,5) # define a numeric vector giving the index of the cell type as targets
  #par(mar=c(5.1,4.1,4.1,2.1))
  p=netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy",title.space = 1,show.legend = TRUE,small.gap =10)
  plot_list[[i]] = p}
dev.off()



pdf("all_pathway_chord_gene_epi_to_dc.pdf", width = 8, height = 8)

for (pw in path) {
  circlize::circos.clear()
  
  tryCatch({
    netVisual_chord_gene(
      cellchat,
      signaling = pw,
      sources.use = c(2,3),
      targets.use = c(1,4,5),
      lab.cex = 0.5,
      legend.pos.x = 1,legend.pos.y = 5,
      small.gap = 10,
      title.name = paste0(pw, " signaling")
    )
  }, error = function(e) {
    message("Skipping ", pw, ": ", e$message)
  })
}

dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf("incom_outgoing.pdf",10,6)

ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",
  font.size = 4,
  font.size.title = 10
)

ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming",
  font.size = 4,
  font.size.title = 10
)
library(ComplexHeatmap)
draw(ht1+ht2, padding = unit(c(2, 2, 2, 2), "mm"))
#draw(ht2, padding = unit(c(8, 8, 8, 8), "mm"))


dev.off()


pdf("bubble.pdf",10,15)
netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(5:10), remove.isolate = FALSE)
#> Comparing communications on a single object

dev.off()


pdf("chord.lr.pdf",10,10)

netVisual_chord_gene(cellchat, sources.use = c(2,3), targets.use = c(1,4,5), 
                     lab.cex = 0.5, legend.pos.x = 1,
                     legend.pos.y = 5)

dev.off()

pdf("all_pathway_chord.pdf", width = 10, height = 10)

for (pw in path) {
  circlize::circos.clear()
  
  tryCatch({
    netVisual_chord_gene(
      cellchat,
      signaling = pw,
      sources.use = c(2,3),
      targets.use = c(1,4,5),
      lab.cex = 0.8,
      small.gap = 6,
      title.name = pw
    )
  }, error = function(e) {
    message("Skipping ", pw, ": ", e$message)
  })
}

dev.off()


saveRDS(cellchat,"dc_epi.cellchat.2026.rds")




df.all <- subsetCommunication(cellchat)




for (pw in unique(df.all$pathway_name)) {
  df.pw <- subset(df.all, pathway_name == pw)
  
  if (nrow(df.pw) == 0) next
  
  pdf(paste0( pw, "_custom_chord.pdf"), width = 10, height = 10)
  circlize::circos.clear()
  
  chordDiagram(df.pw[, c("ligand", "receptor", "prob")])
  
  title(main = pw)
  dev.off()
}

