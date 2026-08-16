dc=readRDS("integrated.DC.only.proportion.rds")
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

integrated=readRDS("rename_epi.dc.330.2026.rds")
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



saveRDS(cc,"dc_epi.input.cellchat.2026.rds")




