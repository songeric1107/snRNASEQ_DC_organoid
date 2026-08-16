
library(CellChat)
library(Seurat)

seu=readRDS("DC_ent_singlet.rds")
#seu=readRDS(path)

DefaultAssay(seu) <- "RNA"

data.input <- GetAssayData(seu, assay = "RNA", layer = "data")
meta <- seu@meta.data[, "celltype.final", drop = FALSE]
colnames(meta) <- "celltype.final"

cc <- createCellChat(object = data.input, meta = meta, group.by = "cluster_annot1")
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
#write.table(df.net,"cellchat.result.txt",sep="\t",quote=F)


cellchat@net$count <- rename_dimnames(cellchat@net$count)
cellchat@net$weight <- rename_dimnames(cellchat@net$weight)


cellchat=readRDS("dc_epi.cellchat.2026.rds")
cellchat=readRDS("/local/projects-t3/BTRAN/analysis/bma_project/snrna_bing_1106/integrate_two_DC_anchor.base/result_813/result_825/result_825/result_2026/result_329/celltype.v3/dc_epi.cellchat.2026.rds")


pdf("bubble.pdf",10,10)
netVisual_bubble(cellchat, sources.use = c(2,4), targets.use = c(1,3,5), remove.isolate = FALSE)
#> Comparing communications on a single object

dev.off()


rename_map <- c(
  "epi_FABP_LYZ_High" = "Epithelial 1",
  "epi_FABP_LYZ_Low" = "Epithelial 2"
)

old_levels <- levels(cellchat@idents)

new_levels <- ifelse(
  old_levels %in% names(rename_map),
  rename_map[old_levels],
  old_levels
)


df_all <- subsetCommunication(
  cellchat,
  sources.use = c(2, 4),
  targets.use = c(1, 3, 5),
  thresh = 1
)

nrow(df_all)                  # total interactions before any filtering
head(df_all)                  # columns: source, target, ligand, receptor,
# interaction_name, pathway_name, prob, pval

summary(df_all$prob)          # distribution of communication strength
summary(df_all$pval)          # distribution of significance

hist(df_all$prob, breaks = 40, main = "Communication probability")

sapply(c(0, 0.01, 0.02, 0.05, 0.1), function(p) sum(df_all$prob >= p))
sapply(c(0.05, 0.01, 0.005, 0.002, 0.001), function(t) sum(df_all$pval <= t))

df_filtered <- subset(df_all, prob >= 0.02 & pval <= 0.05)   # tune 0.02 based on what you saw above

netVisual_chord_gene(
  cellchat,
  net = df_filtered,
  sources.use = c(2, 4),
  targets.use = c(1, 3, 5),
  lab.cex = 1,
  small.gap = 1,
  big.gap = 20,
  legend.pos.x = 1.4,
  legend.pos.y = 5
)                                                   

levels(cellchat@idents)
[1] "Tissue-resident inflammatory cDC2" "epi_1"                            
[3] "Migratory LAMP3hi cDC2"            "epi_2"                            
[5] "Migratory inflammatory cDC2"    

pdf("fig4e.v5.pdf", width = 10, height = 10)

cell.colors <- c(
  "Migratory inflammatory cDC2"      = "orange",  # purple
  "Tissue-resident inflammatory cDC2"  = "#D62F2F",  # red
  "Migratory LAMP3hi cDC2"      = "purple",  # orange
  "epi_1"                = "#45A049",  # green
  "epi_2"                = "#3B75A5"   # blue
)

# Ensure colors follow the identity order in the CellChat object
cell.colors <- cell.colors[levels(cellchat@idents)]

netVisual_chord_gene(
  cellchat,
  sources.use = c(2, 4),
  targets.use = c(1, 3, 5),
  lab.cex = 0.8,thresh=0.01,
  color.use   = cell.colors, legend.pos.x = 1.4,
legend.pos.y = 5)
dev.off()



df_all <- subsetCommunication(
  cellchat,
  sources.use = c(2, 4),
  targets.use = c(1, 3, 5),
  thresh = 1
)



df_filtered <- subset(df_all, prob >= 0.01 & pval <= 0.01)   # tune 0.02 based on what you saw above



pdf("fig4e.v6.pdf", width = 10, height = 10)

cell.colors <- c(
  "Migratory inflammatory cDC2"      = "orange",  # purple
  "Tissue-resident inflammatory cDC2"  = "#D62F2F",  # red
  "Migratory LAMP3hi cDC2"      = "purple",  # orange
  "epi_1"                = "#45A049",  # green
  "epi_2"                = "#3B75A5"   # blue
)

# Ensure colors follow the identity order in the CellChat object
cell.colors <- cell.colors[levels(cellchat@idents)]

netVisual_chord_gene(
  cellchat,net = df_filtered,
  sources.use = c(2, 4),
  targets.use = c(1, 3, 5),
  lab.cex = 1,thresh=0.01,
  color.use   = cell.colors, legend.pos.x = 1.4,
  legend.pos.y = 5)
dev.off()

pdf("fig4e.v7.pdf", width = 10, height = 10)

cell.colors <- c(
  "Migratory inflammatory cDC2"      = "orange",  # purple
  "Tissue-resident inflammatory cDC2"  = "#D62F2F",  # red
  "Migratory LAMP3hi cDC2"      = "purple",  # orange
  "epi_2"                = "#45A049",  # green
  "epi_1"                = "#3B75A5"   # blue
)

# Ensure colors follow the identity order in the CellChat object
cell.colors <- cell.colors[levels(cellchat@idents)]
df_filtered <- subset(df_all, prob >= 0.01 & pval <= 0.01)   # tune 0.02 based on what you saw above


netVisual_chord_gene(
  cellchat,net=df_filtered,
  sources.use = c(2, 4),
  targets.use = c(1, 3, 5),
  lab.cex = 1.5,thresh=0.01,
  color.use   = cell.colors, legend.pos.x = 1.4,
  legend.pos.y = 5)
dev.off()

# pdf("cellchat.fig1.pdf",5,5)
# groupSize <- as.numeric(table(cellchat@idents))
# #par(mfrow = c(1,2))
# netVisual_circle(
#   cellchat@net$weight,
#   vertex.weight = groupSize,
#   weight.scale = TRUE,
#   label.edge = FALSE,vertex.label.cex = 0.5,
#   title.name = "Interaction weights/strength"
# )
# 
# netVisual_circle(
#   cellchat@net$count,
#   vertex.weight = groupSize,
#   weight.scale = TRUE,
#   label.edge = FALSE,vertex.label.cex = 0.5,
#   title.name = "Interaction weights/strength")
# 
# dev.off()

# pdf("pathway.pdf")
# 
# mat <- cellchat@net$weight
# par(mfrow = c(2,2), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize,vertex.label.cex = 0.5,
#                    weight.scale = T, edge.weight.max = max(mat), 
#                    title.name = rownames(mat)[i])
# }
# 
# dev.off()

# saveRDS(cc,"dc_epi.input.cellchat.2026.rds")


pathways.show <- cellchat@netP$pathways
pathways.show
length(pathways.show)

# 
# pdf("all_pathways_circle.pdf", width = 8, height = 8)
# 
# for (pw in pathways.show) {
#   netVisual_aggregate(cellchat, signaling = pw, layout = "circle")
#   title(main = pw)
# }
# 
# 
# dev.off()
# 
# 
# pdf("all_pathways_chord.pdf", width = 8, height = 8)
# 
# for (pw in pathways.show) {
#   circlize::circos.clear()
#   netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
#   #title(main = pw)
# }
# 
# dev.off()


pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
#vertex.receiver = seq(1,4)

# pdf("path.hirachy.829.pdf",15,5)
# 
# path=cellchat@netP$pathways
# plot_list = list()
# for (i in 1:length(path)) {
#   pathways.show=path[i]
#   vertex.receiver <- c(1,4,5) # define a numeric vector giving the index of the cell type as targets
#   #par(mar=c(5.1,4.1,4.1,2.1))
#   p=netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy",title.space = 1,show.legend = TRUE,small.gap =10)
#   plot_list[[i]] = p}
# dev.off()

# 
# 
# pdf("all_pathway_chord_gene_epi_to_dc.pdf", width = 8, height = 8)
# 
# for (pw in path) {
#   circlize::circos.clear()
#   
#   tryCatch({
#     netVisual_chord_gene(
#       cellchat,
#       signaling = pw,
#       sources.use = c(2,3),
#       targets.use = c(1,4,5),
#       lab.cex = 0.5,
#       legend.pos.x = 1,legend.pos.y = 5,
#       small.gap = 10,
#       title.name = paste0(pw, " signaling")
#     )
#   }, error = function(e) {
#     message("Skipping ", pw, ": ", e$message)
#   })
# }
# 
# dev.off()
# 
# cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# 
# pdf("incom_outgoing.pdf",10,6)
# 
# ht1 <- netAnalysis_signalingRole_heatmap(
#   cellchat,
#   pattern = "outgoing",
#   font.size = 4,
#   font.size.title = 10
# )
# 
# ht2 <- netAnalysis_signalingRole_heatmap(
#   cellchat,
#   pattern = "incoming",
#   font.size = 4,
#   font.size.title = 10
# )
# library(ComplexHeatmap)
# draw(ht1+ht2, padding = unit(c(2, 2, 2, 2), "mm"))
# #draw(ht2, padding = unit(c(8, 8, 8, 8), "mm"))
# 
# 
# dev.off()
# 
# 




















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
