
library(CellChat)
library(Seurat)

seu=readRDS("DC_ent_singlet.rds")
#seu=readRDS(path)

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
#write.table(df.net,"cellchat.result.txt",sep="\t",quote=F)



saveRDS(cellchat,"dc_epi.cellchat.2026.rds")



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
# [1] "Tissue-resident inflammatory cDC2" "epi_1"                            
# [3] "Migratory LAMP3hi cDC2"            "epi_2"                            
# [5] "Migratory inflammatory cDC2"    

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

