

library(Seurat)
#Seurat_4.0.2 
library(SeuratData)
library(cowplot)

#library(DoubletFinder)

Sys.setenv("R_C_STACK_LIMIT" = 500000)
suppressMessages(library(SeuratWrappers))
suppressMessages(library("dplyr"))
suppressMessages(library("Matrix"))
#suppressMessages(library("MAST"))
library(future)
options(future.globals.maxSize = 40 * 1024^3)
#options(Seurat.object.assay.version = 'v5')
# Load required libraries
library(Seurat)
library(sctransform)
library(dplyr)

# Set seed for reproducibility
set.seed(1234)



suppressMessages(library("dplyr"))
suppressMessages(library("Matrix"))

suppressMessages(library(SeuratWrappers))
library(future)
library(future.apply)

library(future)

plan(sequential)


options(future.globals.maxSize = 50 * 1024^3)


# load in raw matrix from DC-only and DC-enteroid -------------------------


countFolders <- c("/DC_only/filtered_feature_bc_matrix/","/DC_enteroid/filtered_feature_bc_matrix/")

samples <- c("DC-only","DC_enteroid"  )

Data <- list(NA)

for(i in 1:length(countFolders)){
  
  tmpfile <- countFolders[i]
  
  tmp <- Read10X(gsub("/matrix.mtx","",tmpfile))
  tmp <- as(tmp,Class="dgCMatrix")
  
  tmpcs <- CreateSeuratObject(counts = tmp, project = as.character(samples[i]), min.cells = 3, min.features = 200 )
  tmpcs@meta.data$type=tmpcs@meta.data$orig.ident
  
  Data[[i]] <- tmpcs
  filename <- paste(samples[i], ".rds", sep="")
  #saveRDS(tmpcs,filename)
  
  cat(paste("\n\n",as.character(samples[i])," done \n\n"))
}

##get mt percentage

mito.genes <- grep(pattern = "^MT-", x = rownames(x =GetAssayData(object = Data[[1]], slot = "counts")), value = TRUE)
length(mito.genes)
mito.genes1 <- grep(pattern = "^MT-", x = rownames(x =GetAssayData(object = Data[[2]], slot = "counts")), value = TRUE)
length(mito.genes1)

for(i in 1:length(Data)){
  mito.genes <- grep(pattern = "^MT-", x = rownames(x =GetAssayData(object = Data[[i]], slot = "counts")), value = TRUE)
  length(mito.genes)
  percent.mito <- Matrix::colSums(GetAssayData(object = Data[[i]], slot = "counts")[mito.genes, ]) / Matrix::colSums(GetAssayData(object = Data[[i]], slot = "counts"))*100
  
  Data[[i]]$mt_perc<-percent.mito}



dc_only <- merge(x = Data[[1]], y = Data[[2]], project="scrna")

pdf("nGene_uUMI_mito_diagnostics1.pdf",width=12,height=8)
VlnPlot(object = dc_only,
        features = c("nCount_RNA", "nFeature_RNA", "mt_perc"),ncol=3)

ggplot(dc_only@meta.data, aes(nCount_RNA, nFeature_RNA,color=type)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = c(200,7500), linetype = "dashed", colour = "red")

ggplot(dc_only@meta.data, aes(nCount_RNA, mt_perc,color=type)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0.15, linetype = "dashed", colour = "red")
dev.off()

obj.split=SplitObject(dc_only, split.by = "orig.ident")


# dc only -----------------------------------------------------------------




dc_only <-  obj.split[[1]]

pdf("nGene_uUMI_mito_diagnostics1.pdf",width=12,height=8)
VlnPlot(object = dc_only,
        features = c("nCount_RNA", "nFeature_RNA", "mt_perc"),ncol=3)

ggplot(dc_only@meta.data, aes(nCount_RNA, nFeature_RNA,color=type)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = c(200,7500), linetype = "dashed", colour = "red")

ggplot(dc_only@meta.data, aes(nCount_RNA, mt_perc,color=type)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0.15, linetype = "dashed", colour = "red")
dev.off()

DC <- subset(dc_only, subset = nCount_RNA<60000&nFeature_RNA > 200 & nFeature_RNA < 8500 &mt_perc < 30)

dc=DC%>%Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>% 
  ScaleData(verbose = FALSE)
###52338-->52182 left

C <- naive@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed1 <- order(apply(C, 1, median), decreasing = T)[20:1]

C2 <- dc1@assays$RNA@counts
C2 <- Matrix::t(Matrix::t(C2)/Matrix::colSums(C2)) * 100
most_expressed2 <- order(apply(C2, 1, median), decreasing = T)[20:1]


C3 <- dc@assays$RNA@counts
C3 <- Matrix::t(Matrix::t(C3)/Matrix::colSums(C3)) * 100
most_expressed3 <- order(apply(C3, 1, median), decreasing = T)[20:1]


pdf("top10.express.pdf")
par(mar = c(4, 8, 2, 1))
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell (naive)",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

boxplot(as.matrix(t(C2[most_expressed2, ])), cex = 0.1, las = 1, xlab = "% total count per cell(DC)",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

boxplot(as.matrix(t(C3[most_expressed3, ])), cex = 0.1, las = 1, xlab = "% total count per cell(DC2)",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()



#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pdf("highly_variable_gene_dc.pdf")
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1)
print(plot2)

dev.off()


#remove doublet




dc <- RunPCA(object = dc,
             pc.genes = dc@var.genes,
             do.print = TRUE,
             pcs.print = 1:5,
             genes.print = 5)

pdf("dc_PCA.pdf",width=12,height=8)
VizDimLoadings(dc, dims = 1:2, reduction = "pca")
DimPlot(dc, reduction = "pca")
DimHeatmap(dc, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(dc)
dev.off()

dc<- FindNeighbors(dc, dims = 1:30)
dc <- FindClusters(dc, resolution = 1.5)
dc <- RunUMAP(object = dc,
              dims =  1:30)




saveRDS(dc,"dc.umap.rds")



pdf("umap_dim30.dc.pdf")
#TSNEPlot(object = DataNorm, do.label = T,group.by="orig.ident")
p=DimPlot(object = dc, reduction = "umap")
print(p)
#ggsave(p,paste(samples[i],"_umap_preCCA.pdf",sep=""))
dev.off()



# Retrieve raw EPCAM counts
epcam_counts <- FetchData(
  dc,
  vars = "EPCAM",
  layer = "counts"
)[, "EPCAM"]

# Identify cells with at least one EPCAM UMI
epcam_positive_cells <- names(epcam_counts)[epcam_counts > 0]

length(epcam_positive_cells)

# Remove EPCAM-positive cells
dc.not.no_epcam <- subset(
  dc,
  cells = setdiff(Cells(dc), epcam_positive_cells)
)


library(Seurat)

# Make seurat_clusters the active identity
Idents(dc.not.no_epcam) <- "seurat_clusters"

# Cluster-to-cell-type mapping
cluster_names <- c(
  "0"  = "Tissue-resident regulatory cDC2",
  "1"  = "Migratory LAMP3int cDC2",
  "2"  = "Tissue-resident inflammatory cDC2",
  "3"  = "Tissue-resident inflammatory cDC2",
  "4"  = "Tissue-resident inflammatory cDC2",
  "5"  = "Tissue-resident regulatory cDC2",
  "6"  = "Tissue-resident inflammatory cDC2",
  "7"  = "Tissue-resident inflammatory cDC2",
  "8"  = "Tissue-resident inflammatory cDC2",
  "9"  = "Tissue-resident inflammatory cDC2",
  "10" = "Migratory LAMP3hi cDC2",
  "11" = "Migratory LAMP3hi cDC2",
  "12" = "Tissue-resident inflammatory cDC2",
  "13" = "Tissue-resident inflammatory cDC2",
  "14" = "Migratory LAMP3hi cDC2",
  "15" = "Migratory LAMP3hi cDC2",
  "16" = "Tissue-resident inflammatory cDC2",
  "17" = "Migratory LAMP3lo cDC2",
  "18" = "Migratory LAMP3hi cDC2",
  "19" = "Tissue-resident inflammatory cDC2",
  "20" = "Migratory inflammatory cDC2",
  "21" = "Migratory LAMP3hi cDC2",
  "22" = "Migratory inflammatory cDC2",
  "23" = "Unassigned",
  "24" = "Unassigned",
  "25" = "Migratory LAMP3lo cDC2",
  "26" = "Migratory inflammatory cDC2"
)

# Confirm all active identities have a proposed label


# Confirm all active identities have a proposed label
setdiff(levels(Idents(dc.not.no_epcam)), names(cluster_names))

# Rename active identities
dc.not.no_epcam <- RenameIdents(dc.not.no_epcam, cluster_names)

# Save renamed identities as a metadata column
dc.not.no_epcam$celltype_cluster <- as.character(Idents(dc.not.no_epcam))

# Make the new labels the active identity
Idents(dc.not.no_epcam) <- "celltype_cluster"

saveRDS(dc.not.no_epcam,"dc.singlet.rds")


# dc_enteorid -------------------------------------------------------------


DC_ent=obj.split[[2]]
DC_ent <- subset(DC_ent, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & mt_perc < 5)
#FilterCells(DC,subset.names = c("nFeature_RNA", "mt_perc"),low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.15))
DC_ent<- NormalizeData(object = DC_ent, verbose = FALSE)
#DC <- NormalizeData(DC)
DC_ent <- ScaleData(DC_ent, display.progress = T)
DC_ent <- FindVariableFeatures(DC_ent,selection.method = "vst", nfeatures = 3000)

#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pdf("highl_variable_gene_DC.pdf")
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(DC_ent), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(DC_ent)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1)
print(plot2)

dev.off()


DC_ent <- RunPCA(object = DC_ent,
             pc.genes = DC@var.genes,
             do.print = TRUE,
             pcs.print = 1:5,
             genes.print = 5)

pdf("DC_ent_PCA.pdf",width=12,height=8)
VizDimLoadings(DC_ent, dims = 1:2, reduction = "pca")
DimPlot(DC_ent, reduction = "pca")
DimHeatmap(DC_ent, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(DC_ent)
dev.off()

DC_ent<- FindNeighbors(DC_ent, dims = 1:30)
DC_ent <- FindClusters(DC_ent, renaivelution = 1.5)
DC_ent <- RunUMAP(object = DC_ent,
              dims =  1:30)


library(DoubletFinder)

sweep.res.list <- paramSweep(DC_ent, PCs = 1:20, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

bcmvn <- find.pK(sweep.stats)

print(bcmvn)


library(DoubletFinder)



# choose pK with highest BCmetric
best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
print(best_pK)

# ----------------------------
# 2. Estimate expected doublet rate
# ----------------------------
# Rough rule:
# 10x data often ~0.5% to 1% per 1000 cells
nCells <- ncol(DC_ent)
doublet_rate <- 0.008 * (nCells / 1000)   # adjust for your platform
nExp <- round(doublet_rate * nCells)

print(paste("Cells:", nCells))
print(paste("Expected doublets:", nExp))

# homotypic adjustment
annotations <- Idents(DC_ent)
homotypic.prop <- modelHomotypic(annotations)
nExp.adj <- round(nExp * (1 - homotypic.prop))

print(paste("Homotypic proportion:", round(homotypic.prop, 3)))
print(paste("Adjusted expected doublets:", nExp.adj))



nExp.adj <- as.numeric(nExp.adj)
nExp.adj <- round(nExp.adj)

nExp.adj
class(nExp.adj)


DC_ent <- doubletFinder(
  DC_ent,
  PCs = 1:20,
  pN = 0.25,
  pK = best_pK,
  nExp = nExp.adj,
  reuse.pANN = NULL,
  sct = FALSE
)

# ----------------------------
# 3. Run DoubletFinder
# ----------------------------

# ----------------------------
# 4. Find result columns
# ----------------------------

df_col <- grep("^DF.classifications", colnames(DC_ent@meta.data), value = TRUE)[1]
table(DC_ent@meta.data[[df_col]])

# optional: store in clean metadata column
DC_ent$doublet_status <- DC_ent@meta.data[[df_col]]

# ----------------------------
# 5. Visualize
# ----------------------------

pdf("DC_ent.double.pdf")
DimPlot(DC_ent, group.by = "doublet_status", cols = c("grey70", "red"))
dev.off()


# inspect doublet proportion by cluster
table(Idents(DC_ent), DC_ent$doublet_status)

# keep singlets only if wanted
DC_ent_singlet <- subset(DC_ent, subset = doublet_status == "Singlet")

saveRDS(DC_ent_singlet,"DC_ent_singlet.rds")



#annotate celltype
