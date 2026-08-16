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

monodc=readRDS("dc.singlet.rds")

Idents(monodc)=monodc$celltype_cluster


cols_dc <- c(
  "Tissue−resident inflammatory cDC2" = "#C17C9B",
  "Migratory inflammatory cDC2" = "#9E9E9E",
  "Undefined DCs" = "black",
  "Migratory LAMP3hi cDC2" = "#C9D7E3"
  # "Tissue-resident regulatory cDC2" = "#6FA8DC",
  #"Migratory LAMP3lo cDC2" = "#6B3F4F",
  #"Migratory LAMP3int cDC2" = "#2F6DB3"
)

pdf("monodc.2026.pdf",6,5)
DimPlot_scCustom(monodc,colors_use = cols_dc,label=F)
dev.off()

monodc$cluster_annot.new3=Idents(monodc)

# 
# monodc.split.meta <- all.not@meta.data[
#   which(all.not@meta.data$type == "Monol.+DC"),]
# rownames( monodc.split.meta)=gsub("_1","",rownames( monodc.split.meta))
# rownames( monodc.split.meta)=gsub("_2","",rownames( monodc.split.meta))
# 
# match(rownames(monodc@meta.data),rownames(monodc.split.meta))
# 
# m1=monodc@meta.data
# m2=monodc.split.meta
# 
# m12=merge(m1,m2[,c(1,12)],by=0,all=T)


pdf("fig4b.gene.monoDC.pdf", width = 8, height = 8)

DefaultAssay(monodc) <- "RNA"

FeaturePlot_scCustom(
  monodc,
  features = c("CD1A","CD1C","SIRPA","HLA-DRA","CLEC10A","CD14"),
  min.cutoff = 0,
  max.cutoff = 2,
  keep.scale = "all"
)

dev.off()




monodc$cluster_annot.new3=monodc$celltype_cluster

# Find neighbors and clusters
Idents(monodc)=monodc$cluster_annot.new3
expr <- GetAssayData(monodc, layer = "data", assay = "RNA")  # normalized log data


cell_types <- monodc$cluster_annot.new3


# Filter expression data
exprf <- expr[Matrix::rowSums(expr) > 0, ]
exprf <- as.matrix(exprf)  # Ensure matrix format
library(dorothea)
# Load human regulons
regulons <- dorothea_hs  # Use human regulons
#regulons <- regulons[regulons$confidence %in% c("A", "B"), ]

# Convert to regulon object and add target field
regulons_list <- df2regulon(regulons)
regulons_list <- lapply(regulons_list, function(x) {
  x$target <- names(x$tfmode)
  return(x)
})

# Extract target genes
target_genes <- unique(unlist(lapply(regulons_list, function(x) names(x$tfmode))))
head(target_genes, 10)
length(target_genes)

# Find common genes (case-insensitive)
common_genes <- intersect(toupper(rownames(exprf)), toupper(target_genes))
length(common_genes)
head(common_genes)

# Update case if matches found
if (length(common_genes) > 0) {
  rownames(exprf) <- toupper(rownames(exprf))
  regulons_list <- lapply(regulons_list, function(x) {
    names(x$tfmode) <- toupper(names(x$tfmode))
    x$target <- toupper(names(x$tfmode))
    return(x)
  })
} else {
  stop("No common genes found. Check species or gene identifiers.")
}

# Subset exprf
exprf <- exprf[common_genes, ]
exprf <- as.matrix(exprf)  # Ensure matrix format

library(viper)
# Run VIPER
tf_activities <- viper(
  exprf,
  regulons_list,
  verbose = TRUE,
  minsize = 4,
  eset.filter = FALSE
)

# Inspect results
head(tf_activities)

saveRDS(tf_activities,"tf_monol.dc.rds")

tf_activities=readRDS("tf_monol.dc.rds")

unique_cell_types <- unique(cell_types)
cell_type_matrices <- list()
for (ct in unique_cell_types) {
  cell_indices <- which(cell_types == ct)
  cell_type_matrices[[ct]] <- exprf[, cell_indices, drop = FALSE]
}

# Run VIPER per cell type
tf_activities_list <- list()
for (ct in names(cell_type_matrices)) {
  cat("Running VIPER for cell type:", ct, "\n")
  expr_ct <- cell_type_matrices[[ct]]
  if (nrow(expr_ct) < 4 || ncol(expr_ct) == 0) {
    warning(paste("Skipping cell type", ct, ": insufficient genes or cells"))
    next
  }
  tf_activities_list[[ct]] <- viper(
    expr_ct,
    regulons_list,
    verbose = TRUE,
    minsize = 4,
    eset.filter = FALSE
  )
}

# Aggregate mean TF activities
mean_tf_activities <- lapply(tf_activities_list, function(x) {
  if (!is.null(x)) rowMeans(x, na.rm = TRUE) else NULL
})
mean_tf_matrix <- do.call(cbind, mean_tf_activities)
colnames(mean_tf_matrix) <- names(mean_tf_activities)
head(mean_tf_matrix)




tf_var <- apply(mean_tf_matrix, 1, var, na.rm = TRUE)
# Rank and select the top N most variable TFs
top_n <- 20   # <-- change this number as you need (e.g., top 20, top 50)
top_tfs <- names(sort(tf_var, decreasing = TRUE))[1:top_n]
other_gene=c("ID2", "IRF4", "KLF4", "ZBTB46", "MAFB", "IRF8", "PU.1", "NOTCH2", "ZEB")

gene=unique(c(top_tfs,other_gene))
# Subset the matrix
mean_tf_top <- mean_tf_matrix[which(rownames(mean_tf_matrix)%in%gene), ]

library(pheatmap)
library(pheatmap)

# Define color palette
breaks <- seq(-2, 2, length.out = 101)

colors <- colorRampPalette(c("blue", "grey", "red"))(100)

pdf("fig4c.tf.monoDC.topvar.pdf",6,7)
pheatmap(
  mean_tf_top,
  scale = "row",
  color = colors,
  breaks = breaks,
  main = paste("Top25", "Variable TF Activities Across Cell Types")
)
dev.off()



