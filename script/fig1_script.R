library(Seurat)
#Seurat_4.0.2 
#library(SeuratData)
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


dc.not=readRDS("input_fig1.rds")
library(scCustomize)
celltype <- levels(Idents(dc.not))

cols_use <- c(
  "#A76E93",
  "#5D3A4D",
  "#3B71AF",
  "#79AAD2",
  "#CBDCEB",
  "#A8AAAA"
)

names(cols_use) <- celltype

pdf("fig1a.umap.dc.noT.pdf",8,5)
library(scCustomize)
library(scCustomize)

DimPlot_scCustom(
  seurat_object = dc.not,
  colors_use = cols_use,
  label = FALSE
)
dev.off()


pdf("fig1b.gene.DC.pdf", width = 12, height = 6)

DefaultAssay(dc.not) <- "RNA"

p=FeaturePlot_scCustom(
  dc.not,
  features = c("CD1A","CD1C","SIRPA","HLA-DRA","CLEC10A","CD14"),
  min.cutoff = 0,
  max.cutoff = 2,
  keep.scale = "all"
)
p+ patchwork::plot_layout(ncol = 3)
dev.off()



# fig1e -------------------------------------------------------------------



DefaultAssay(dc.not)="RNA"

expr <- GetAssayData(dc.not, layer = "data", assay = "RNA")  # normalized log data

# If using a Seurat object
library(Seurat)
cell_types <- dc.not$cluster_annot.new3


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

pdf("fig1e.tf.DC.topvar.top25.pdf",6,7)
pheatmap(
  mean_tf_top,
  scale = "row",
  color = colors,
  breaks = breaks,
  main = paste("Top25", "Variable TF Activities Across Cell Types")
)
dev.off()




# fig1i -------------------------------------------------------------------



cytokin=c("IL1B", "IL6", "IL10", "IL12A","IL12B", 
          "IL23A", "CSF1", "CXCL10", "CCL2", "CCL3",
          "CCL4", "CCL13", "CCL17", "CCL22")

genes_use <- intersect(cytokin, rownames(dc.not))
setdiff(cytokin, genes_use)   # check missing genes

# Option 1: scale these genes, then heatmap

DefaultAssay(dc.not) <- "RNA"


library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)


cols_use <- c(
  "#A76E93",
  "#5D3A4D",
  "#3B71AF",
  "#79AAD2",
  "#CBDCEB",
  "#A8AAAA"
)

cols_dc <- c(
  "Tissue−resident inflammatory cDC2" = "#C17C9B",
  "Migratory inflammatory cDC2" = "#9E9E9E",
  #"Undefined DCs" = "black",
  "Migratory LAMP3hi cDC2" = "#C9D7E3",
  "Tissue-resident regulatory cDC2" = "#6FA8DC",
  "Migratory LAMP3lo cDC2" = "#6B3F4F",
  "Migratory LAMP3int cDC2" = "#2F6DB3"
)


levels(Idents(dc.not)) <- gsub("\u2212", "-", levels(Idents(dc.not)))

pdf("dotplot_cytokin.pdf",10,6)
p=DotPlot_scCustom(seurat_object = dc.not, features = genes_use)+coord_flip()
p+theme(
  axis.text.y = element_text(size = 8),
  axis.text.x = element_text(size = 8))+theme(
    axis.text.x = element_text(angle = 45, hjust = 1))



dev.off()

pdf("figs1b.marker.featureplot.pdf", width = 8, height = 8)

DefaultAssay(dc.not) <- "RNA"

FeaturePlot_scCustom(
  dc.not,
  features = c("CD40","CD80","CD86","CD83","CCR7","LAMP3","S100A9","FCGR3A","ITGAX"),
  #  min.cutoff = 0,
  # max.cutoff = 2,
  keep.scale = "all"
)

dev.off()



# figs1 -------------------------------------------------------------------




dc.not=readRDS("dc.input.fig1.rds")

library(scCustomize);library(ggplot2)

gene.dc2 <- c("CLEC10A", "SIGLEC10", "CD1A", "CD33", "CLEC4G", "S100B", "IL1R2", 
              "SPP1", "TNFRSF11A", "RAB33A", "CD40", "HLA-DRB1", "HLA-DRA", 
              "CSF1", "CCL2", "CCL3", "CCL4", "CCL4L2", "TNFSF13B", "CD80","CD83","CD86","LAMP3","CCR7",
              "ITGAX","ITGAM",
              "FABP5", "CCL22", "LYZ","BACH2", "INPP4B", "CAMK4", "SKAP1", "THEMIS", "ANK3")

gene.dc3=c("CD81", "S100A9", "S100A8", "MRC1", "CD1C", "CD1B", "CD1A", "ANXA1", 
           "CCL22", "IL18", "APOC1",  "FABP4", "MMP9", "SOD2", "UPP1", 
           "RALA", "EMP3", "JUN", "KLF6", "ATF3", "CXCL8")

pdf("fig.s1.gene.dotplot.pdf",15,5)
DotPlot_scCustom(dc.not,features=unique(c(gene.dc2,gene.dc3)))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#+ coord_flip()

dev.off()


gene_to_plot=c(
               "CD40","CD80","CD86","CD83","CCR7","LAMP3","S100A9","FCGR3A","ITGAX")


pdf("figs1b.gene.featureplot.pdf",15,5)
DefaultAssay(dc.not)="RNA"

#DimPlot(obj1,label=T)
FeaturePlot_scCustom(
  dc.not,
  features = gene_to_plot,
  # cols = c("lightyellow", "red"), 
  order = TRUE,combine=T
  #ncol = 4   # place it here as its own argument, not in ...
) 

#wrap_plots(p.list, ncol = 4) & theme(legend.position = "right")

dev.off()

