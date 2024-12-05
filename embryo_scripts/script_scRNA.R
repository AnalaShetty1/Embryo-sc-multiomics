library(Seurat)
library(gdata)
library(tidyverse)
library(ggplot2)
library(scales)
library(viridis)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(clusterProfiler)
library(ggpubr)
library(org.Mm.eg.db)
library(SeuratWrappers)
library(monocle3)

setwd("/home/data/BISP-666")

# Load read counts
x <- read.table("zhan2142/rna_seq/Counts/subread_counts.txt", header=T, row.names=1, sep="\t", stringsAsFactors = F)
meta <- read.table("zhan2142/rnaseq.group.csv", header=T, row.names=1, sep=",", stringsAsFactors = F)

# Select gene counts and gene length columns
cnt <- x[,5:ncol(x)]

# Remove leading X
new_names = c()
for (c in colnames(cnt)){
  if (startsWith(c,'X')){
    c <- substr(c, 2, nchar(c))
  }
  new_names<- c(new_names, c)
}
colnames(cnt) <- new_names

# Gene length normalization and log1p transform
cnt.norm = cnt[,2:ncol(cnt)] / cnt[,1]
cnt.norm = t(t(cnt.norm) * 1e6 / colSums(cnt.norm))
cnt.norm.1p <- log1p(cnt.norm)

# Remove genes with 0 counts across samples
keep <- rowSums(cnt.norm) > 0

# Remove gene length column in count data and subset genes
rnacnt <- cnt[keep, 2:ncol(cnt)]

# Reorder columns
rnacnt <- cnt[, rownames(meta)]
cnt.norm.1p <- cnt.norm.1p[, rownames(meta)]
cnt.norm <- cnt.norm[, rownames(meta)]

# Rename all objects
newnms <- paste0(sapply((strsplit(colnames(rnacnt), split = "_")), '[', 1), "_", sapply((strsplit(colnames(rnacnt), split = "_")), '[', 5))
rownames(meta) <- newnms
colnames(rnacnt) <- newnms
colnames(cnt.norm.1p) <- newnms
colnames(cnt.norm) <- newnms

# Load gene mapping table
genes <- read.table("zhan2142/gene_id_gene_name_map.txt", stringsAsFactors = FALSE)
names(genes) <- c("ensembl_gene_id", "mgi_symbol")
mgisymbol <- genes$mgi_symbol
names(mgisymbol) <- genes$ensembl_gene_id

# Select mitochondrial genes
mitochondrial_genes <- genes[substr(genes$mgi_symbol, 1, 3) == 'mt-', ]$ensembl_gene_id

# Add counts
sce <- CreateSeuratObject(counts=rnacnt, names.field=1)

# Add normalized slot
sce <- SetAssayData(object=sce, slot="data", new.data=cnt.norm.1p)

# Rename meatadata fields
sce$stage <- meta$Group
sce$nUMI = sce$nCount_RNA
sce$nGene = sce$nFeature_RNA   
sce$mito_ratio <-  PercentageFeatureSet(object = sce, features = mitochondrial_genes) / 100

# Remove cells with low genes/UMIs
filtered_sce <- subset(
  x = sce,
  subset =
    (nUMI >= 700) &
    (nGene >= 500)
)

# Process with the standard Seurat pipeline
norm_sce <- NormalizeData(object = filtered_sce)
norm_sce <- FindVariableFeatures(object = norm_sce)
norm_sce <- ScaleData(object = norm_sce)
norm_sce <- RunPCA(object = norm_sce)

# Select maximal number of PCs
pct <- norm_sce[["pca"]]@stdev / sum(norm_sce[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
max_pc <- min(co1, co2)
norm_sce$informative_pcs <- max_pc

# UMAP, find neighbors and clustering with that number of PCs
norm_sce <- RunUMAP(object = norm_sce, dims = 1:max_pc)
norm_sce <- FindNeighbors(object = norm_sce, reduction = 'umap', dims = 1:2)
norm_sce <- FindClusters(object = norm_sce, resolution = 0.2)
norm_sce <- RunTSNE(object = norm_sce, dims = 1:max_pc)

# Plot UMAP
umap_plot <- DimPlot(norm_sce,
                     reduction = 'umap', 
                     group.by = 'stage', pt.size = 3) +
  theme_minimal() +
  scale_color_brewer(palette="Paired") +
  labs(title='UMAP plot')

ggsave('plots/umap.png', plot = umap_plot, dpi = 400 , bg = "white")
ggsave('plots/umap.pdf', plot = umap_plot, dpi = 400 , bg = "white")

# Compute DE genes between stages pairwise
Idents(norm_sce) <- 'stage'

# Create this object for all combinations between groups
zygote_2_cell <- FindMarkers(norm_sce, ident.1 = 'Zygote', ident.2 = '4cell', test.use = 'DESeq2')
cells.subset <- subset(norm_sce, subset = (stage %in% c('4cell', 'Zygote')))
significant <- zygote_2_cell[zygote_2_cell$p_val_adj < 0.05,]

# Generate heatmap for the significant DE genes
heatmap <- DoHeatmap(cells.subset, features = rownames(significant), label = F) +
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  theme(axis.text.y = element_text(size=0)) +
  labs(title = 'Zygote vs 4 cell') +
  scale_fill_viridis(na.value='white') 

ggsave('plots/Zygote_4_cell_heatmap.png', heatmap)
ggsave('plots/Zygote_4_cell_heatmap.pdf', heatmap)

# GO-term enrichment on the significant DE genes
r_go <- clusterProfiler::enrichGO(
  gene = rownames(significant),
  universe = rownames(zygote_2_cell),
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = 'ALL',
  pAdjustMethod = "BH",
  readable = TRUE,
  pool = TRUE
)

write.table(r_go, 'plots/GO_all.tsv')
ggsave('plots/GO_all.png', barplot(r_go))
ggsave('plots/GO_all.pdf', barplot(r_go))

################# Trajectory inference #########################################

# Create CDS oject
norm_sce@active.assay = 'RNA'
cds <- as.cell_data_set(norm_sce)
cds <- preprocess_cds(cds, num_dim = 24)
cds <- reduce_dimension(cds, reduction_method = 'UMAP')

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

cds <- learn_graph(cds, use_partition = FALSE)

## Step 6: Order cells
cds <- order_cells(cds, reduction_method = 'UMAP')

# Plot the trajectory
trajectory <- plot_cells(cds, 
                         color_cells_by = 'stage',
                         group_cells_by = 'cluster',
                         group_label_size = 7,
                         trajectory_graph_segment_size = 1,
                         trajectory_graph_color = 'black',
                         alpha=0.5,
                         cell_size = 1.5)

ggsave('plots/trajectory.png', trajectory)
ggsave('plots/trajectory.pdf', trajectory)
