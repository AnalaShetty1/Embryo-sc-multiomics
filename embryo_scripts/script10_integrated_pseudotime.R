# Combine S-phase cells between groups and extract CN_bg values per bin per cell

mat_data_sphase <- rbind(SingleCell_zygote$SPhase, SingleCell_2cell$SPhase, SingleCell_4cell$SPhase) %>%
  mutate(bin = paste0(chr, ":", start, "-", end)) %>%
  dplyr::select(Cell, CN_bg, bin) %>%
  mutate(CN_bg = as.numeric(CN_bg)) %>%
  pivot_wider(names_from = Cell, values_from = CN_bg)

# Fill in NAs
mat_data_sphase[is.na(mat_data_sphase)] <- 0

# Rename cells
newnms <- paste0(sapply((strsplit(colnames(mat_data_sphase), split = "_")), '[', 1), 
                 "_", sapply((strsplit(colnames(mat_data_sphase), split = "_")), '[', 5))
newnms[1] = 'bin'
colnames(mat_data_sphase) = newnms

# Select only cells present in both RT and RNA
intersection = Reduce(intersect, list(colnames(mat_data_sphase), colnames(rnacnt)))

# Create RT assay
RT_assay = mat_data_sphase[, intersection]

# Subset RNA
rnacnt_subset = rnacnt[, intersection]
cnt.norm.1p_subset = cnt.norm.1p[, intersection]

# Create Seurat object and Add counts
sce_integrated <- CreateSeuratObject(counts=rnacnt_subset, names.field=1)

# Add normalized slot
sce_integrated <- SetAssayData(object=sce_integrated, slot="data", new.data=cnt.norm.1p_subset)
sce_integrated$nUMI = sce_integrated$nCount_RNA
sce_integrated$nGene = sce_integrated$nFeature_RNA   

# Add RT assay
RT_assay <- CreateAssayObject(counts = RT_assay)
sce_integrated[["RT"]] <- RT_assay

# Filter cells with low number of genes/UMIs
filtered_sce <- subset(
  x = sce_integrated,
  subset =
    (nUMI >= 700) &
    (nGene >= 500)
)

# Set RNA as default assay
DefaultAssay(filtered_sce) = 'RNA'

# Run standard Seurat pipeline
norm_sce <- NormalizeData(object = filtered_sce)
norm_sce <- FindVariableFeatures(object = norm_sce, nfeatures = 200)
norm_sce <- ScaleData(object = norm_sce)
norm_sce <- RunPCA(object = norm_sce, reduction.name = 'pca.rna', reduction.key = 'rnaPCA_')

# Set RT as default assay
DefaultAssay(norm_sce) = 'RT'

# Run standard Seurat pipeline
norm_sce <- FindVariableFeatures(norm_sce, nfeatures = 1000)
norm_sce <- ScaleData(object = norm_sce)
norm_sce <- RunPCA(object = norm_sce, reduction.name = 'pca.RT', reduction.key = 'RT_PCA_')

# Run Weighted Nearest Neighbor Analysis on both RNA and RT PCA embeddings
norm_sce <- FindMultiModalNeighbors(norm_sce, reduction.list = list("pca.rna", "pca.RT"), dims.list = list(1:20, 1:20), knn.range = 50)

# Create weighted nearest neighbour UMAP
norm_sce <- RunUMAP(norm_sce, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Create cell data set using the WNN embeddings
cds <- as.cell_data_set(norm_sce, reductions = 'wnn.umap')
cds@int_colData@listData$reducedDims$UMAP <- norm_sce@reductions$wnn.umap@cell.embeddings

# Run trajectory inference
cds <- cluster_cells(cds, resolution = 0.2)
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds, reduction_method = 'UMAP')

# Plot trajectory
trajectory <- plot_cells(cds, 
                         color_cells_by = 'ident',
                         group_cells_by = 'cluster',
                         group_label_size = 5,
                         trajectory_graph_segment_size = 1,
                         trajectory_graph_color = 'black',
                         alpha=0.5,
                         cell_size = 1.5)

ggsave('../plots/integrated_trajectory/trajectory.png', trajectory)
ggsave('../plots/integrated_trajectory/trajectory.pdf', trajectory)

# Color cells by pseudotime
pseudo <- plot_cells(cds, 
                     color_cells_by = 'pseudotime',
                     group_cells_by = 'cluster',
                     group_label_size = 7,
                     trajectory_graph_segment_size = 1,
                     trajectory_graph_color = 'black',
                     alpha=1,
                     cell_size = 2)

ggsave('../plots/integrated_trajectory/pseudotime.png', pseudo)
ggsave('../plots/integrated_trajectory/pseudotime.pdf', pseudo)



