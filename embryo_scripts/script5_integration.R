# Set seed for heatmaps
set.seed(123)

# Select oocyte genes - post ZGA genes
ref = read.table('../references/GSE66582_stage_FPKM.txt')
oocyte_genes = data.frame(ref[ref$V2 <= 0.5, 1])
colnames(oocyte_genes) <- c('gene_name')

# Create RNA FPKM object
fpkm <- rownames_to_column(as.data.frame(cnt.norm), 'gene_id')

# Map gene names
fpkm$gene_name = merge(fpkm, genes,
                       by.x='gene_id', 
                       by.y='ensembl_gene_id', sort=F)$mgi_symbol

fpkm = fpkm[!duplicated(fpkm$gene_name),]

# Merge stages and select post ZGA genes
fpkm = fpkm %>%
  mutate(zyg_exp = rowMeans(fpkm[,grepl('Zygote', colnames(fpkm))])) %>%
  mutate(cell2_exp = rowMeans(fpkm[,grepl('2cell', colnames(fpkm))])) %>% 
  mutate(cell4_exp = rowMeans(fpkm[,grepl('4cell', colnames(fpkm))])) %>%
  filter(gene_name %in% oocyte_genes$gene_name) %>%
  remove_rownames() %>%
  column_to_rownames('gene_id') %>%
  select(zyg_exp, cell2_exp, cell4_exp)

# Remove non-expressed genes and set FPKM of 0 to -3
fpkm = fpkm[rowSums(fpkm) > 0,]
fpkm[fpkm == 0] = -3

# Create RT object across stages
zygote <- read.table(paste0('../RT_times_per_gene_zygote.tsv'), row.names = 1, header = 1)
zygote = rownames_to_column(zygote, 'Gene_ID')
colnames(zygote) = c('Gene_ID', 'Zygote')
zygote = zygote[zygote$Gene_ID != '_Unknown_',]

cell2 <- read.table(paste0('../RT_times_per_gene_2cell.tsv'), row.names = 1, header = 1)
cell2 = rownames_to_column(cell2, 'Gene_ID')
colnames(cell2) = c('Gene_ID', '2_cell')
cell2 = cell2[cell2$Gene_ID != '_Unknown_',]

cell4 <- read.table(paste0('../RT_times_per_gene_4cell.tsv'), row.names = 1, header = 1)
cell4 = rownames_to_column(cell4, 'Gene_ID')
colnames(cell4) = c('Gene_ID', '4_cell')
cell4 = cell4[cell4$Gene_ID != '_Unknown_',]

l <- list(zygote, cell2, cell4)

df_rt = purrr::reduce(.x = l, merge, by = c('Gene_ID'), all = T)

# Map to gene name
df_rt$gene_name = merge(df_rt, genes,
                       by.x='Gene_ID', 
                       by.y='ensembl_gene_id', all.x=T, sort=F)$mgi_symbol

df_rt = df_rt[!duplicated(df_rt$gene_name),]

# Select post ZGA genes
df_rt = df_rt %>%
 filter(gene_name %in% oocyte_genes$gene_name)

# Merge RNA and RT using both gene name and ID
fpkm = rownames_to_column(fpkm, 'Gene_ID')
df = inner_join(fpkm, df_rt)

# Remove NAs
df[is.na(df)] = 0

# Plot RT as heatmap and cluster the rows using 8 clusters
plot = Heatmap(as.matrix(df[,5:7]),
               cluster_columns = F,
               show_row_names = F,
               col = colorRamp2(c(0, 0.45, 0.5, 0.55, 1), c("blue", "blue", "black", "red", "red")),
               column_names_rot=45,
               heatmap_legend_param = list(title = "RT"),
               km=8)

png(file="../output/RT_clustered.png",  width = 480, height = 960, units = "px")
plot
dev.off()
dev.off()

pdf(file="../output/RT_clustered.pdf",  width = 8, height = 12)
plot
dev.off()

# Keep the row order of the clusters
ht = draw(plot)  
rcl.list <- row_order(ht)

# Start ATAC-seq analysis
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Load the ATAC peak files per group
# Annotate peaks and select Ensembl IDs
# Merge peaks by gene
# Filter outliers and rescale values from 0 to 1
zyg_atac <- annotatePeak('../references/ATAC-seq/2cell_early_peaks.bed', tssRegion=c(-5000, 5000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

zyg_atac = data.frame(zyg_atac@anno)[,c('ENSEMBL', 'V5')] %>% 
  group_by(ENSEMBL) %>%
  summarise('ATAC_early_2cell' = mean(V5)) %>% 
  rename('Gene_ID' = ENSEMBL) %>%
  filter(ATAC_early_2cell < min(boxplot(ATAC_early_2cell)$out)) %>% 
  mutate(ATAC_early_2cell=rescale(ATAC_early_2cell, to=c(0,1)))

cell2_atac <- annotatePeak('../references/ATAC-seq/2cell_peaks.bed', tssRegion=c(-5000, 5000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")

cell2_atac = data.frame(cell2_atac@anno)[,c('ENSEMBL', 'V5')] %>% 
  group_by(ENSEMBL) %>%
  summarise('ATAC_2cell' = mean(V5)) %>% 
  rename('Gene_ID' = ENSEMBL) %>%
  filter(ATAC_2cell < min(boxplot(ATAC_2cell)$out)) %>% 
  mutate(ATAC_2cell=rescale(ATAC_2cell, to=c(0,1)))

cell4_atac <- annotatePeak('../references/ATAC-seq/4cell_peaks.bed', tssRegion=c(-5000, 5000),
                           TxDb=txdb, annoDb="org.Mm.eg.db")

cell4_atac = data.frame(cell4_atac@anno)[,c('ENSEMBL', 'V5')] %>% 
  group_by(ENSEMBL) %>%
  summarise('ATAC_4cell' = mean(V5)) %>%
  rename('Gene_ID' = ENSEMBL) %>%
  filter(ATAC_4cell < min(boxplot(ATAC_4cell)$out)) %>% 
  mutate(ATAC_4cell=rescale(ATAC_4cell, to=c(0,1)))

# Merge ATAC-seq objects
l <- list(zyg_atac, cell2_atac, cell4_atac)
df_atac = purrr::reduce(.x = l, merge, by = c('Gene_ID'), all = T)

# Remove NAs
df_atac[is.na(df_atac)] = 0

# Map gene names
df_atac$gene_name = merge(df_atac, genes,
                        by.x='Gene_ID', 
                        by.y='ensembl_gene_id', all.x=T, sort=F)$mgi_symbol

df_atac = df_atac[!duplicated(df_atac$gene_name),]

# Select post ZGA genes
df_atac = df_atac %>%
  filter(gene_name %in% oocyte_genes$gene_name)

# Add ATAC-seq data to the RNA/RT object
df_final = left_join(df, df_atac, by='gene_name')
df_final$Gene_ID = df_final$Gene_ID.x
df_final$Gene_ID.x = NULL
df_final[is.na(df_final)] = 0

# Extract clusters
library(magrittr)

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(Gene_ID = df_final[rcl.list[[i]],]$Gene_ID,
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

# Add cluster IDs
df_final = inner_join(df_final, clu_df)

# Create expression object per cluster
exp = df_final %>%
  pivot_longer(c('zyg_exp', 'cell2_exp', 'cell4_exp'), names_to='cell', values_to='Expression') %>%
  mutate(cell = recode_factor(cell, "zyg_exp" = "zygote", 'cell2_exp' = '2cell', 'cell4_exp' = '4cell')) %>%
  group_by(Cluster, cell) %>%
  summarise(Expression=mean(Expression))

# Create RT object per cluster
rt = df_final %>%
  pivot_longer(c('Zygote', '2_cell', '4_cell'), names_to='cell', values_to='RT') %>%
  mutate(cell = recode_factor(cell, "Zygote" = "zygote", '2_cell' = '2cell', '4_cell' = '4cell')) %>%
  group_by(Cluster, cell) %>%
  summarise(RT=mean(RT))

# Create ATAC-seq object per cluster
atac = df_final %>%
  pivot_longer(c('ATAC_early_2cell', 'ATAC_2cell', 'ATAC_4cell'), names_to='cell', values_to='ATAC') %>%
  mutate(cell = recode_factor(cell, "ATAC_early_2cell" = "zygote", 'ATAC_2cell' = '2cell', 'ATAC_4cell' = '4cell')) %>%
  group_by(Cluster, cell) %>%
  summarise(ATAC=mean(ATAC))

# Combine cluster means
cluster_means = inner_join(exp, rt)
cluster_means = inner_join(cluster_means, atac)

# Use RNA scaling coefficient
coeff <- max(cluster_means$Expression)

# Generate pairwise plots for RNA/RT and ATAC-seq for the cluster means
per_cluster_comparison = ggplot(cluster_means, aes(x=cell, group=1)) +
  geom_point(aes(y=RT),size=2) +
  geom_point(aes(y=Expression / coeff),size=2) +
  geom_line( aes(y=RT, linetype='RT')) + 
  facet_wrap(~ Cluster) +
  geom_line(aes(y=Expression / coeff, linetype='Expression')) +
  scale_y_continuous(
    name = "RT",
    sec.axis = sec_axis(~.*coeff, name="RNA Expression")
  ) + scale_linetype_manual("Dataset",values=c("RT" = 1, "Expression" = 2)) +
  xlab('Group')

ggsave('../output/per_cluster_comparison_RT_RNA_inner_join.png', per_cluster_comparison)
ggsave('../output/per_cluster_comparison_RT_RNA_inner_join.pdf', per_cluster_comparison)

per_cluster_comparison = ggplot(cluster_means, aes(x=cell, group=1)) +
  geom_point(aes(y=ATAC),size=2) +
  geom_point(aes(y=Expression / coeff),size=2) +
  geom_line( aes(y=ATAC, linetype='ATAC')) + 
  facet_wrap(~ Cluster) +
  geom_line(aes(y=Expression / coeff, linetype='Expression')) +
  scale_y_continuous(
    name = "ATAC",
    sec.axis = sec_axis(~.*coeff, name="RNA Expression")
  ) + scale_linetype_manual("Dataset",values=c("ATAC" = 1, "Expression" = 2)) +
  xlab('Group')

ggsave('../output/per_cluster_comparison_ATAC_RNA_inner_join.png', per_cluster_comparison)
ggsave('../output/per_cluster_comparison_ATAC_RNA_inner_join.pdf', per_cluster_comparison)

per_cluster_comparison = ggplot(cluster_means, aes(x=cell, group=1)) +
  geom_point(aes(y=ATAC),size=2) +
  geom_point(aes(y=RT),size=2) +
  geom_line( aes(y=ATAC, linetype='ATAC')) + 
  geom_line(aes(y=RT, linetype='RT')) +
  facet_wrap(~ Cluster) +
  scale_linetype_manual("Dataset",values=c("ATAC" = 1, "RT" = 2)) +
  xlab('Group')

ggsave('../output/per_cluster_comparison_ATAC_RT_inner_join.png', per_cluster_comparison)
ggsave('../output/per_cluster_comparison_ATAC_RT_inner_join.pdf', per_cluster_comparison)

write.table(cluster_means, '../output/cluster_means.tsv', row.names = F)

# Compute GO term enrichment per cluster
for (i in seq(1, 8)){
  genes = clu_df[clu_df$Cluster == paste0('cluster', i),]$Gene_ID
  
  r_go <- clusterProfiler::enrichGO(
    gene = genes ,
    universe = clu_df$GeneID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENSEMBL",
    ont = 'ALL',
    pAdjustMethod = "BH",
    readable = TRUE,
    pool = TRUE
  )
  
  write.table(r_go@result, paste0('heatmaps/new_heatmaps/cluster_', i, '_GO.tsv'))
  if (!is.null(r_go)){
    ggsave(paste0('heatmaps/new_heatmaps/cluster_', i, '_GO.png'), barplot(r_go))
    ggsave(paste0('heatmaps/new_heatmaps/cluster_', i, '_GO.pdf'), barplot(r_go))
  }
}

write.table(df_final, 'heatmaps/new_heatmaps/data.tsv')
