# select post ZGA genes

ref = read.table('../references/GSE66582_stage_FPKM.txt')
oocyte_genes = data.frame(ref[ref$V2 <= 0.5, 1])
colnames(oocyte_genes) <- c('gene_name')

# Map to Ensembl
oocyte_genes = merge(oocyte_genes, genes,
                     by.x='gene_name', 
                     by.y='mgi_symbol', sort=F)$ensembl_gene_id

# This part of the code is ran per-group
################################################################################
cell = 'zygote'

# Load pseudobulk RT per group
rt <- read.table(paste0('../RT_times_per_gene_', cell, '.tsv'), row.names = 1, header = 1)
rt <- rownames_to_column(rt, 'gene_id')

# Subset RNA per group
subset = cnt.norm[,startsWith(colnames(cnt.norm), str_to_title(cell))]

# Binarize expression ad set as expressed genes with FPKM > 5 in at least 1 cell
genes_expressed = data.frame(is_exp = rowSums(log2(cnt.norm) > 5) > 1,
                             exp = rowSums(subset))

genes_expressed = rownames_to_column(genes_expressed, 'gene_id')

# Merge RT and expression
zyg = merge(genes_expressed, rt, by.x='gene_id', by.y='gene_id', sort=FALSE) 

# Select post ZGA genes
zyg = zyg[zyg$gene_id %in% oocyte_genes,]

# Group genes in groups of 50 ordered by RT
# and compute probability of expression per group
barplot_data <- zyg %>% arrange(RT) %>%
  group_by(x = ceiling(row_number()/50)) %>%
  summarise(prob = sum(is_exp)/50, RT=mean(RT))

# Plot RT and probability of expression
rt_exp_plot <- ggplot(barplot_data, aes(x=RT, y=prob)) +
  geom_area(fill='darkgray') +
  ylab('Probability of expression per bin (log FPKM > 1 for > 75% of cells)') +
  xlab('Mean RT per bin (rounded)') +
  labs(title=cell) +
  theme_classic() +
  geom_smooth(method='glm', color='black', linewidth=0.8)

ggsave(paste0('heatmaps/new_heatmaps/', cell, '_expression_vs_RT.png'), rt_exp_plot)
ggsave(paste0('heatmaps/new_heatmaps/', cell, '_expression_vs_RT.pdf'), rt_exp_plot)

# Binarise genes into early and late RT
zyg$category <- ifelse(zyg$RT > 0.5, 'Early', 'Late')

# Plot boxplot of early and late RT
rt_exp_boxplot <- ggplot(zyg, aes(x=category, y=exp, fill=category)) +
  geom_boxplot(outlier.size=2, weight=2) +
  scale_y_continuous(trans='log10') +
  scale_fill_brewer(palette = 'Set1') +
  ylab('Gene Expression (FPKM)') +
  xlab('RT category (early is > 0.5)') +
  labs(title=paste0(cell, '')) +
  stat_compare_means(label.x.npc = 'center') 

ggsave(paste0('plots/RT_RNA_correlation/', cell, '_late_early_boxplot.png'), rt_exp_boxplot)
ggsave(paste0('plots/RT_RNA_correlation/', cell, '_late_early_boxplot.pdf'), rt_exp_boxplot)

# Compare pseudotime and replication percentage
###################################################

# Load Kronos replication percentage per cell and parse the data
r_percent <- read.table('per_group/ReplicationPercent_all_cells.tsv', row.names = 1, header = 1)
split_cell_names = sapply(rownames(r_percent), function(x) strsplit(x, "_")[[1]], USE.NAMES=FALSE)
r_percent$cell = paste(split_cell_names[1,], split_cell_names[5,], sep= '_')
r_percent$group = split_cell_names[1,]

# Load pseudotime from the Monocle3 object in script_scRNA
pseudotime <- data.frame(pseudotime(cds))
pseudotime <- rownames_to_column(pseudotime, 'cell')
df <- merge(r_percent, pseudotime, by.x='cell', by.y='cell', sort=FALSE)
colnames(df) <- c('cell', 'percent_replication', 'group', 'pseudotime')

# Plot pseudotime and replication % per cell
corr <- ggplot(df, aes(x=pseudotime, y=percent_replication, color=group)) +
  geom_point(size=3) +
  theme_light() +
  scale_color_brewer(palette="Paired") +
  ylab('Percent Replication') +
  xlab('Pseudotime') +
  labs(title = paste('Correlation:', round(cor(df$pseudotime, df$percent_replication), 2))) 

ggsave('kronos_output/per_group/pseudotime_RT_correlation.png', corr)
ggsave('kronos_output/per_group/pseudotime_RT_correlation.pdf', corr)

