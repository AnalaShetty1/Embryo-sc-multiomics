# Load DengQ reference file and parse data
deng_q = read.table('../references/RNA/GSE45719/RPKM.txt', header = T)

deng_q_zygote = data.frame('Gene_symbol' = deng_q$Gene_symbol, 'RPKM' = rowMeans(deng_q[,grepl('zy', colnames(deng_q))]))
deng_q_4cell = data.frame('Gene_symbol' = deng_q$Gene_symbol, 'RPKM' = rowMeans(deng_q[,grepl('4cell', colnames(deng_q))]))
deng_q_early2cell = data.frame('Gene_symbol' = deng_q$Gene_symbol, 'RPKM' = rowMeans(deng_q[,grepl('early2cell', colnames(deng_q))]))
deng_q_mid2cell= data.frame('Gene_symbol' = deng_q$Gene_symbol, 'RPKM' = rowMeans(deng_q[,grepl('mid2cell', colnames(deng_q))]))
deng_q_late2cell = data.frame('Gene_symbol' = deng_q$Gene_symbol, 'RPKM' = rowMeans(deng_q[,grepl('late2cell', colnames(deng_q))]))

# Load Wu refence and parse
wu_j = read.table('../references/RNA/GSE66582_stage_FPKM.txt', header = T)[,c(1,3,4,5,6)]
colnames(wu_j) = c('Gene_symbol', 'Wu_J_zygote', 'Wu_J_early2cell', 'Wu_J_2cell', 'Wu_J_4cell')

# Create data object and map gene names
fpkm <- rownames_to_column(as.data.frame(cnt.norm), 'gene_id')
fpkm$gene_name = merge(fpkm, genes,
                     by.x='gene_id', 
                     by.y='ensembl_gene_id', sort=F)$mgi_symbol

# Select mean expression per group
fpkm_zygote = data.frame('Gene_symbol' = fpkm$gene_name, 'RPKM' = rowMeans(fpkm[,grepl('Zygote', colnames(fpkm))]))
fpkm_2cell = data.frame('Gene_symbol' = fpkm$gene_name, 'RPKM' = rowMeans(fpkm[,grepl('2cell', colnames(fpkm))]))
fpkm_4cell = data.frame('Gene_symbol' = fpkm$gene_name, 'RPKM' = rowMeans(fpkm[,grepl('4cell', colnames(fpkm))]))

# Create zygote object
df_zygote = merge(fpkm_zygote, deng_q_zygote, by.x='Gene_symbol', by.y='Gene_symbol', sort=F) %>% 
  merge(wu_j) %>% select(Zygote=RPKM.x, Deng_Q_zygote=RPKM.y, Wu_J_zygote)

# Create 2-cell object
deng_q_2cell = data.frame('Gene_symbol' = deng_q_early2cell$Gene_symbol, 
                          'Deng_Q_early_2cell' = deng_q_early2cell$RPKM,
                          'Deng_Q_mid_2cell' = deng_q_mid2cell$RPKM,
                          'Deng_Q_late_2cell' = deng_q_late2cell$RPKM)

df_2cell = merge(fpkm_2cell, deng_q_2cell, by.x='Gene_symbol', by.y='Gene_symbol', sort=F) %>%
  merge(wu_j) %>% select(FPKM_2cell=RPKM, Deng_Q_early_2cell, Deng_Q_mid_2cell,
                         Deng_Q_late_2cell, Wu_J_early2cell, Wu_J_2cell)

# Create 4-cell object
df_4cell = merge(fpkm_4cell, deng_q_4cell, by.x='Gene_symbol', by.y='Gene_symbol', sort=F) %>% 
  merge(wu_j) %>% select(FPKM_4cell=RPKM.x, Deng_Q_4cell=RPKM.y, Wu_J_4cell)

# Compute and save zygote correlation
correlation = cor(df_zygote)
write_tsv(as.data.frame(correlation), '../plots/reference_correlation/zygote_correlation.tsv')

plot = ggcorrplot::ggcorrplot(correlation)
ggsave('../plots/reference_correlation/zygote_heatmap.png', plot)
ggsave('../plots/reference_correlation/zygote_heatmap.pdf', plot)

# Compute and save 2-cell correlation
correlation = cor(df_2cell)
write_tsv(as.data.frame(correlation), '../plots/reference_correlation/2cell_correlation.tsv')

plot = ggcorrplot::ggcorrplot(correlation)
ggsave('../plots/reference_correlation/2cell_heatmap.png', plot)
ggsave('../plots/reference_correlation/2cell_heatmap.pdf', plot)

# Compute and save 4-cell correlation
correlation = cor(df_4cell)
write_tsv(as.data.frame(correlation), '../plots/reference_correlation/4cell_correlation.tsv')

plot = ggcorrplot::ggcorrplot(correlation)
ggsave('../plots/reference_correlation/4cell_heatmap.png', plot)
ggsave('../plots/reference_correlation/4cell_heatmap.pdf', plot)

write_tsv(as.data.frame(fpkm), '../plots/reference_correlation/FPKM_data.tsv')
