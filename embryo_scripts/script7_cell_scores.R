# Load annotation
annotation = read.table('kronos_output/annotation_organism_genes.bed')
colnames(annotation) <- c('chr', 'start', 'end', 'annotation')

# Run analysis for all groups
for (group in c('zygote', '2cell', '4cell')){
  
  #  Generate object for cell score
  cell_scores = data.frame()
  
  # Load replication values per cell
  file = read.table(paste0('Rep_per_gene_per_cell_', group, '.tsv'), header = T)
  colnames(file) = gsub('X', '', colnames(file))
  
  # Loop across cells
  for (cell in colnames(file[,2:ncol(file)])){
    
    # Extract per-cell RT
    rt = data.frame('Gene_ID' = file$gene_id, 'RT' = file[,cell])
    rt = rt %>% drop_na()
    
    # Map cell name to RNA
    rna_name = ''
    for (rna in colnames(cnt.norm)){
      if(grepl(cell, rna)){
        rna_name = rna
      }
    }
    
    if(rna_name != ''){
      
      # Select RNA for cell
      rna = data.frame('Expression' = cnt.norm[,rna_name])
      rna = rownames_to_column(rna, 'Gene_ID')
      
      # Join RNA and RT
      df = inner_join(rna, as.data.frame(rt))
      
      # Select genes with FPKM > 1 as expressed
      # Bin genes in bins of 50 and compute
      # probability of expression
      df$is_exp = log2(df$Expression) > 1
      df = df %>% arrange(RT) %>%
        group_by(x = ceiling(row_number()/50)) %>%
        summarise(prob = sum(is_exp)/50, RT=mean(RT))
      
      # Compute slope of regressiion between 
      # probability of expression and RT
      mod <- lm(df$prob ~ df$RT)
      cf <- coef(mod)
      
      # Save slope
      cell_scores = rbind(cell_scores, c(cell, cf[[2]], mean(df$RT)))
    }
  }
  
  # Modify cell scores
  colnames(cell_scores) = c('Cell', 'Slope', 'RT')
  cell_scores$Slope = as.numeric(cell_scores$Slope)
  cell_scores$RT = as.numeric(cell_scores$RT)
  
  write.table(cell_scores, paste0('output/cell_scores_embryos_', group, '.tsv'), quote=F, row.names = F)
  
  # Create cell score plot
  plot = ggplot(cell_scores, aes(x = RT, y=Slope, color=RT)) +
    geom_jitter(size=2.5) +
    xlab('RT') +
    scale_colour_gradient(low = 'red',high = 'blue') +
    labs(title = 'Correlation slope per cell vs RT') +
    geom_smooth(method='lm', color='black', linewidth=0.8)
  
  ggsave(paste0('output/correlation_embryos_', group, '.png'), plot)
  ggsave(paste0('output/correlation_embryos_', group, '.pdf'), plot)
  
  # Create RT/RNA correlation plot
  plot = ggplot(cell_scores, aes(x = tools::toTitleCase(group), y=Slope, color=RT)) +
    geom_jitter(size=2.5) +
    xlab('') +
    scale_colour_gradient(low = 'red',high = 'blue') +
    labs(title = 'RT/RNA correlation line slope per cell')
  
  ggsave(paste0('output/cellscore_embryos_', group, '.png'), plot)
  ggsave(paste0('output/cellscore_embryos_', group, '.pdf'), plot)
  
  # Pseudobulk score computation
  rt = read.table(paste0('RT_times_per_gene_', group, '.tsv'), header = T)
  
  # Select cells per group
  subset = cnt.norm[,startsWith(colnames(cnt.norm), str_to_title(group))]
  
  # Binarize expression: FPKM > 1 in at least 1 sample
  genes_expressed = data.frame(is_exp = rowSums(log2(subset) > 1) > 1,
                               exp = rowMeans(subset))
  genes_expressed = rownames_to_column(genes_expressed, 'category')
  
  # Join RNA and RT on gene ID
  zyg = inner_join(genes_expressed, rt)
  
  # Bin genes in groups of 50 and compute probability of expression
  barplot_data <- zyg %>% arrange(RT) %>%
    group_by(x = ceiling(row_number()/50)) %>%
    summarise(prob = sum(is_exp)/50, RT=mean(RT))
  
  # Compute correlation line
  mod <- lm(barplot_data$prob ~ barplot_data$RT)
  cf_pseudobulk <- coef(mod)
  
  # Print slope of correlation line of pseudobulk
  print(cf_pseudobulk[[2]])
}
