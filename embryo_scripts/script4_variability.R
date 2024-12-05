# This part of the code is run for each group separately ########################
# Recompute pseudobulk for the variability computation
current_group='zygote'
SingleCell <- SingleCell_zygote

# Rebin to 20kb bins
Bins = Kronos.scRT::GenomeBinning(
  Chr_size = Chromsize,
  Chr_filter =paste0('', 1:19),
  size = 20000,
  Cores = 6
)

# Compute S-phase
SingleCell$SPhase = Kronos.scRT::Rebin(PerCell = SingleCell$PerCell,
                                       scCN = SingleCell$CNV,
                                       Bins = Bins, Sphase = T)

# Rebin the control cells from script1
rebin_G1G2 = Kronos.scRT::Rebin(PerCell = SingleCell_controls$PerCell,
                                scCN = SingleCell_controls$CNV,
                                Bins = Bins, 
                                Sphase = F)

# Compute background
rebin_MedianG1G2 = Kronos.scRT::BackGround(G1_scCN = rebin_G1G2)

# Use G1 controls from script 1 to compute RT 
SingleCell$SPhase = Kronos.scRT::Replication_state(
  Samples = SingleCell$SPhase,
  background = rebin_MedianG1G2 %>%
    mutate(basename=current_group) %>% 
    mutate(group=current_group),
  Chr_filter = paste0('', 1:19),
  cores = 6
)

# Compute pseudobulk
SingleCell$pseudobulk=Kronos.scRT::pseudoBulkRT(S_scCN = SingleCell$SPhase)

# Compute variability
Var=Kronos.scRT::Variability(S_scCN=SingleCell$SPhase, scRT=SingleCell$pseudobulk)
Var=Kronos.scRT::TW_RTAnnotation(Variability=Var, RT_Groups = 5)

# Add gene annotation
annotation = read.table('annotation_organism_genes.bed')
colnames(annotation) <- c('chr', 'start', 'end', 'annotation')

Var = Kronos.scRT::TW_GenomeAnnotation(Variability = SingleCell$pseudobulk,
                                       GenomeAnnotation = annotation)

# Compute RT per gene and save data
rt_times_zygote <- Var %>% 
  group_by(category) %>%
  summarise(RT = mean(RT))

write.table(rt_times_zygote, paste0('../RT_times_per_gene_', current_group, '.tsv'), row.names = FALSE, quote = F)

# Compute percentage of replication
s_phase <- SingleCell$SPhase %>% group_by(Cell) %>%
  summarise(PercentageReplication = dplyr::first(PercentageReplication))
write.table(s_phase, paste0('Replication_Percent_', current_group, '.tsv'), row.names = FALSE, quote = F)

# Compute Twidth using 5 categories
Var$category <- Var$category.x
Fit_Data=Kronos.scRT::Twidth_fit_data(
  df = Var,
  ncores = 6
)

# Generate and save T-width barplot
Twidth = Kronos.scRT::Twidth(Fit_Data)

twidth_barplot = Kronos.scRT::Twidth_barplot(Variability = Var,Twidth = Twidth)
ggsave(paste0('twidth_plots/', current_group, '_twidth_barplot.png'), twidth_barplot)
ggsave(paste0('twidth_plots/', current_group, '_twidth_barplot.pdf'), twidth_barplot)
twidth_plot = Kronos.scRT::Twidth_extended_plot(Variability = Var,Fitted_data = Fit_Data,Twidth = Twidth)
ggsave(paste0('twidth_plots/', current_group, '_twidth.png'), twidth_plot)
ggsave(paste0('twidth_plots/', current_group, '_twidth.pdf'), twidth_plot)

# Generate and save bin probability plot
BinProbS=Kronos.scRT::Prepare_S_phase_cells_forBinRepProb(S = SingleCell$SPhase,RT = SingleCell$pseudobulk)
bin_rep_plot <- Kronos.scRT::BinRepProbPlot(Variability = BinProbS)
ggsave(paste0('twidth_plots/', current_group, '_bin_rep_prob.png'), bin_rep_plot)
ggsave(paste0('twidth_plots/', current_group, '_bin_rep_prob.pdf'), bin_rep_plot)

# Extract RT values per gene and per cell for integration
Var = Kronos.scRT::TW_GenomeAnnotation(Variability = SingleCell$SPhase,
                                       GenomeAnnotation = annotation)

Var = Var %>%
  filter(category != '_Unknown_') %>%
  group_by(Cell, category) %>%
  summarise(RT = mean(Rep)) %>%
  pivot_wider(names_from = Cell, values_from = RT)

split_cell_names = sapply(colnames(Var)[2:length(colnames(Var))], function(x) strsplit(x, "_")[[1]], USE.NAMES=FALSE)
colnames(Var) = c('gene_id', paste(split_cell_names[1,], split_cell_names[5,], sep= '_'))

write.table(Var, paste0('Rep_per_gene_per_cell_', current_group, '.tsv'), row.names = FALSE, quote = F)




