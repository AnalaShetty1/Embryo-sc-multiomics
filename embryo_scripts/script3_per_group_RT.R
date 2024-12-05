library(zoo)

# This part of the code is run for each group separately ########################
current_group='4cell'

# load fixed binning function
bins = binning(
  RefGenome = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107.fa',
  bowtie2_index = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107',
  black_list = '/home/data/BISP-666/zhan2142/src/mm10-blacklist.v2.bed',
  directory_to_bamfiles = paste0('/home/data/BISP-666/kronos_output/BAM/', current_group, '/'),
  cores = 10,
  bin_size = 20000
)

# Create object
SingleCell = CallCNV(
  directory = paste0('/home/data/BISP-666/kronos_output/BAM/', current_group, '/'),
  chrom_size = Chromsize,
  bins = bins,
  basename = current_group,
  chr_prefix = '',
  chr_range=1:19,
  ploidy = 2,
  cores = 6
)

# Set all cells as S-phase cells
whoiswho = data.frame(Cell=SingleCell$PerCell$Cell, S_Phase = TRUE)
SingleCell$PerCell = Kronos.scRT::WhoIsWho(PerCell = SingleCell$PerCell, WhoIsWho = whoiswho)

# Use combined object's parameters from script1
Diagnostic_output$Settings$basename <- current_group
Diagnostic_output$Settings$group <- current_group

# Shift cells
SingleCell$PerCell = Kronos.scRT::AdjustPerCell(PerCell = SingleCell$PerCell,
                                                Settings = Diagnostic_output$Settings)

SingleCell$CNV = Kronos.scRT::AdjustCN(PerCell = SingleCell$PerCell, scCN = SingleCell$CNV)

# Rebin to 200k
Bins = Kronos.scRT::GenomeBinning(
  Chr_size = Chromsize,
  Chr_filter =paste0('', 1:19),
  size = 200000,
  Cores = 6
)

# Compute S-phase
SingleCell$SPhase = Kronos.scRT::Rebin(PerCell = SingleCell$PerCell,
                                       scCN = SingleCell$CNV,
                                       Bins = Bins, Sphase = T)

# Use G1 controls from script2 to compute RT 
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

# Smooth RT
SingleCell$pseudobulk = SingleCell$pseudobulk %>% 
  group_by(chr) %>%
  mutate(RT=rollapply(RT, 4, mean, fill=0)) %>%
  ungroup()

# Create per-chromosome heatmaps
for (c in 1:nrow(Chromsize)){
  
  plot <- scRTplot(SingleCell$pseudobulk,
                   S_scCN = SingleCell$SPhase,
                   Coordinates = list(chr = Chromsize$chr[c], 
                                      start = 0,
                                      end = Chromsize$size[c]),
                   rasterized_heatmap = T)
  
  ggsave(paste0('per_group/', current_group, '/',Chromsize$chr[c], '.png'), plot)
  ggsave(paste0('per_group/', current_group, '/',Chromsize$chr[c], '.svg'), plot)
}

# Save CNV per gene for each group
CNV_annotated = Kronos.scRT::TW_GenomeAnnotation(SingleCell$CNV,
                                       GenomeAnnotation = annotation) %>%
                filter(category != '_Unknown_') %>%
                select(category, copy_number, copy_number_corrected) %>% 
                arrange(desc(copy_number_corrected))

write_tsv(CNV_annotated, paste0('kronos_output/CNV/', current_group, '.tsv'))

# Compute correlation between groups
all_corr = Kronos.scRT::KCorr_plot(
  SingleCell_zygote$pseudobulk,
  SingleCell_4cell$pseudobulk,
  SingleCell_2cell$pseudobulk,
  method = 'pearson')

# Create heatmap
heatmap = ggcorrplot::ggcorrplot(cor(all_corr$data))
ggsave('kronos_output/correlation_heatmap.pdf', heatmap)
