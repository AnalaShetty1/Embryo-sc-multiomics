# Data alignment and QC  #######################################################

setwd("/home/data/BISP-666/kronos_output/")

# Alignment call
Kronos.scRT::FastqToBam(
  bowtie2_index = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107',
  File1 = list.files(path='trimmed', pattern="_R1_001.fastq"),
  File2 = list.files(path='trimmed', pattern="_R2_001.fastq"),
  outputdir = '/home/data/BISP-666/kronos_output',
  trim = F,
  cores = 10)

# Collect alignment stats:
bm_all = data.frame()
for (file in list.files(path='BAM/', pattern='.bam$')){
  print(file)
  bm <- Kronos.scRT::BamMetrics(paste0('BAM/', file), isPE = T)
  bm_all <- rbind(bm_all, bm)
}

write.table(bm_all,
            file=paste0("alignment_stats/full_stats.tsv"), 
            row.names = FALSE,
            quote=FALSE,
            sep='\t')

# To visualize stats:
Kronos.scRT::Pre_processing()

# Process data jointly to obtain joint S phase shift settings ##################
# Input data includes the G1 controls and the experiment data

# Compute bins
bins = binning(
  RefGenome = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107.fa',
  bowtie2_index = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107',
  black_list = '/home/data/BISP-666/zhan2142/src/mm10-blacklist.v2.bed',
  directory_to_bamfiles = paste0('/home/data/BISP-666/kronos_output/BAM/all/'),
  cores = 10,
  bin_size = 20000
)

# Create combined object
SingleCell = Kronos.scRT::CallCNV(
  directory = paste0('/home/data/BISP-666/kronos_output/BAM/all/'),
  chrom_size = Chromsize,
  bins = bins,
  basename = 'all',
  chr_prefix = '',
  chr_range=1:19,
  ploidy = 2,
  cores = 6
)

# Set all project cells as S-phase and all control cells as G1/G2
whoiswho = data.frame(Cell=SingleCell$PerCell$Cell, S_Phase = as.logical(str_detect(SingleCell$PerCell$Cell, '_gDNA_')))
SingleCell$PerCell = Kronos.scRT::WhoIsWho(PerCell = SingleCell$PerCell, WhoIsWho = whoiswho)

# Create diagnostic settings
Diagnostic_output = Kronos.scRT::diagnostic(SingleCell$PerCell)

# Save diagnostic plots
ggsave('plots/all_cells_plot.png',Diagnostic_output$all_cells_plot)
ggsave('plots/first_filtering_plot.png',Diagnostic_output$first_filtering_plot)
ggsave('plots/elected_G1_S_cells_plot.png',Diagnostic_output$selected_G1_S_cells_plot)
ggsave('plots/S_phase_cell_distribution_plot.png',Diagnostic_output$S_phase_cell_distribution_plot)
Diagnostic_output$Settings

# A tibble: 1 × 9
#threshold_Sphase threshold_G1G2phase Sphase_first_part Sphase_second_part Ploidy RPMPH_TH RPM_TH basename group
#<lgl>            <lgl>                           <dbl>              <dbl>  <dbl>    <int>  <dbl> <chr>    <chr>
#  1 NA             NA                            1.03              0.535   2.05       75    154   all      all  

