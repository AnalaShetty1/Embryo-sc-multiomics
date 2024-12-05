library(Kronos.scRT)

# Controls analysis and G1G2 processing ########################################

# Alignment of files
setwd("/home/data/BISP-666/references/G1_controls/")
Kronos.scRT::FastqToBam(
  bowtie2_index = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107',
  File1 = list.files(path='.', pattern=".fastq"),
  outputdir = '/home/data/BISP-666/references/G1_controls/',
  cores = 10)

# Compute alignment statistics
bm_all = data.frame()
for (file in list.files(path='BAM/', pattern='.bam$')){
  print(file)
  bm <- Kronos.scRT::BamMetrics(paste0('BAM/', file), isPE = F)
  bm_all <- rbind(bm_all, bm)
}

write.table(bm_all,file=paste0("full_stats.tsv"), row.names = FALSE, quote=FALSE, sep='\t')
Kronos.scRT::Pre_processing()

# Generation of bins
bins_controls = binning(
  RefGenome = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107.fa',
  bowtie2_index = '/home/data/BISP-666/bowtie_index/mmu_GRCm39.107',
  black_list = '/home/data/BISP-666/zhan2142/src/mm10-blacklist.v2.bed',
  directory_to_bamfiles = '/home/data/BISP-666/references/G1_controls/BAM/',
  cores = 10,
  bin_size = 20000
)

# Load chromsize file
Chromsize = readr::read_tsv(url('http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'), col_names = c('chr', 'size'))
Chromsize$chr = str_replace(Chromsize$chr, 'chr', '')
Chromsize = Chromsize[Chromsize$chr %in% seq(1,19),]

# Generate controls object
SingleCell_controls = Kronos.scRT::CallCNV(
  directory = '/home/data/BISP-666/references/G1_controls/BAM/',
  chrom_size = Chromsize,
  bins = bins_controls,
  basename = 'control',
  chr_prefix = '',
  chr_range=1:19,
  ploidy = 2,
  cores = 6
)

# Set all cells to G1/G2
whoiswho = data.frame(Cell=SingleCell_controls$PerCell$Cell, S_Phase = FALSE)
SingleCell_controls$PerCel$PerCell = Kronos.scRT::WhoIsWho(PerCell = SingleCell_controls$PerCel, WhoIsWho = whoiswho)

# Load previous settings and use them for shifting of the controls
Diagnostic_output$Settings$basename <- 'control'
Diagnostic_output$Settings$group <- 'control'

# Shift cells
SingleCell_controls$PerCell = Kronos.scRT::AdjustPerCell(PerCell = SingleCell_controls$PerCell,
                                                         Settings = Diagnostic_output$Settings)

SingleCell_controls$CNV = Kronos.scRT::AdjustCN(PerCell = SingleCell_controls$PerCell, 
                                                scCN = SingleCell_controls$CNV)

# Create genome bins
Bins = Kronos.scRT::GenomeBinning(
  Chr_size = Chromsize,
  Chr_filter =paste0('', 1:19),
  size = 20000,
  Cores = 6
)

# Re-bin using the genome bins
SingleCell_controls$G1G2 = Kronos.scRT::Rebin(PerCell = SingleCell_controls$PerCell,
                                              scCN = SingleCell_controls$CNV,
                                              Bins = Bins, 
                                              Sphase = F)

# Compute background
SingleCell_controls$MedianG1G2 = Kronos.scRT::BackGround(G1_scCN = SingleCell_controls$G1G2)

# Compute replication
SingleCell_controls$G1G2 = Kronos.scRT::Replication_state(
  Samples = SingleCell_controls$G1G2,
  background = SingleCell_controls$MedianG1G2,
  cores = 6
)
