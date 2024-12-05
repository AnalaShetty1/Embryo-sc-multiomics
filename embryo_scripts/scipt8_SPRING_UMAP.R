# Generate input data for the SPRING algorithm
###########################################################

# Merge S-phase cells across stages and select CN_bg values per cell per bin
mat_data_sphase <- rbind(SingleCell_zygote$SPhase, SingleCell_2cell$SPhase, SingleCell_4cell$SPhase) %>%
  mutate(bin = paste0(chr, ":", start, "-", end)) %>%
  dplyr::select(Cell, CN_bg, bin) %>%
  mutate(CN_bg = as.numeric(CN_bg)) %>%
  pivot_wider(names_from = Cell, values_from = CN_bg)

# Fill in NAs
mat_data_sphase[is.na(mat_data_sphase)] <- 0

# Extract controls values per cell per bin
mat_data_gphase <- SingleCell_controls$G1G2 %>%
  mutate(bin = paste0(chr, ":", start, "-", end)) %>%
  dplyr::select(Cell, CN_bg, bin) %>%
  mutate(CN_bg = as.numeric(CN_bg)) %>%
  pivot_wider(names_from = Cell, values_from = CN_bg)

# Fill in NAs
mat_data_gphase[is.na(mat_data_gphase)] <- 0

# Merge S and G1/2 phase cells
mat_data <- merge(mat_data_sphase, mat_data_gphase, by.x="bin", by.y="bin", sort=FALSE)

# Save files for SRPING algorithm
# Plot and cell embeddings were generated using Spring viewer:
# https://kleintools.hms.harvard.edu/tools/spring.html
####################################################

# Save final matrix [referenced as Expression matrix]
write.table(mat_data[,2:ncol(mat_data)], file='umap/SPRING_input.tsv', row.names = FALSE, col.names = FALSE, sep = '\t')

# Save bins [referenced as gene list]
write.table(mat_data$bin, file='umap/SPRING_input_bins.tsv', row.names = FALSE, col.names = FALSE, sep = '\t')

# Cell grouping file
c <- str_split(colnames(mat_data[,2:ncol(mat_data)]), '_')
phase = c(ifelse(colnames(mat_data[,2:ncol(mat_data)]) %in% colnames(mat_data_sphase), unlist(c)[8*(1:nrow(umap$layout))-7], 'G1/G2'))
write.table(paste0('Phase,',paste(phase, collapse = ',')), file='umap/SPRING_input_cells.tsv', row.names = FALSE, col.names = FALSE, quote = FALSE)

# Compute UMAP embeddings for the combined matrix
umap <- umap::umap(t(mat_data[,2:ncol(mat_data)]))
umap$layout <- as.data.frame(umap$layout)
c <- str_split(rownames(umap$layout), '_')
umap$layout$stage <- ifelse(str_detect(rownames(umap$layout), '_trimmed.bam'), 'G1 control', unlist(c)[8*(1:nrow(umap$layout))-7])
umap$layout$RT <- colMeans(mat_data[,2:ncol(mat_data)])

# Plot RT UMAP
plot <- ggplot(umap$layout, aes(x=V1, y=V2, color=RT)) +
  geom_density_2d(colour='lightgrey', bins=6) +
  geom_point(size=2.2, alpha= 0.8) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  labs(title='RT UMAP') 
scale_color_brewer(palette="Paired") 

ggsave('umap/RT.pdf', plot)
