# For each group, generate dataframe merged with the references

#zygote
df <- merge(SingleCell_zygote$pseudobulk %>% mutate(bin=paste0(chr, ':', start, '-', end)),
            ref_xu_zygote %>% mutate(bin=paste0(chr, ':', start, '-', end)), by.x='bin', by.y='bin') %>% 
  merge(ref_natakani_zygote  %>% mutate(bin=paste0(chr, ':', start, '-', end))) %>% arrange(chr.x, start.x) %>%
  select(chr = chr.x,RT_zygote = RT.x, RT_xu = RT.y, RT_natakani = RT)

# 2cell
df <- merge(SingleCell_2cell$pseudobulk %>% mutate(bin=paste0(chr, ':', start, '-', end)),
            ref_xu_2cell %>% mutate(bin=paste0(chr, ':', start, '-', end)), by.x='bin', by.y='bin') %>% 
  merge(ref_natakani_2cell %>% mutate(bin=paste0(chr, ':', start, '-', end))) %>% arrange(chr.x, start.x) %>%
  select(chr = chr.x,RT_zygote = RT.x, RT_xu = RT.y, RT_natakani = RT)

# 4cell
df <- merge(SingleCell_4cell$pseudobulk %>% mutate(bin=paste0(chr, ':', start, '-', end)),
            ref_xu_4cell %>% mutate(bin=paste0(chr, ':', start, '-', end)), by.x='bin', by.y='bin') %>% 
  merge(ref_natakani_4cell %>% mutate(bin=paste0(chr, ':', start, '-', end))) %>% arrange(chr.x, start.x) %>%
  select(chr = chr.x,RT_zygote = RT.x, RT_xu = RT.y, RT_natakani = RT)

# Compute early/late bins between data and Natakani reference
print(nrow(df[df$RT_zygote > 0.5 & df$RT_natakani > 0.5,]))
print(nrow(df[df$RT_zygote <= 0.5 & df$RT_natakani <= 0.5,]))
print(nrow(df[df$RT_zygote <= 0.5 & df$RT_natakani > 0.5,]))
print(nrow(df[df$RT_zygote > 0.5 & df$RT_natakani <= 0.5,]))

# Compute early/late bins between data and Xu reference
print(nrow(df[df$RT_zygote > 0.5 & df$RT_xu > 0.5,]))
print(nrow(df[df$RT_zygote <= 0.5 & df$RT_xu <= 0.5,]))
print(nrow(df[df$RT_zygote <= 0.5 & df$RT_xu > 0.5,]))
print(nrow(df[df$RT_zygote > 0.5 & df$RT_xu <= 0.5,]))

# Compute early/late bins between both references
print(nrow(df[df$RT_natakani > 0.5 & df$RT_xu > 0.5,]))
print(nrow(df[df$RT_natakani <= 0.5 & df$RT_xu <= 0.5,]))
print(nrow(df[df$RT_natakani <= 0.5 & df$RT_xu > 0.5,]))
print(nrow(df[df$RT_natakani > 0.5 & df$RT_xu <= 0.5,]))

# Final piecharts were created in Excel