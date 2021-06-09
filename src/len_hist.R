install.packages('dplyr')
library(ggplot2)
library(dplyr)
# library(tidyr)   # replace_na
# library(tibble)  # column_to_rownames

###

#NAME <- 'H3K36me3.ENCFF905GSB.hg19'
#NAME <- 'H3K36me3.ENCFF905GSB.hg38'
#NAME <- 'H3K36me3.ENCFF398TIJ.hg19'
#NAME <- 'H3K36me3.ENCFF398TIJ.hg38'
#NAME <- 'H3K36me3.ENCFF398TIJ.hg19.filtered'
#NAME <- 'H3K36me3.ENCFF905GSB.hg19.filtered'
#NAME <- 'DeepZ'
NAME <- 'H3K36me3.merge.hg19'
OUT_DIR <- '/cloud/project/results'

###

bed_df <- read.delim(paste0('/cloud/project/', NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start
head(bed_df)

# hist(bed_df$len)

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.pdf'), path = OUT_DIR)

