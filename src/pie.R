

#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("PChIp")
#BiocManager::install("annotatePeak")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#BiocManager::install("GenomicFeatures")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
###

NAME <- 'H3K36me3.ENCFF398TIJ.hg19.filtered'
#NAME <- 'H3K36me3.ENCFF905GSB.hg19.filtered'
#NAME <- 'DeepZ'

DATA_DIR <- '/cloud/project/results/'
BED_FN <- paste0(DATA_DIR, NAME, '.bed')
OUT_DIR <- '/cloud/project/data/'

###

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(readPeakFile(BED_FN), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

png(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.png'))
plotAnnoPie(peakAnno)
dev.off()

# peak <- readPeakFile(BED_FN)
# pdf(paste0(OUT_DIR, 'chip_seeker.', NAME, '.covplot.pdf'))
# covplot(peak, weightCol="V5")
# dev.off()
# 

