if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

install.packages("uuid")
install.packages("shadowtext")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)

getwd()
setwd("C:/Users/msc2.MSC2-15-226/Desktop/4thsem_practicals/Advanced_Genomics")
peaks <- read.table("peaks.bed", header = FALSE)
peaks

# Convert to GRanges
library(GenomicRanges)

gr <- makeGRangesFromDataFrame(
  peaks,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3"
)

# Load annotation databases
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks
peakAnno <- annotatePeak(
  gr,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

# Convert results
peakAnno_DF <- as.data.frame(peakAnno)
peakAnno_G  <- as.GRanges(peakAnno)

# Plot results
plotAnnoBar(peakAnno)

# Basic pie chart with proper labels inside
plotAnnoPie(peakAnno)

anno <- as.data.frame(peakAnno@annoStat)

library(ggplot2)

ggplot(anno, aes(x = "", y = Frequency, fill = Feature)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  theme(legend.position = "right") +
  labs(fill = "Genomic Annotation")

write.table(peakAnno_DF, file = "filename.txt", sep = "\t", quote = FALSE, row.names = FALSE)

