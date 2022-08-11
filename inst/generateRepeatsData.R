# The following two files
# - repeatmaskerhg38.bed.gz 
# - repeatmaskermm10.bed.gz
# were downloaded from http://genome.ucsc.edu/cgi-bin/hgTables
# and must be in a directory ./repeats
library(readr)
library(devtools)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)

# For human:
data <- read_tsv("repeats/repeatmaskerhg38.bed.gz",
                 col_names=FALSE)
data <- as.data.frame(data)
colnames(data) <- c("chr", "start", "end",
                    "repeat", "score", "strand")
gr.repeats <- GRanges(data$chr,
                      IRanges(start=data$start, end=data$end),
                      strand=data$strand)
gr.repeats$type  <- data[["repeat"]]
gr.repeats$score <- data[["score"]]
genome(gr.repeats) <- "hg38"
gr.repeats <- gr.repeats[seqnames(gr.repeats) %in% standardChromosomes(gr.repeats)]
gr.repeats <- keepStandardChromosomes(gr.repeats)
#Add proper GenomeInfoDb info:
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
info <- seqinfo(bsgenome)
info <- keepStandardChromosomes(info)
lens  <- seqlengths(info)[as.character(seqnames(seqinfo(gr.repeats)))]
circs <- isCircular(info)[as.character(seqnames(seqinfo(gr.repeats)))]
seqlengths(seqinfo(gr.repeats)) <- lens
isCircular(seqinfo(gr.repeats)) <- circs
gr.repeats.hg38  <- gr.repeats


# For mouse:
data <- read_tsv("repeats/repeatmaskermm10.bed.gz", col_names=FALSE)
data <- as.data.frame(data)
colnames(data) <- c("chr", "start", "end","repeat", "score", "strand")
gr.repeats <- GRanges(data$chr,
                      IRanges(start=data$start, end=data$end),
                      strand=data$strand)
gr.repeats$type  <- data[["repeat"]]
gr.repeats$score <- data[["score"]]
genome(gr.repeats) <- "mm10"
gr.repeats <- gr.repeats[seqnames(gr.repeats) %in% standardChromosomes(gr.repeats)]
gr.repeats <- keepStandardChromosomes(gr.repeats)
#Add proper GenomeInfoDb info:
bsgenome <- BSgenome.Mmusculus.UCSC.mm10
info <- seqinfo(bsgenome)
info <- keepStandardChromosomes(info)
lens  <- seqlengths(info)[as.character(seqnames(seqinfo(gr.repeats)))]
circs <- isCircular(info)[as.character(seqnames(seqinfo(gr.repeats)))]
seqlengths(seqinfo(gr.repeats)) <- lens
isCircular(seqinfo(gr.repeats)) <- circs
gr.repeats.mm10  <- gr.repeats

# Saving data:
use_data(gr.repeats.hg38,
         compress="xz")
use_data(gr.repeats.mm10,
         compress="xz")



