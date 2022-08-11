# The file hg38_liftover+new_CAGE_peaks_phase1and2.bed.gz must be 
# in a directory called /fantom5
library(devtools)
library(GenomicRanges)
library(stringr)
library(dplyr)
data <- read.table("fantom5/human/bed/hg38_liftover+new_CAGE_peaks_phase1and2.bed.gz")
cols <- c(1,2,3,4,5,6,7)
data <- data[,cols]
colnames(data) <- c("chr", "start", "end","cageid", "score","strand", "tss_pos")
rownames(data) <- paste0("peak_", seq_len(nrow(data)))



createTranscriptTable <- function(ann){
    cols <- c(1,2)
    ann <- ann[,cols]
    colnames(ann) <- c("cageid", "txid")
    ann <- ann[ann$txid!="",]
    txs <- ann$txid
    txs <- lapply(txs, function(x) {
      strsplit(x, split=" ")[[1]]
    })
    txs <- lapply(txs, function(tx){
        tx <- tx[grepl("ENS", tx)]
        tx
    })
    ns <- vapply(txs, length, FUN.VALUE=1)
    good <- which(ns>0)
    ann <- ann[good,]
    txs <- txs[good]
    txs <- lapply(txs, function(x){
      data.frame(txid=x)
    })
    for (k in 1:length(txs)){
      txs[[k]]$cageid <- ann[k, "cageid"]
    }
    out <- dplyr::bind_rows(txs)
    out <- out[,2:1]
    out$txid <- str_extract(out$txid,
                            "ENST[0-9]+|ENSMUST[0-9]+")
    return(out)
}
# Let's get the annotation
ann <- read.csv("fantom5/human/annot/hg38_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz",
                sep="\t")
txTable <- createTranscriptTable(ann)

data <- data[data$cageid %in% txTable$cageid,]
data <- full_join(data, txTable, by="cageid")



gr <- GRanges(data$chr,
              IRanges(start=data$tss_pos,
                      end=data$tss_pos),
              strand=data$strand)
gr$score <- data$score                                
gr$peak_start <- data$start
gr$peak_end <- data$end
gr$txid <- data$txid
gr$ID <- paste0("cagePeak_", 1:length(gr))
names(gr) <- gr$ID
humanCagePeaks <- gr
use_data(humanCagePeaks,
         compress="xz",
         overwrite=TRUE)



# For mouse
data <- read.table("fantom5/mouse/bed/mm10_liftover+new_CAGE_peaks_phase1and2.bed.gz")
cols <- c(1,2,3,4,5,6,7)
data <- data[,cols]
colnames(data) <- c("chr", "start", "end","cageid", "score","strand", "tss_pos")


# Let's get the annotation
ann <- read.csv("fantom5/mouse/annot/mm10_liftover+new_CAGE_peaks_phase1and2_annot.txt.gz",
                sep="\t")
txTable <- createTranscriptTable(ann)

data <- data[data$cageid %in% txTable$cageid,]
data <- full_join(data, txTable, by="cageid")



gr <- GRanges(data$chr,
              IRanges(start=data$tss_pos,
                      end=data$tss_pos),
              strand=data$strand)
gr$score <- data$score                                
gr$peak_start <- data$start
gr$peak_end <- data$end
gr$txid <- data$txid
gr$ID <- paste0("cagePeak_", 1:length(gr))
names(gr) <- gr$ID
mouseCagePeaks <- gr
use_data(mouseCagePeaks,
         compress="xz",
         overwrite=TRUE)








