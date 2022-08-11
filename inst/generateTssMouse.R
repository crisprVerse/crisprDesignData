library(crisprDesign)
library(devtools)
score_threshold <- 4
dist_treshold <- 500

# First, getting Ensembl-based TSSs:
load("../data/txdb_mouse.rda")
tss <- getTssObjectFromTxObject(txdb_mouse)
# Getting FANTOM5 peaks:
load("../data/mouseCagePeaks.rda")
cage <- mouseCagePeaks
cage$tx_id <- cage$txid
cage$txid <- NULL
cage <- cage[cage$tx_id %in% tss$tx_id]
wh <- match(cage$tx_id, tss$tx_id)
cage$gene_id <- tss$gene_id[wh]
cage$score <- log10(cage$score+1)

cage <- as.data.frame(cage)
cage <- dplyr::rename(cage, chr=seqnames)
cage$coordid <- paste0(cage$chr, ":",cage$start)
dfs <- split(cage, f=cage$gene_id)
dfs <- lapply(dfs, function(df){
    df[!duplicated(df$coordid),,drop=FALSE]
})
dfs <- lapply(dfs, function(df){
    df <- df[order(-df$score),, drop=FALSE]
    df
})

dfs <- lapply(dfs, function(df){
    n <- sum(df$score>=score_threshold)
    if (n==0){
        n <- 1
    }
    df <- df[1:n,,drop=FALSE]
    return(df)
})
dfs <- lapply(dfs, function(df){
    n <- nrow(df)
    if (n>1){
        good <- rep(TRUE, n)
        for (k in 2:n){
            if (df[k,"tx_id"] %in% df[1:(k-1), "tx_id"]){
                good[k] <- FALSE
            }
        }
        df <- df[good,,drop=FALSE]
    }
    df
})


# Removing TSSs that are too close:
dfs <- lapply(dfs, function(df){
    n <- nrow(df)
    if (n>1){
        good <- rep(TRUE, n)
        for (k in 2:n){
            pos <- df[k, "start"]
            diffs <- abs(pos- df[1:(k-1), "start"])
            if (any(diffs<=dist_treshold)){
                good[k] <- FALSE
            }
        }
        df <- df[good,,drop=FALSE]
    }
    return(df)
})

# Adding ID
dfs <- lapply(dfs, function(df){
    df$promoter <-paste0("P", 1:nrow(df))
    df$ID <- paste0(df$gene_id, "_", df$promoter)
    df
})


final <- dplyr::bind_rows(dfs)
final$coordid <- NULL

# Now getting remaining TSSs from Ensembl
# Need principal tx information, etc:
tss <- tss[!tss$gene_id %in% final$gene_id]
appris <- read.table("appris/appris_data.principal_mouse.txt.gz",
                     sep="\t", head=TRUE)
appris <- appris[,c(2,3, 5)]
colnames(appris) <- c("gene_id", "tx_id", "label")
appris <- appris[appris$tx_id %in% tss$tx_id,]
tss <- tss[tss$tx_id %in% appris$tx_id]
wh <- match(tss$tx_id, appris$tx_id)
tss$appris <- appris$label[wh]
df <- data.frame(tss)
df$appris_label <- stringr::str_extract(df$appris, "PRINCIPAL|ALTERNATIVE")
df$appris_number <- stringr::str_extract(df$appris, "[0-9]+")
df$appris_priority <- ifelse(df$appris_label=="PRINCIPAL",1,2)
dfs <- split(df, f=df$gene_id)
dfs <- lapply(dfs, function(df){
    df <- df[order(df$appris_priority, df$appris_number),,drop=FALSE]
    df
})
# Removing TSSs that are too close:
dfs <- lapply(dfs, function(df){
    n <- nrow(df)
    if (n>1){
        good <- rep(TRUE, n)
        for (k in 2:n){
            pos <- df[k, "start"]
            diffs <- abs(pos- df[1:(k-1), "start"])
            if (any(diffs<=dist_treshold)){
                good[k] <- FALSE
            }
        }
        df <- df[good,,drop=FALSE]
    }
    return(df)
})
dfs <- lapply(dfs, function(df){
    df$promoter <- paste0("P", 1:nrow(df))
    df$ID <- paste0(df$gene_id, "_", df$promoter)
    df
})
final_ensembl <- dplyr::bind_rows(dfs)
final_ensembl <- dplyr::rename(final_ensembl, chr=seqnames)
final_ensembl$score <- NA
final_ensembl$peak_start <- NA
final_ensembl$peak_end <- NA
final_ensembl$gene_symbol <- NULL
final_ensembl <- final_ensembl[, colnames(final)]
final_ensembl$source <- "appris"
final$source <- "fantom5"

finalTss <- rbind(final, final_ensembl)

# Are we missing any genes at all?
ids <- unique(txdb_mouse$cds$gene_id)
missing <- setdiff(ids, finalTss$gene_id)
tss <- getTssObjectFromTxObject(txdb_mouse)
tss <- tss[tss$gene_id %in% missing]
tss <- as.data.frame(tss)
dfs <- split(tss, f=tss$gene_id)
dfs <- lapply(dfs, function(df){
    n <- nrow(df)
    if (n>1){
        good <- rep(TRUE, n)
        for (k in 2:n){
            pos <- df[k, "start"]
            diffs <- abs(pos- df[1:(k-1), "start"])
            if (any(diffs<=dist_treshold)){
                good[k] <- FALSE
            }
        }
        df <- df[good,,drop=FALSE]
    }
    return(df)
})
dfs <- lapply(dfs, function(df){
    df$promoter <- paste0("P", 1:nrow(df))
    df$ID <- paste0(df$gene_id, "_", df$promoter)
    df
})
final_missing <- dplyr::bind_rows(dfs)
final_missing$source <- "ensembl"
final_missing <- dplyr::rename(final_missing, chr=seqnames)
final_missing$score <- NA
final_missing$peak_start <- NA
final_missing$peak_end <- NA
final_missing <- final_missing[,colnames(finalTss)]

finalTss <- rbind(finalTss, final_missing)


library(GenomicRanges)
gr <- GRanges(finalTss$chr,
              IRanges(start=finalTss$start,
                      end=finalTss$start),
              strand=finalTss$strand)
gr$score <- finalTss$score                                
gr$peak_start <- finalTss$start
gr$peak_end <- finalTss$end
gr$tx_id <- finalTss$tx_id
gr$gene_id <- finalTss$gene_id
gr$source <- finalTss$source
gr$promoter <- finalTss$promoter
gr$ID <- finalTss$ID
# Let's add gene symbol
wh <- match(gr$gene_id, txdb_mouse$cds$gene_id)
gr$gene_symbol <- txdb_mouse$cds$gene_symbol[wh]
names(gr) <- gr$ID
tss_mouse <- gr
use_data(tss_mouse, compress="xz", overwrite=TRUE)

