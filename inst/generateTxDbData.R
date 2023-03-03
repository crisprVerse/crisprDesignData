library(crisprDesign)
library(devtools)
library(biomaRt)


### Human
txdb <- getTxDb(organism="Homo sapiens",
                release=104)
txdb <- TxDb2GRangesList(txdb)
GenomeInfoDb::genome(txdb) <- "hg38"
txdb_human <- txdb
use_data(txdb_human,
         compress="xz",
         overwrite=TRUE)


### Mouse
txdb <- getTxDb(organism="Mus musculus",
                release=102) # latest release for GRCm38
txdb <- TxDb2GRangesList(txdb)
GenomeInfoDb::genome(txdb) <- "mm10"
txdb_mouse <- txdb
use_data(txdb_mouse, compress="xz", overwrite=TRUE)





txdb <- getTxDb(organism="Homo sapiens",
                release=106)
txdb <- TxDb2GRangesList(txdb)
GenomeInfoDb::genome(txdb) <- "hg38"
txdb_human_106 <- txdb
use_data(txdb_human_106, compress="xz", overwrite=TRUE)


# Cyno 6.0:
txdb <- getTxDb(organism="Macaca fascicularis",
                tx_attrib=NULL)
txdb <- TxDb2GRangesList(txdb,
                         standardChromOnly=FALSE,
                         seqlevelsStyle="NCBI")
GenomeInfoDb::genome(txdb) <- "Macaca_fascicularis_6.0"
txdb_cyno <- txdb
#use_data(txdb_cyno, compress="xz", overwrite=TRUE)

library(BSgenome.Mfascicularis.NCBI.6.0)
bsgenome <- BSgenome.Mfascicularis.NCBI.6.0
levels1 <- GenomeInfoDb::seqlevels(txdb_cyno)
levels2 <- GenomeInfoDb::seqlevels(bsgenome)
common <- intersect(levels1, levels2)
tx <- txdb_cyno
GenomeInfoDb::genome(tx) <- "Macaca_fascicularis_6.0"
for (i in 1:length(tx)){
    tx[[i]] <- tx[[i]][seqnames(tx[[i]]) %in% common]
}
seqlevels(tx) <- seqlevelsInUse(tx)
txdb_cyno <- tx
use_data(txdb_cyno, compress="xz", overwrite=TRUE)


# Cyno 5.0:
txdb <- getTxDb(organism="Macaca fascicularis",
                release=102,
                tx_attrib=NULL)
txdb <- TxDb2GRangesList(txdb,
                         standardChromOnly=FALSE,
                         seqlevelsStyle="NCBI")
GenomeInfoDb::genome(txdb) <- "Macaca_fascicularis_5.0"
txdb_cyno_ensembl102 <- txdb
#use_data(txdb_cyno_ensembl102, compress="xz", overwrite=TRUE)
library(BSgenome.Mfascicularis.NCBI.5.0)
bsgenome <- BSgenome.Mfascicularis.NCBI.5.0
levels1 <- GenomeInfoDb::seqlevels(txdb_cyno_ensembl102)
levels1_standard <- GenomeInfoDb::standardChromosomes(txdb_cyno_ensembl102)
levels1_standard <- levels1_standard[levels1_standard != "MT"]
levels1[levels1 %in% levels1_standard] <- paste0("MFA", levels1_standard)
GenomeInfoDb::seqlevels(txdb_cyno_ensembl102) <- levels1
levels2 <- GenomeInfoDb::seqlevels(bsgenome)
common <- intersect(levels1, levels2)
tx <- txdb_cyno_ensembl102
GenomeInfoDb::genome(tx) <- "Macaca_fascicularis_5.0"
for (i in 1:length(tx)){
    tx[[i]] <- tx[[i]][seqnames(tx[[i]]) %in% common]
}
seqlevels(tx) <- seqlevelsInUse(tx)
txdb_cyno_ensembl102 <- tx
use_data(txdb_cyno_ensembl102, compress="xz", overwrite=TRUE)





getCanonicalTranscripts <- function(gene_ids,
                                    organism=c("hsapiens",
                                               "mmusculus",
                                               "mfascicularis")
){
    organism <- match.arg(organism)
    mart <- biomaRt::useMart("ensembl",
                             dataset=paste0(organism,
                                            "_gene_ensembl"))
    filters <- c("ensembl_gene_id")
    x1 <- "ensembl_transcript_id"
    x2 <- "transcript_is_canonical"
    x3 <- "ensembl_gene_id"
    attributes <- c(x1,x2,x3)
    df <- getBM(attributes=attributes,
                filters=filters,
                values=gene_ids,
                mart=mart)
    df <- df[which(df$transcript_is_canonical==1),,drop=FALSE]
    df$transcript_is_canonical <- NULL
    rownames(df) <- NULL
    colnames(df) <- c("tx_id", "gene_id")
    return(df)
}


library(crisprDesignData)
ids <- unique(txdb_human$cds$gene_id)
canonicalHuman <- getCanonicalTranscripts(ids)
setdiff(ids, canonicalHuman$gene_id)
ids <- unique(txdb_mouse$cds$gene_id)
canonicalMouse <- getCanonicalTranscripts(ids,
                                          organism="mmusculus")

ids <- unique(txdb_cyno$cds$gene_id)
canonicalCyno <- getCanonicalTranscripts(ids,
                                         organism="mfascicularis")


cleanCanonical <- function(canonicalTxs,
                           appris=NULL,
                           txObject
){
    cols <- c("gene_id", "tx_id", "gene_symbol")
    tx2Gene <- mcols(txObject[["cds"]])[, cols, drop = FALSE]
    tx2Gene <- as.data.frame(tx2Gene)
    tx2Gene <- tx2Gene[!duplicated(tx2Gene), ]
    rownames(tx2Gene) <- NULL

    # Only keeping transcripts that exist in current annotation:
    good <- canonicalTxs$tx_id %in% tx2Gene$tx_id
    canonicalTxs <- canonicalTxs[good,]

    if (!is.null(appris)){
        # Let's add MANE select transcripts:
        missing <- setdiff(tx2Gene$gene_id, canonicalTxs$gene_id)
        appris <- appris[appris$gene_id %in% missing,]
        mane <- appris[appris$mane_select,]
        stopifnot(!any(duplicated(mane$gene_id)))
        mane <- mane[, c("gene_id", "tx_id")]
        canonicalTxs <- rbind(canonicalTxs, mane)

        # Let's add Appris principal:
        missing <- setdiff(tx2Gene$gene_id, canonicalTxs$gene_id)
        appris <- appris[appris$gene_id %in% missing,]
        principal <- appris[appris$appris_label=="PRINCIPAL",]
        principal <- principal[order(principal$gene_id, principal$appris_number),]
        missing <- intersect(missing, principal$gene_id)
        wh <- match(missing, principal$gene_id)
        principal <- principal[wh,]
        principal <- principal[, c("gene_id", "tx_id")]
        canonicalTxs <- rbind(canonicalTxs, principal)

        # Let's add Appris alternative:
        missing <- setdiff(tx2Gene$gene_id, canonicalTxs$gene_id)
        appris <- appris[appris$gene_id %in% missing,]
        missing <- intersect(missing, appris$gene_id)
        wh <- match(missing, appris$gene_id)
        appris <- appris[wh,]
        appris <- appris[, c("gene_id", "tx_id")]
        canonicalTxs <- rbind(canonicalTxs, appris)
    }

    # Are there any still missing?
    # We will select the longest isoform for missing genes:
    missing <- setdiff(tx2Gene$gene_id, canonicalTxs$gene_id)
    tx2Gene <- tx2Gene[tx2Gene$gene_id %in% missing,]
    cds <- txObject[["cds"]]
    cds <- cds[cds$gene_id %in% missing]
    dfs <- split(cds, f=cds$tx_id)
    ns <- vapply(dfs, function(x) sum(BiocGenerics::width(x)), FUN.VALUE=1)
    wh <- match(tx2Gene$tx_id, names(ns))
    tx2Gene$len <- ns[wh]
    dfs <- split(tx2Gene, f=tx2Gene$gene_id)
    dfs <- lapply(dfs, function(df){
        df <- df[order(-df$len),,drop=FALSE]
        df[1,,drop=FALSE]
    })
    df <- do.call(rbind, dfs)
    df <- df[, c("gene_id", "tx_id")]
    canonicalTxs <- rbind(canonicalTxs, df)
    return(canonicalTxs)
}

# Going to clean the canonical transcripts for current release:
#load("../data/canonicalHuman.rda")
load("../data/txdb_human.rda")
load("../data/apprisHuman.rda")
canonicalHuman <- cleanCanonical(canonicalHuman,
                                 appris=apprisHuman,
                                 txObject=txdb_human)
#load("../data/canonicalMouse.rda")
load("../data/txdb_mouse.rda")
load("../data/apprisMouse.rda")
canonicalMouse <- cleanCanonical(canonicalMouse,
                                 appris=apprisMouse,
                                 txObject=txdb_mouse)



#load("../data/canonicalMouse.rda")
load("../data/txdb_cyno.rda")
canonicalCyno <- cleanCanonical(canonicalCyno,
                                appris=NULL,
                                txObject=txdb_cyno)


use_data(canonicalHuman, compress="xz", overwrite=TRUE)
use_data(canonicalMouse, compress="xz", overwrite=TRUE)
use_data(canonicalCyno, compress="xz", overwrite=TRUE)



