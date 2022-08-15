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






getCanonicalTranscripts <- function(gene_ids,
                                    organism=c("hsapiens", "mmusculus")
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
                values=ids,
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


use_data(canonicalHuman, compress="xz", overwrite=TRUE)
use_data(canonicalMouse, compress="xz", overwrite=TRUE)




