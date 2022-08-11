library(crisprDesign)
library(devtools)
txdb <- getTxDb(organism="Homo sapiens")
txdb <- TxDb2GRangesList(txdb)
GenomeInfoDb::genome(txdb) <- "hg38"
txdb_human <- txdb
use_data(txdb_human, compress="xz", overwrite=TRUE)

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
