library(crisprDesign)
library(crisprDesignData)
library(devtools)
library(BSgenome.Hsapiens.UCSC.hg38)
data(txdb_human, package="crisprDesignData")
txObject <- txdb_human
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
ids <- unique(txObject$exons$tx_id)
mrnasHuman <- getMrnaSequences(ids,
                               bsgenome=bsgenome,
                               txObject=txObject)
use_data(mrnasHuman, compress="xz")
#library(Biostrings)
#writeXStringSet(mrnasHuman,
#                file="ensembl_human_104.fasta",
#                format="fasta")


library(BSgenome.Mmusculus.UCSC.mm10)
data(txdb_mouse, package="crisprDesignData")
txObject <- txdb_mouse
bsgenome <- BSgenome.Mmusculus.UCSC.mm10
ids <- unique(txObject$exons$tx_id)
mrnasMouse <- getMrnaSequences(ids,
                               bsgenome=bsgenome,
                               txObject=txObject)
use_data(mrnasMouse, compress="xz")
#library(Biostrings)
#writeXStringSet(mrnasMouse,
#                file="ensembl_mouse_102.fasta",
#                format="fasta")


