library(crisprDesign)
library(devtools)
load("../data/txdb_human.rda")
load("../data/txdb_mouse.rda")
pfamTableHuman <- preparePfamTable(txdb_human, mart_dataset="hsapiens_gene_ensembl")
pfamTableMouse <- preparePfamTable(txdb_mouse, mart_dataset="mmusculus_gene_ensembl")
pfamTableHuman <- pfamTableHuman[pfamTableHuman$pfam!="",]
pfamTableMouse <- pfamTableMouse[pfamTableMouse$pfam!="",]

use_data(pfamTableHuman, compress="xz")
use_data(pfamTableMouse, compress="xz")
