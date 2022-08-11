library(dplyr)
library(crisprDesign)
library(devtools)
data <- read.csv("appris/appris_data.principal_human.txt.gz", sep="\t")
data <- data[, c("Gene.ID", "Transcript.ID", "APPRIS.Annotation", "MANE")]
data$MANE[data$MANE=="MANE_Plus_Clinical"] <- ""
data <- dplyr::rename(data, gene_id="Gene.ID")
data <- dplyr::rename(data, tx_id="Transcript.ID")
data <- dplyr::rename(data, mane_select="MANE")
data <- dplyr::rename(data, appris="APPRIS.Annotation")
data$mane_select <- data$mane_select!=""

load("../data/txdb_human.rda")
key <- crisprDesign:::.getTx2GeneTable(txdb_human)
key$id <- paste0(key$gene_id, "_", key$tx_id)
data$id <- paste0(data$gene_id, "_", data$tx_id)
data <- data[data$id %in% key$id,]
data$id <- NULL
data$appris_label  <- stringr::str_extract(data$appris, "PRINCIPAL|ALTERNATIVE")
data$appris_number <- as.integer(stringr::str_extract(data$appris, "[0-9]+"))
apprisHuman <- data
use_data(apprisHuman, compress="xz")


data <- read.table("appris/appris_data.principal_mouse.txt.gz", head=FALSE)
data <- data[,c(2,3,5)]
colnames(data) <- c("gene_id", "tx_id", "appris")

library(devtools)
load("../data/txdb_mouse.rda")
key <- crisprDesign:::.getTx2GeneTable(txdb_mouse)
key$id <- paste0(key$gene_id, "_", key$tx_id)
data$id <- paste0(data$gene_id, "_", data$tx_id)
data <- data[data$id %in% key$id,]
data$id <- NULL
data$appris_label  <- stringr::str_extract(data$appris, "PRINCIPAL|ALTERNATIVE")
data$appris_number <- as.integer(stringr::str_extract(data$appris, "[0-9]+"))
apprisMouse <- data
use_data(apprisMouse, compress="xz")

