#' Repeat elements coordinates found in the human genome (hg38/GRCh38)
#'
#' Repeat elements coordinates for hg38 obtained from the RepeatMasker track
#' found at the UCSC genome browser. 
#' 
#' @format GRanges
"gr.repeats.hg38"


#' Repeat elements coordinates found in the mouse genome (mm10/GRCm38)
#'
#' Repeat elements coordinates for mm10 obtained from the RepeatMasker track
#' found at the UCSC genome browser. 
#' 
#' @format GRanges
"gr.repeats.mm10"


#' GRangesList representing Ensembl gene model coordinates for human
#'
#' GRangesList representing Ensembl gene model coordinates for human in 
#' hg38/GRCh38 coordinates (Ensembl release 104)
#' 
#' @format GRangesList
"txdb_human"


#' GRangesList representing Ensembl gene model coordinates for mouse
#'
#' GRangesList representing Ensembl gene model coordinates for mouse in 
#' mm10/GRCm38 coordinates (Ensembl release 102)
#' 
#' @format GRangesList
"txdb_mouse"


#' GRangesList representing Ensembl gene model coordinates for human
#'
#' GRanges representing Ensembl-based TSS coordinates for human in 
#' hg38/GRCh38 coordinates (Ensembl release 104)
#' 
#' @format GRanges
"tss_human"


#' GRanges representing Ensembl-based TSS coordinates for mouse
#'
#' GRanges representing Ensembl-based TSS coordinates for mouse in mm10/GRCm38
#' coordinates (Ensembl release 102)
#' 
#' @format GRanges
"tss_mouse"



#' Genomic coordinates of FANTOM5 CAGE Peaks in hg38 coordinates
#'
#' Genomic coordinates of FANTOM5 CAGE Peaks in hg38 coordinates.
#' 
#' @format GRanges
"humanCagePeaks"



#' Genomic coordinates of FANTOM5 CAGE Peaks in mm10 coordinates
#'
#' Genomic coordinates of FANTOM5 CAGE Peaks in mm10 coordinates.
#' 
#' @format GRanges
"mouseCagePeaks"




#' DNAStringSet of mRNA sequences of human transcripts (Ensembl)
#'
#' DNAStringSet of mRNA sequences of human transcripts
#' (Ensembl release 104) extracted from the txdb_human object.
#' 
#' @format DNAStringSet
"mrnasHuman"


#' DNAStringSet of mRNA sequences of mouse transcripts (Ensembl)
#'
#' DNAStringSet of mRNA sequences of mouse transcripts
#' (Ensembl release 102) extracted from the txdb_human object.
#' 
#' @format DNAStringSet
"mrnasMouse"



#' APPRIS annotation for human genes (Ensembl release 104)
#'
#' data.frame containing APPRIS annotation for human genes
#' using the APPRIS release 47 and the Ensembl release 104
#' gene annotation.
#' 
#' @format data.frame
"apprisHuman"


#' APPRIS annotation for mouse genes (Ensembl release 102)
#'
#' data.frame containing APPRIS annotation for mouse genes
#' using the APPRIS release 47 and the Ensembl release 102
#' gene annotation.
#' 
#' @format data.frame
"apprisMouse"


#' Canonical Ensembl transcripts for human genes (Ensembl release 104)
#'
#' data.frame containing canonical Ensembl transcripts for human genes.
#' First column contains transcript ID, second column contains gene ID.
#' 
#' @format data.frame
"canonicalHuman"


#' Canonical Ensembl transcripts for mouse genes (Ensembl release 102)
#'
#' data.frame containing canonical Ensembl transcripts for mouse genes.
#' First column contains transcript ID, second column contains gene ID.
#' 
#' @format data.frame
"canonicalMouse"
