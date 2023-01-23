#' Repeat elements coordinates found in the human genome (hg38/GRCh38)
#'
#' Repeat elements coordinates for hg38 obtained from the RepeatMasker track
#' found at the UCSC genome browser. 
#' 
#' @format GRanges
#' @usage data(gr.repeats.hg38, package="crisprDesignData")
"gr.repeats.hg38"


#' Repeat elements coordinates found in the mouse genome (mm10/GRCm38)
#'
#' Repeat elements coordinates for mm10 obtained from the RepeatMasker track
#' found at the UCSC genome browser. 
#' 
#' @format GRanges
#' @usage data(gr.repeats.mm10, package="crisprDesignData")
"gr.repeats.mm10"


#' GRangesList representing Ensembl gene model coordinates for human
#'
#' GRangesList representing Ensembl gene model coordinates for human in 
#' hg38/GRCh38 coordinates (Ensembl release 104)
#' 
#' @format GRangesList
#' @usage data(txdb_human, package="crisprDesignData")
"txdb_human"


#' GRangesList representing Ensembl gene model coordinates for mouse
#'
#' GRangesList representing Ensembl gene model coordinates for mouse in 
#' mm10/GRCm38 coordinates (Ensembl release 102)
#' 
#' @format GRangesList
#' @usage data(txdb_mouse, package="crisprDesignData")
"txdb_mouse"

#' GRangesList representing Ensembl gene model coordinates for crab-eating macaque
#'
#' GRangesList representing Ensembl gene model coordinates for the crab-eating macaque
#' in macFas6 coordinates (Ensembl release 108)
#' 
#' @format GRangesList
#' @usage data(txdb_cyno, package="crisprDesignData")
"txdb_cyno"




#' GRangesList representing Ensembl gene model coordinates for human
#'
#' GRanges representing Ensembl-based TSS coordinates for human in 
#' hg38/GRCh38 coordinates (Ensembl release 104)
#' 
#' @format GRanges
#' @usage data(tss_human, package="crisprDesignData")
"tss_human"


#' GRanges representing Ensembl-based TSS coordinates for mouse
#'
#' GRanges representing Ensembl-based TSS coordinates for mouse in mm10/GRCm38
#' coordinates (Ensembl release 102)
#' 
#' @format GRanges
#' @usage data(tss_mouse, package="crisprDesignData")
"tss_mouse"



#' Genomic coordinates of FANTOM5 CAGE Peaks in hg38 coordinates
#'
#' Genomic coordinates of FANTOM5 CAGE Peaks in hg38 coordinates.
#' 
#' @format GRanges
#' @usage data(humanCagePeaks, package="crisprDesignData")
"humanCagePeaks"



#' Genomic coordinates of FANTOM5 CAGE Peaks in mm10 coordinates
#'
#' Genomic coordinates of FANTOM5 CAGE Peaks in mm10 coordinates.
#' 
#' @format GRanges
#' @usage data(mouseCagePeaks, package="crisprDesignData")
"mouseCagePeaks"




#' DNAStringSet of mRNA sequences of human transcripts (Ensembl)
#'
#' DNAStringSet of mRNA sequences of human transcripts
#' (Ensembl release 104) extracted from the txdb_human object.
#' 
#' @format DNAStringSet
#' @usage data(mrnasHuman, package="crisprDesignData")
"mrnasHuman"


#' DNAStringSet of mRNA sequences of mouse transcripts (Ensembl)
#'
#' DNAStringSet of mRNA sequences of mouse transcripts
#' (Ensembl release 102) extracted from the txdb_human object.
#' 
#' @format DNAStringSet
#' @usage data(mrnasMouse, package="crisprDesignData")
"mrnasMouse"



#' APPRIS annotation for human genes (Ensembl release 104)
#'
#' data.frame containing APPRIS annotation for human genes
#' using the APPRIS release 47 and the Ensembl release 104
#' gene annotation.
#' 
#' @format data.frame
#' @usage data(apprisHuman, package="crisprDesignData")
"apprisHuman"


#' APPRIS annotation for mouse genes (Ensembl release 102)
#'
#' data.frame containing APPRIS annotation for mouse genes
#' using the APPRIS release 47 and the Ensembl release 102
#' gene annotation.
#' 
#' @format data.frame
#' @usage data(apprisMouse, package="crisprDesignData")
"apprisMouse"


#' Canonical Ensembl transcripts for human genes (Ensembl release 104)
#'
#' data.frame containing canonical Ensembl transcripts for human genes.
#' First column contains transcript ID, second column contains gene ID.
#' 
#' @format data.frame
#' @usage data(canonicalHuman, package="crisprDesignData")
"canonicalHuman"


#' Canonical Ensembl transcripts for mouse genes (Ensembl release 102)
#'
#' data.frame containing canonical Ensembl transcripts for mouse genes.
#' First column contains transcript ID, second column contains gene ID.
#' 
#' @format data.frame
#' @usage data(canonicalMouse, package="crisprDesignData")
"canonicalMouse"

#' Canonical Ensembl transcripts for crab-eating macaque genes (Ensembl release 108)
#'
#' data.frame containing canonical Ensembl transcripts for crab-eating macaque 
#' genes. First column contains transcript ID, second column contains gene ID.
#' 
#' @format data.frame
#' @usage data(canonicalCyno, package="crisprDesignData")
"canonicalCyno"



#' Table of Pfam doains for human transcripts
#'
#' DataFrame containing Pfam domain information for human Ensembl transcripts.
#' (Ensembl release 104)).
#' 
#' @format DataFrame
#' @usage data(pfamTableHuman, package="crisprDesignData")
"pfamTableHuman"

#' Table of Pfam doains for mouse transcripts
#'
#' DataFrame containing Pfam domain information for mouse Ensembl transcripts.
#' (Ensembl release 102)).
#' 
#' @format DataFrame
#' @usage data(pfamTableMouse, package="crisprDesignData")
"pfamTableMouse"




