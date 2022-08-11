crisprDesignData: example data for the crisprDesign ecosystem
================

-   [Overview](#overview)
-   [Installation](#installation)
    -   [Software requirements](#software-requirements)
        -   [OS Requirements](#os-requirements)
    -   [Installation](#installation-1)
        -   [Getting started](#getting-started)
-   [Datasets](#datasets)
-   [TxDb datasets](#txdb-datasets)
-   [TSS datasets](#tss-datasets)
-   [mRNA datasets](#mrna-datasets)
-   [Repeats datasets](#repeats-datasets)
-   [License](#license)
-   [Reproducibility](#reproducibility)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: July 16, 2022

# Overview

The `crisprDesignData` package provides ready-to-use annotation data
needed needed for the `crisprDesign` ecosystem, for both human and
human.

# Installation

## Software requirements

### OS Requirements

This package is supported for macOS, Linux and Windows machines. It was
developed and tested on R version 4.2.

## Installation

`crisprDesignData` can be installed by typing the following commands
inside of an R session:

``` r
install.packages("devtools")
devtools::install_github("Jfortin1/crisprDesignData")
```

### Getting started

`crisprDesignData` can be loaded into an R session in the usual way:

``` r
library(crisprDesignData)
```

# Datasets

| Object name       | Object class   | Version     | Description                                                     |
|-------------------|----------------|-------------|-----------------------------------------------------------------|
| `txdb_human`      | `GRangesList`  | Release 104 | Ensembl gene model for human (hg38/GRCh38)                      |
| `txdb_mouse`      | `GRangesList`  | Release 102 | Ensembl gene model for mouse (mm10/GRCm38)                      |
| `tss_human`       | `GRanges`      | Release 104 | Ensembl-based TSS coordinates for human (hg38/GRCh38)           |
| `tss_mouse`       | `GRanges`      | Release 102 | Ensembl-based TSS coordinates for human (mm10/GRCm38)           |
| `mrnasHuman`      | `DNAStringSet` | Release 104 | Ensembl-based mRNA nucleotide sequences for human (hg38/GRCh38) |
| `mrnasMouse`      | `DNAStringSet` | Release 102 | Ensembl-based mRNA nucleotide sequences for mouse (mm10/GRCm38) |
| `gr.repeats.hg38` | `GRanges`      |             | RepeatMasker data from UCSC genome browser (hg38/GRCh38)        |
| `gr.repeats.mm10` | `GRanges`      |             | RepeatMasker data from UCSC genome browser (mm10/GRCm38)        |

# TxDb datasets

The `txdb_human` and `txdb_mouse` objects are `GRangesList` representing
gene models for human and mouse, respectively, from Ensembl. They were
constructed using the function `getTxDb` in `crisprDesign`. See the
script `generateTxDbData.R`in the `inst` folder to see how to generate
such data for other organisms (internet connection needed).

Let’s look at the `txdb_human` object. We first load the data:

``` r
data(txdb_human, package="crisprDesignData")
```

We can look at metadata information about the gene model by using the
`metadata` function from the `S4Vectors` package:

``` r
head(S4Vectors::metadata(txdb_human))
```

    ##                 name                    value
    ## 1            Db type                     TxDb
    ## 2 Supporting package          GenomicFeatures
    ## 3        Data source                  Ensembl
    ## 4           Organism             Homo sapiens
    ## 5    Ensembl release                      104
    ## 6   Ensembl database homo_sapiens_core_104_38

The object is a `GRangesList` with 7 elements that contain genomic
coordinates for different levels of the gene model:

``` r
names(txdb_human)
```

    ## [1] "transcripts" "exons"       "cds"         "fiveUTRs"    "threeUTRs"  
    ## [6] "introns"     "tss"

As an example, let’s look at the `GRanges` containing genomic
coordinates for all exons represented in the gene model:

``` r
txdb_human$exons
```

    ## GRanges object with 796644 ranges and 14 metadata columns:
    ##     seqnames      ranges strand |           tx_id         gene_id
    ##        <Rle>   <IRanges>  <Rle> |     <character>     <character>
    ##         chr1 11869-12227      + | ENST00000456328 ENSG00000223972
    ##         chr1 12613-12721      + | ENST00000456328 ENSG00000223972
    ##         chr1 13221-14409      + | ENST00000456328 ENSG00000223972
    ##         chr1 12010-12057      + | ENST00000450305 ENSG00000223972
    ##         chr1 12179-12227      + | ENST00000450305 ENSG00000223972
    ##   .      ...         ...    ... .             ...             ...
    ##         chrM   5826-5891      - | ENST00000387409 ENSG00000210144
    ##         chrM   7446-7514      - | ENST00000387416 ENSG00000210151
    ##         chrM 14149-14673      - | ENST00000361681 ENSG00000198695
    ##         chrM 14674-14742      - | ENST00000387459 ENSG00000210194
    ##         chrM 15956-16023      - | ENST00000387461 ENSG00000210196
    ##          protein_id                tx_type gene_symbol         exon_id
    ##         <character>            <character> <character>     <character>
    ##                <NA>   processed_transcript     DDX11L1 ENSE00002234944
    ##                <NA>   processed_transcript     DDX11L1 ENSE00003582793
    ##                <NA>   processed_transcript     DDX11L1 ENSE00002312635
    ##                <NA> transcribed_unproces..     DDX11L1 ENSE00001948541
    ##                <NA> transcribed_unproces..     DDX11L1 ENSE00001671638
    ##   .             ...                    ...         ...             ...
    ##                <NA>                Mt_tRNA       MT-TY ENSE00001544488
    ##                <NA>                Mt_tRNA      MT-TS1 ENSE00001544487
    ##     ENSP00000354665         protein_coding      MT-ND6 ENSE00001434974
    ##                <NA>                Mt_tRNA       MT-TE ENSE00001544476
    ##                <NA>                Mt_tRNA       MT-TP ENSE00001544473
    ##     exon_rank cds_start   cds_end  tx_start    tx_end   cds_len exon_start
    ##     <integer> <integer> <integer> <integer> <integer> <integer>  <integer>
    ##             1      <NA>      <NA>     11869     14409         0       <NA>
    ##             2      <NA>      <NA>     11869     14409         0       <NA>
    ##             3      <NA>      <NA>     11869     14409         0       <NA>
    ##             1      <NA>      <NA>     12010     13670         0       <NA>
    ##             2      <NA>      <NA>     12010     13670         0       <NA>
    ##   .       ...       ...       ...       ...       ...       ...        ...
    ##             1      <NA>      <NA>      5826      5891         0       <NA>
    ##             1      <NA>      <NA>      7446      7514         0       <NA>
    ##             1     14149     14673     14149     14673       525       <NA>
    ##             1      <NA>      <NA>     14674     14742         0       <NA>
    ##             1      <NA>      <NA>     15956     16023         0       <NA>
    ##      exon_end
    ##     <integer>
    ##          <NA>
    ##          <NA>
    ##          <NA>
    ##          <NA>
    ##          <NA>
    ##   .       ...
    ##          <NA>
    ##          <NA>
    ##          <NA>
    ##          <NA>
    ##          <NA>
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

The function `queryTxObject` in `crisprDesign` is a user-friendly
function to work with such objects, for instance once can return the CDS
coordinates for the KRAS transcripts using the following lines of code:

``` r
library(crisprDesign)
cds <- queryTxObject(txdb_human,
                     featureType="cds",
                     queryColumn="gene_symbol",
                     queryValue="KRAS")
head(cds)
```

    ## GRanges object with 6 ranges and 14 metadata columns:
    ##            seqnames            ranges strand |           tx_id         gene_id
    ##               <Rle>         <IRanges>  <Rle> |     <character>     <character>
    ##   region_1    chr12 25245274-25245384      - | ENST00000256078 ENSG00000133703
    ##   region_2    chr12 25227234-25227412      - | ENST00000256078 ENSG00000133703
    ##   region_3    chr12 25225614-25225773      - | ENST00000256078 ENSG00000133703
    ##   region_4    chr12 25215441-25215560      - | ENST00000256078 ENSG00000133703
    ##   region_5    chr12 25245274-25245384      - | ENST00000311936 ENSG00000133703
    ##   region_6    chr12 25227234-25227412      - | ENST00000311936 ENSG00000133703
    ##                 protein_id        tx_type gene_symbol         exon_id exon_rank
    ##                <character>    <character> <character>     <character> <integer>
    ##   region_1 ENSP00000256078 protein_coding        KRAS ENSE00000936617         2
    ##   region_2 ENSP00000256078 protein_coding        KRAS ENSE00001719809         3
    ##   region_3 ENSP00000256078 protein_coding        KRAS ENSE00001644818         4
    ##   region_4 ENSP00000256078 protein_coding        KRAS ENSE00001189807         5
    ##   region_5 ENSP00000256078 protein_coding        KRAS ENSE00000936617         2
    ##   region_6 ENSP00000256078 protein_coding        KRAS ENSE00001719809         3
    ##            cds_start   cds_end  tx_start    tx_end   cds_len exon_start
    ##            <integer> <integer> <integer> <integer> <integer>  <integer>
    ##   region_1      <NA>      <NA>  25205246  25250929       570   25245274
    ##   region_2      <NA>      <NA>  25205246  25250929       570   25227234
    ##   region_3      <NA>      <NA>  25205246  25250929       570   25225614
    ##   region_4      <NA>      <NA>  25205246  25250929       570   25215437
    ##   region_5      <NA>      <NA>  25205246  25250929       567   25245274
    ##   region_6      <NA>      <NA>  25205246  25250929       567   25227234
    ##             exon_end
    ##            <integer>
    ##   region_1  25245395
    ##   region_2  25227412
    ##   region_3  25225773
    ##   region_4  25215560
    ##   region_5  25245395
    ##   region_6  25227412
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

# TSS datasets

The `tss_human` and `tss_mouse` objects are `GRanges` representing the
transcription starting sites (TSSs) coordinates for human and mouse,
respectively. The coordinates were extracted from the transcripts stored
in the Ensembl-based models `txdb_human` and `txdb_mouse` using the
function `getTssObjectFromTxObject` from `crisprDesign`. See the script
`generateTssObjects.R`in the `inst` folder to see how to generate such
data.

Let’s take a look at `tss_human`:

``` r
data(tss_human, package="crisprDesignData")
head(tss_human)
```

    ## GRanges object with 6 ranges and 5 metadata columns:
    ##         seqnames    ranges strand |           tx_id         gene_id gene_symbol
    ##            <Rle> <IRanges>  <Rle> |     <character>     <character> <character>
    ##   11402     chr1     65419      + | ENST00000641515 ENSG00000186092       OR4F5
    ##   11442     chr1    923923      + | ENST00000616016 ENSG00000187634      SAMD11
    ##   11444     chr1    925731      + | ENST00000342066 ENSG00000187634      SAMD11
    ##   11445     chr1    960584      + | ENST00000338591 ENSG00000187961      KLHL17
    ##   11446     chr1    960639      + | ENST00000622660 ENSG00000187961      KLHL17
    ##   11447     chr1    966482      + | ENST00000379410 ENSG00000187583     PLEKHN1
    ##                promoter                     ID
    ##             <character>            <character>
    ##   11402 ENST00000641515  OR4F5_ENST00000641515
    ##   11442 ENST00000616016 SAMD11_ENST00000616016
    ##   11444 ENST00000342066 SAMD11_ENST00000342066
    ##   11445 ENST00000338591 KLHL17_ENST00000338591
    ##   11446 ENST00000622660 KLHL17_ENST00000622660
    ##   11447 ENST00000379410 PLEKHN1_ENST00000379..
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

The function `queryTss` in `crisprDesign` is a user-friendly function to
work with such objects, accepting an argument called `tss_window` to
specify a number of nucleotides upstream and downstream of the TSS. This
is particularly useful to return genomic regions to target for CRISPRa
and CRISPRi.

For instance, if we want to target the region 500 nucleotides upstream
of any of the KRAS TSSs, one can use the following lines of code:

``` r
library(crisprDesign)
tss <- queryTss(tss_human,
                queryColumn="gene_symbol",
                queryValue="KRAS",
                tss_window=c(-500,0))
head(tss)
```

    ## GRanges object with 3 ranges and 5 metadata columns:
    ##            seqnames            ranges strand |           tx_id         gene_id
    ##               <Rle>         <IRanges>  <Rle> |     <character>     <character>
    ##   region_1    chr12 25250752-25251251      - | ENST00000256078 ENSG00000133703
    ##   region_2    chr12 25250752-25251251      - | ENST00000557334 ENSG00000133703
    ##   region_3    chr12 25250765-25251264      - | ENST00000556131 ENSG00000133703
    ##            gene_symbol        promoter                   ID
    ##            <character>     <character>          <character>
    ##   region_1        KRAS ENST00000256078 KRAS_ENST00000256078
    ##   region_2        KRAS ENST00000557334 KRAS_ENST00000557334
    ##   region_3        KRAS ENST00000556131 KRAS_ENST00000556131
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

# mRNA datasets

The `mrnasHuman` and `mrnasMouse` objects are `DNAStringSet` storing the
nucleotide sequence of mRNAs derived from the `txdb_human` and
`txdb_mouse` gene models, respectively. It was obtained using the
function `getMrnaSequences` from `crisprDesign`. See the script
`generateMrnaData.R`in the `inst` folder to see how to generate such
data. The names of the `DNAStringSet` are Ensembl transcript IDs:

``` r
data(mrnasHuman, package="crisprDesignData")
data(mrnasMouse, package="crisprDesignData")
head(mrnasHuman)
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]  1032 CTGCTGCTGCTGCGCCCCATCCC...TAATAAATTTGCTGTGGTTTGTA ENST00000000233
    ## [2]  2450 AGAGTGGGGCACAGCGAGGCGCT...GATTAAAAAACAAACAAAACATA ENST00000000412
    ## [3]  2274 GTCAGCTGGAGGAAGCGGAGTAG...ATATATAATACCGAGCTCAAAAA ENST00000000442
    ## [4]  3715 CCTACCCCAGCTCTCGCGCCGCG...CTAGTGAGGATGTTTTGTTAAAA ENST00000001008
    ## [5]  4556 ACAGCCAATCCCCCGAGCGGCCG...TAGAATAAACCGTGGGGACCCGC ENST00000001146
    ## [6]  2184 GTCTAAGCCCCAGCTCCTGGCGG...ACGAGTAATTTCATAGAAACGAA ENST00000002125

``` r
head(mrnasMouse)
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]  3262 CACACATCCGGTTCTTCCGGGAG...TTCTTCACTTCTTTGTCTCTGCA ENSMUST00000000001
    ## [2]   902 GTCAGTGCACAACTGCCAACTGG...TTAAATAAATTTATTTTACTTGC ENSMUST00000000003
    ## [3]  2574 GGTCCGTGTGCCACCTTTTCCCT...AATATACATATCACTCTAGAAAA ENSMUST00000000010
    ## [4]  2143 TGGAAACACATTCAAATAATGTG...AAAATGTTTGATGTTTTATCCCC ENSMUST00000000028
    ## [5]  3708 GGCACTGACCAGTTCGCAAACTG...AATAAAGCATTTAAAATACTATT ENSMUST00000000033
    ## [6]  1190 TCTCTTCAGCAGAAGACACCACT...AAGTTACTGGATTGCCTCAGTTA ENSMUST00000000049

Those objects are particularly useful for gRNA design for RNA-targeting
nucleases such as RfxCas13d (CasRx).

# Repeats datasets

The objects `gr.repeats.hg38` and `gr.repeats.mm10` objects are
`GRanges` representing the genomic coordinates of repeat elements in the
human and mouse genomes, as defined by the RepeatMasker tracks in the
UCSC genome browser.

Let’s look at the repeats elements in the human genome:

``` r
data(gr.repeats.hg38, package="crisprDesignData")
head(gr.repeats.hg38)
```

    ## GRanges object with 6 ranges and 2 metadata columns:
    ##       seqnames            ranges strand |        type     score
    ##          <Rle>         <IRanges>  <Rle> | <character> <numeric>
    ##   [1]     chr1 67108753-67109046      + |        L1P5      1892
    ##   [2]     chr1   8388315-8388618      - |        AluY      2582
    ##   [3]     chr1 25165803-25166380      + |       L1MB5      4085
    ##   [4]     chr1 33554185-33554483      - |       AluSc      2285
    ##   [5]     chr1 41942894-41943205      - |        AluY      2451
    ##   [6]     chr1 50331336-50332274      + |        HAL1      1587
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

# License

The package is licensed under the MIT license.

# Reproducibility

``` r
sessionInfo()
```

    ## R Under development (unstable) (2022-03-21 r81954)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] crisprDesign_0.99.99    crisprBase_1.1.2        GenomicRanges_1.47.6   
    ## [4] GenomeInfoDb_1.31.6     IRanges_2.29.1          S4Vectors_0.33.11      
    ## [7] BiocGenerics_0.41.2     crisprDesignData_0.99.8
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] bitops_1.0-7                  matrixStats_0.61.0           
    ##   [3] bit64_4.0.5                   filelock_1.0.2               
    ##   [5] progress_1.2.2                httr_1.4.2                   
    ##   [7] tools_4.2.0                   utf8_1.2.2                   
    ##   [9] R6_2.5.1                      DBI_1.1.2                    
    ##  [11] tidyselect_1.1.2              prettyunits_1.1.1            
    ##  [13] bit_4.0.4                     curl_4.3.2                   
    ##  [15] compiler_4.2.0                crisprBowtie_1.1.1           
    ##  [17] cli_3.3.0                     Biobase_2.55.0               
    ##  [19] basilisk.utils_1.9.1          crisprScoreData_1.1.3        
    ##  [21] xml2_1.3.3                    DelayedArray_0.21.2          
    ##  [23] rtracklayer_1.55.4            randomForest_4.7-1           
    ##  [25] readr_2.1.2                   rappdirs_0.3.3               
    ##  [27] stringr_1.4.0                 digest_0.6.29                
    ##  [29] Rsamtools_2.11.0              rmarkdown_2.13               
    ##  [31] crisprScore_1.1.9             basilisk_1.9.2               
    ##  [33] XVector_0.35.0                pkgconfig_2.0.3              
    ##  [35] htmltools_0.5.2               MatrixGenerics_1.7.0         
    ##  [37] dbplyr_2.1.1                  fastmap_1.1.0                
    ##  [39] BSgenome_1.63.5               rlang_1.0.2                  
    ##  [41] rstudioapi_0.13               RSQLite_2.2.12               
    ##  [43] shiny_1.7.1                   BiocIO_1.5.0                 
    ##  [45] generics_0.1.2                jsonlite_1.8.0               
    ##  [47] BiocParallel_1.29.18          dplyr_1.0.8                  
    ##  [49] VariantAnnotation_1.41.3      RCurl_1.98-1.6               
    ##  [51] magrittr_2.0.2                GenomeInfoDbData_1.2.7       
    ##  [53] Matrix_1.4-0                  Rcpp_1.0.8.3                 
    ##  [55] fansi_1.0.2                   reticulate_1.24              
    ##  [57] Rbowtie_1.35.0                lifecycle_1.0.1              
    ##  [59] stringi_1.7.6                 yaml_2.3.5                   
    ##  [61] SummarizedExperiment_1.25.3   zlibbioc_1.41.0              
    ##  [63] AnnotationHub_3.3.9           BiocFileCache_2.3.4          
    ##  [65] grid_4.2.0                    blob_1.2.2                   
    ##  [67] promises_1.2.0.1              parallel_4.2.0               
    ##  [69] ExperimentHub_2.3.5           crayon_1.5.0                 
    ##  [71] dir.expiry_1.3.0              lattice_0.20-45              
    ##  [73] Biostrings_2.63.2             GenomicFeatures_1.47.13      
    ##  [75] hms_1.1.1                     KEGGREST_1.35.0              
    ##  [77] knitr_1.37                    pillar_1.7.0                 
    ##  [79] rjson_0.2.21                  biomaRt_2.51.3               
    ##  [81] BiocVersion_3.15.0            XML_3.99-0.9                 
    ##  [83] glue_1.6.2                    evaluate_0.15                
    ##  [85] BiocManager_1.30.16           httpuv_1.6.5                 
    ##  [87] png_0.1-7                     vctrs_0.3.8                  
    ##  [89] tzdb_0.2.0                    purrr_0.3.4                  
    ##  [91] assertthat_0.2.1              cachem_1.0.6                 
    ##  [93] xfun_0.30                     mime_0.12                    
    ##  [95] xtable_1.8-4                  restfulr_0.0.13              
    ##  [97] later_1.3.0                   tibble_3.1.6                 
    ##  [99] GenomicAlignments_1.31.2      AnnotationDbi_1.57.1         
    ## [101] memoise_2.0.1                 interactiveDisplayBase_1.33.0
    ## [103] ellipsis_0.3.2
