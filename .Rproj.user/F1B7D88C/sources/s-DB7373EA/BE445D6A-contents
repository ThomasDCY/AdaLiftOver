## AdaLiftOver
**AdaLiftOver** is a handy R package for adaptively identifying orthologous regions across different species. For given query genomic regions, AdaLiftOver incorporates epigenomic signals as well as sequence grammars to priorize and pinpoint candidate target genomic regions.


## Installation

**AdaLiftOver** can be downloaded and installed in R by: 

```r
## install.packages("devtools")
devtools::install_github("ThomasDCY/AdaLiftOver")
```

If the installation fails, make sure you can install the following R packages:

```r
## data.table
install.packages("data.table")

## Matrix
install.packages("Matrix")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## motifmatchr
BiocManager::install("motifmatchr")

## rtracklayer
BiocManager::install("rtracklayer")

## GenomicRanges
BiocManager::install("GenomicRanges")
```

## A quick start

Download the UCSC chain file from [mm10.hg38.rbest.chain.gz](http://hgdownload.cse.ucsc.edu/goldenpath/mm10/vsHg38/reciprocalBest/mm10.hg38.rbest.chain.gz) and unzip it.

```r
library(AdaLiftOver)

data("data_example")
data("epigenome_mm10")
data("epigenome_hg38")

## load the UCSC chain file
# chain <- rtracklayer::import.chain("mm10.hg38.rbest.chain")

## map the query regions
# gr_list <- adaptive_liftover(gr, chain)

## compute epigenome signal similarity
gr_list <- compute_similarity_epigenome(gr, gr_list, epigenome_mm10, epigenome_hg38)

## compute sequence grammar similarity
data("jaspar_pfm_list")
gr_list <- compute_similarity_grammar(gr, gr_list, "mm10", "hg38", jaspar_pfm_list)

## filter target candidate regions
gr_list_filter <- gr_candidate_filter(
    gr_list,
    best_k = 1L,
    threshold = 0.5
)
```


See the vignette for more information!

### Reference

**C. Dong**, and **S. Keles**, "AdaLiftOver: High-resolution identification of orthologous regulatory elements with adaptive liftOver".
