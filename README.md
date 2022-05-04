## AdaLiftOver
**AdaLiftOver** is a handy R package for adaptively identifying orthologous regions across different species. For each query region, AdaLiftOver outputs a scored and filtered list of candidate target regions that are most similar to the query region in terms of regulatory information.


## Installation

**AdaLiftOver** can be downloaded and installed in R by: 

```r
## install.packages("devtools")
devtools::install_github("ThomasDCY/AdaLiftOver", build_vignettes = TRUE)
```

If the installation fails, make sure you can install the following R packages:

```r
## data.table
install.packages("data.table")

## Matrix
install.packages("Matrix")

## PRROC
install.packages("PRROC")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## motifmatchr
## See also https://github.com/GreenleafLab/motifmatchr
BiocManager::install("motifmatchr")

## rtracklayer
BiocManager::install("rtracklayer")

## GenomicRanges
BiocManager::install("GenomicRanges")

## The BSgenome packages required to run the examples
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

## A quick start

Download the UCSC chain file from [mm10.hg38.rbest.chain.gz](http://hgdownload.cse.ucsc.edu/goldenpath/mm10/vsHg38/reciprocalBest/mm10.hg38.rbest.chain.gz) and unzip it.

```r
library(AdaLiftOver)

data("data_example")

## load the ENCODE repertoire
data("epigenome_mm10")
data("epigenome_hg38")

## load the UCSC chain file
chain <- rtracklayer::import.chain("mm10.hg38.rbest.chain")

## map the query regions
gr_list <- adaptive_liftover(gr, chain)

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

We might need to learn the parameters for a pair of matched epigenome datasets other than the ENCODE repertoire we provide, especially for model organisms other than mice. AdaLiftOver provides a handy training module to estimate the logistic regression parameters and suggest an optimal score threshold without filtering out too many candidate target regions.

```r
data("training_module_example")
training_module(gr_candidate, gr_true)
```


See the vignette for more information!

```r
browseVignettes("AdaLiftOver")
```

### Reference

**C. Dong**, and **S. Keles**, "AdaLiftOver: High-resolution identification of orthologous regulatory elements with adaptive liftOver".
