if (!require(BiocManager)){ install.packages('BiocManager', repos='https://www.stats.bris.ac.uk/R/')}

if (!require(DirichletMultinomial)){ BiocManager::install("DirichletMultinomial")}

if (!require(Biobase)){ BiocManager::install("Biobase")}

if (!require(optparse)){ install.packages('optparse', repos='https://www.stats.bris.ac.uk/R/')}

if (!require(tidyverse)){ install.packages('tidyverse', repos='https://www.stats.bris.ac.uk/R/')}

if (!require(devtools)){ install.packages("devtools", repos='https://www.stats.bris.ac.uk/R/')}

if (!require(lazyeval)){ install.packages('lazyeval', repos='https://www.stats.bris.ac.uk/R/')}

# library(devtools)
if (!requireNamespace("RNAseqProcessing", quietly = TRUE)){ devtools::install_github("RHReynolds/RNAseqProcessing")}

if (!requireNamespace("Biostrings", quietly = TRUE)){ BiocManager::install("Biostrings")}

if (!requireNamespace("tximport", quietly = TRUE)){ BiocManager::install("tximport")}

if (!requireNamespace("DESeq2", quietly = TRUE)){ BiocManager::install("DESeq2")}

if (!requireNamespace("leafcutter", quietly = TRUE)){ devtools::install_github("davidaknowles/leafcutter/leafcutter")}

if (!requireNamespace("rtracklayer", quietly = TRUE)){ BiocManager::install("rtracklayer")}

if (!requireNamespace("GenomicRanges", quietly = TRUE)){ BiocManager::install("GenomicRanges")}

if (!requireNamespace("megadepth", quietly = TRUE)) {
   BiocManager::install("megadepth") }

library("megadepth")

## Install Megadepth's pre-compiled binary on your system
install_megadepth()
