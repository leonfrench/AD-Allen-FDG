#Bioconductor libraries needs different code now
#this is not a complete list

# load csv using readr read_csv
if (! require(readr, quietly=TRUE)) {
  install.packages("readr")
  library(readr)
}

if (! require(data.table, quietly=TRUE)) {
  install.packages("data.table")
  library(data.table)
}
if (! require(metap, quietly=TRUE)) {
  install.packages("metap")
  library(metap)
}

if (! require(ggsignif, quietly=TRUE)) {
  install.packages("ggsignif")
  library(ggsignif)
}
if (! require(cowplot, quietly=TRUE)) {
  install.packages("cowplot")
  library(cowplot)
}
if (! require(readr, quietly=TRUE)) {
  install.packages("readr")
  library(readr)
}
if (! require(ggplot2, quietly=TRUE)) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (! require(tidyr, quietly=TRUE)) {
  install.packages("tidyr")
  library(tidyr)
}
if (! require(magrittr, quietly=TRUE)) {
  install.packages("magrittr")
  library(magrittr)
}

# need devtools to use install_github
if(! require(devtools, quiety=TRUE)) {
  install.packages("devtools")
  library(devtools)
}
if (! require(homologene, quietly=TRUE)) {
  install_github('oganm/homologene')
  library(homologene)
}

# libraries for entrez
if (! require(org.Hs.eg.db, quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("org.Hs.eg.db")
}
if (! require(AnnotationDbi, quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("AnnotationDbi")
}
if (! require(annotate, quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("annotate")
}
if (! require(GO.db, quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GO.db")
}
if (! require(tmod, quietly=TRUE)) {
  install.packages("tmod")
  library(tmod)
}
if (! require(metap, quietly=TRUE)) {
  install.packages("metap")
  library(metap)
}
if (! require(reshape2, quietly=TRUE)) {
  install.packages("reshape2")
  library(reshape2)
}

if (! require(markerGeneProfile, quietly=TRUE)) {
  devtools::install_github('oganm/markerGeneProfile')
}

f (! require(biomaRt, quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("biomaRt")
}

if (! require(dplyr, quietly=TRUE)) {
  install.packages("dplyr")
  library(dplyr)
}
