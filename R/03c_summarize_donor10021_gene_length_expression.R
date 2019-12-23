library(here)
library(readr)
library(metap)
#detach("package:dplyr", unload=TRUE)
library(dplyr)
library(tidyr)
library(corrr)
library(Hmisc)
library(tidyr)
library(tibble)
library(magrittr)
library(homologene)
library(biomaRt)


exp_matrix <- 'microarray'

##Stats for donor10021##
#number of over and under expression of genes after adjusted p<0.025 (same number as in script 03b)
donor_10021_geneStatistics <- read_csv(here("results", exp_matrix, "Allen_Schwartz_NiftiValue", "donors_10021", "geneStatistics.csv"))
donor_10021_geneStatistics %<>% mutate(isSignificant = (metaP.pos.adj < 0.025) | (metaP.neg.adj < 0.025))
D10021_OverExp <- filter(donor_10021_geneStatistics, isSignificant=="TRUE" & medianDirection == 1)
print (paste0("Number of overexpressed genes for donor10021  ", length(unique(D10021_OverExp$gene_symbol))))
D10021_UnderExp <- filter(donor_10021_geneStatistics, isSignificant=="TRUE" & medianDirection == -1) 
print (paste0("Number of underexpressed genes for donor10021  ", length(unique(D10021_UnderExp$gene_symbol))))

#find length of all gene symbols 
all_genelist <-  donor_10021_geneStatistics%>%distinct(gene_symbol)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords <-getBM(attributes=c("hgnc_symbol","start_position","end_position"), filters="hgnc_symbol", values=all_genelist$gene_symbol, mart=human)
gene_coords$size <- gene_coords$end_position - gene_coords$start_position
#merge gene size data with geneStatistics
D10021_geneSize <- merge(donor_10021_geneStatistics, gene_coords, by.x="gene_symbol", by.y="hgnc_symbol")

              
#wilcox test of rank and size between over and under expression of genes and median of size:
D10021_OverExp_size <- filter(D10021_geneSize, isSignificant=="TRUE" & medianDirection == 1)
D10021_UnderExp_size <- filter(D10021_geneSize, isSignificant=="TRUE" & medianDirection == -1) 
WilCox_Test_D10021 <- wilcox.test(D10021_OverExp_size$size, D10021_UnderExp_size$size)
WilCox_Test_D10021
D10021_OverExp_size_median <- median(D10021_OverExp_size$size)
print (paste0("Median length of overrexpressed genes ", D10021_OverExp_size_median, "bp")) 
D10021_UnderExp_size_median <- median(D10021_UnderExp_size$size)
print (paste0("Median length of underexpressed genes ", D10021_UnderExp_size_median, "bp")) 
