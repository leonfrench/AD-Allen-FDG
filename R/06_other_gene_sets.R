library(here)
library(readr)
library(metap)
library(dplyr)
library(tidyr)
library(magrittr)
library(homologene)
library(tmod)

source(here("R", "GeneSetBuilders.R"))

#single runs for testing
#run_single_folder <- "donors_12876_14380_15496_15697_9861"
#run_single_folder <- "donors_10021"

use_rnaseq <- FALSE
if (use_rnaseq) {
  expression_dir <- 'rnaseq'
} else {
  expression_dir <- 'microarray'
}


prefixed_lists <- c("GO", "Darmanis")

for(result_single_folder in list.files(here("results", expression_dir), pattern = "_NiftiValue", full.names = F)) {
  print(result_single_folder)
  for(run_single_folder in list.files(here("results", expression_dir, result_single_folder), pattern = "donor", full.names = F)) {
    print(run_single_folder)
    
    run_gene_stats <- read_csv(here("results", expression_dir, result_single_folder, run_single_folder, "geneStatistics.csv"))
    run_gene_stats %<>% arrange(rank)
    sortedGenes <- run_gene_stats$gene_symbol
    
    for (prefix in prefixed_lists) {
      isMouse = prefix == "NeuroExpresso.Cortex"
      gene_sets <- loadFileSets(prefix, isMouse)
      filterGenes <- FALSE
      result <- tbl_df(tmodUtest(c(sortedGenes), mset=gene_sets, qval = 1.01, filter = filterGenes))
      
      #double p-value because tmod doesn't 
      result %<>% mutate(P.Value=P.Value*2) 

      #adjust again
      result %<>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value)
      #write out file for that prefix table
      result %>% write_csv(here("results", expression_dir, result_single_folder, run_single_folder, paste0(prefix,"_custom_genelists_results.csv")))
    }
  }
}


