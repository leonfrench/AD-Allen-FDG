library(here)
library(readr)
library(tidyr)
library(magrittr)
library(homologene)
library(GO.db) #Leon is using April 24 source date (GO.db reports this info)
library(tmod)

source(here("R", "GeneSetBuilders.R"))
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)

#single runs for testing
#run_single_folder <- "donors_12876_14380_15496_15697_9861"
#run_single_folder <- "donors_10021"

use_rnaseq <- TRUE
for(use_rnaseq in c(TRUE, FALSE)) {
  if (use_rnaseq) {
    expression_dir <- 'rnaseq'
  } else {
    expression_dir <- 'microarray'
  }
  
  
  for(result_single_folder in list.files(here("results", expression_dir), pattern = "_NiftiValue", full.names = F)) {
    print(result_single_folder)
    for(run_single_folder in list.files(here("results", expression_dir, result_single_folder), pattern = "donor", full.names = F)) {
      print(run_single_folder)
      
      run_gene_stats <- read_csv(here("results", expression_dir, result_single_folder, run_single_folder, "geneStatistics.csv"))
      run_gene_stats %<>% arrange(rank)
      sortedGenes <- run_gene_stats$gene_symbol
      
      if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { 
      } else {
        geneSetsGO <- loadGOSets(sortedGenes)
      }
      
      filterGenes <- TRUE
      result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1.01, filter = filterGenes))
      result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
      
      #merge same GO groups using the genes in the set - these are sets that have the same set of genes
      result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsGO$MODULES2GENES[ID])), collapse = " "))
      result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
        summarize(MainTitle = first(Title), ID=first(ID), all_IDs=paste(ID, collapse=","), 
                  AUC = first(AUC), P.Value= first(P.Value), aspect= first(aspect), 
                  otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
      result %<>% ungroup() %>% dplyr::select(-genes)
      
      #double p-value because tmod doesn't 
      result %<>% mutate(P.Value=P.Value*2) 
      
      #adjust again now that identical groups have been merged
      result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value)
      result %<>% mutate(rank = rank(P.Value))
      result %<>% dplyr::select(MainTitle, ID, geneCount = N1, AUC, P.Value, adj.P.Value, everything()) 
      
      result %>% write_csv(here("results", expression_dir, result_single_folder, run_single_folder, "GO_results.csv"))
    }
  }
  
}
