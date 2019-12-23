#test if there is overlap with GWAS hits - nothing of note
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

gene_sets <- loadFileSets("alz", FALSE)

universe <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% .$gene_symbol
#target_genes <- read_csv(here("/data/gene_lists/GO.SRPdependentTranslationalProteinTargetingMembrane.txt"), col_names = F) %>% .$X1

target_genes <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% 
  filter(metaP.pos.adj < 0.025) %>% .$gene_symbol


gene_stats <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv"))
#cell killing
gene_stats %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES["GO:0001906"]))
#chronic inflammatory response
gene_stats %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES["GO:0002544"])) %>% head(20)


for(target_set in rownames(gene_sets$MODULES)) {
  print(target_set)
  target_set_genes <- unlist(gene_sets$MODULES2GENES[target_set])
  print(intersect(target_set_genes,target_genes))
  print(length(target_genes))
  print(length(universe))
  print(length(intersect(universe, target_set_genes)))
  
}
# "CHRNE" "MS4A2" are in JasenEtAl.GWAS.TableS13, p = 0.27
