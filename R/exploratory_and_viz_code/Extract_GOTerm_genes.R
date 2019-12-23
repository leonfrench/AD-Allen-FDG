#this code is for visualizing our results for specific GO groups
library(here)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readr)

source(here("R", "GeneSetBuilders.R"))
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)

#get genes used (a bit on the slow side)
gene_universe <- read_csv(here("data/R_processed_expression_data/","microarray_expression_probe_averaged_qcNames.csv")) %>% 
  select(gene_symbol) %>% distinct() %>% .$gene_symbol

if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { 
} else {
  geneSetsGO <- loadGOSets(gene_universe)
}
go_name_to_genes <- function(go_name) {
  targetGroupID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == go_name)$ID
  unlist(geneSetsGO$MODULES2GENES[targetGroupID])
}
go_name_to_genes("HOPS complex")
length(go_name_to_genes("nuclear-transcribed mRNA catabolic process, nonsense-mediated decay"))
length(intersect(go_name_to_genes("nuclear-transcribed mRNA catabolic process, nonsense-mediated decay"), go_name_to_genes("SRP-dependent cotranslational protein targeting to membrane")))
setdiff(go_name_to_genes("nuclear-transcribed mRNA catabolic process, nonsense-mediated decay"), go_name_to_genes("SRP-dependent cotranslational protein targeting to membrane"))

read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% filter(gene_symbol %in% go_name_to_genes("HOPS complex"))
as.data.frame(read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% filter(gene_symbol %in% setdiff(go_name_to_genes("nuclear-transcribed mRNA catabolic process, nonsense-mediated decay"), go_name_to_genes("SRP-dependent cotranslational protein targeting to membrane"))))
read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% filter(gene_symbol %in% c("UBB", "UBC", "UBA52", "RPS27A"))

#ubiquitin ligase complex
read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% filter(gene_symbol %in% go_name_to_genes("ubiquitin ligase complex")) %>% write_csv("/Users/lfrench/Desktop/results/AlzAllenPET/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/ubiq.csv")

#exploratory
setdiff(go_name_to_genes("protein targeting to ER"), go_name_to_genes("SRP-dependent cotranslational protein targeting to membrane"))
setdiff(go_name_to_genes("SRP-dependent cotranslational protein targeting to membrane"), go_name_to_genes("protein targeting to ER"))

print_and_write <- function(targetGroup) {
  targetGroupID <- dplyr::filter(tbl_df(geneSetsGO$MODULES), Title == targetGroup)$ID
  print(targetGroupID)
  
  unlist(geneSetsGO$MODULES2GENES[targetGroupID])
}
print_and_write("SRP-dependent cotranslational protein targeting to membrane")

genes <- setdiff(read_csv(here("/data/gene_lists/GO.RPS and RPL genes.txt"), col_names = F)$X1,
        read_csv(here("/data/gene_lists/GO.SRPdependentTranslationalProteinTargetingMembrane.txt"), col_names = F)$X1) 
as.data.frame(x=genes) %>% write.table(here("/data/gene_lists/GO.RPLorSminusSRP.txt"), row.names = F, col.names = F, quote = F)
  
