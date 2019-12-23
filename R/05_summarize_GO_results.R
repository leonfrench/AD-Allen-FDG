library(here)
library(readr)
library(tidyr)
library(magrittr)
detach("package:XML", unload = TRUE) #handles count()/n() conflict (sometimes)
detach("package:dplyr", unload = TRUE)
library(dplyr)

use_rnaseq <- TRUE
if (use_rnaseq) {
  expression_dir <- 'rnaseq'
} else {
  expression_dir <- 'microarray'
}

allGOResults <- as_tibble()
for(result_single_folder in list.files(here("results", expression_dir), pattern = "_NiftiValue", full.names = F)) {
  print(result_single_folder)
  for(run_single_folder in list.files(here("results", expression_dir, result_single_folder), pattern = "donor", full.names = F)) {
    print(run_single_folder)
    single_GO <- read_csv(here("results", expression_dir, result_single_folder, run_single_folder, "GO_results.csv"))
    single_GO %<>% mutate(donor_set = run_single_folder, result_set = result_single_folder)
    allGOResults <- bind_rows(allGOResults, single_GO)
  }
}

allGOResults %<>% separate(result_set, into = c("mapping", "PET_MAP"), remove = FALSE, extra = 'drop', sep="_")
allGOResults %>% dplyr::select(mapping, PET_MAP) %>% distinct()

#focus on one MRI mapping and target map
allGOResults %<>% filter(mapping == "Allen", PET_MAP == "Schwartz")

#signifcant group counts
allGOResults %<>% mutate(isSignificant = adj.P.Value < 0.05)
allGOResults %<>% mutate(isSignificantPos = (adj.P.Value < 0.05 & AUC > 0.5))
allGOResults %<>% mutate(isSignificantNeg = (adj.P.Value < 0.05 & AUC < 0.5))

allGOResults %>% filter(mapping == "Allen", PET_MAP == "Schwartz") %>% group_by(donor_set, isSignificant) %>% 
  filter(isSignificantPos | isSignificantNeg) %>% summarize(n = length(unique(ID)))

allGOResults %>% filter(mapping == "Allen", PET_MAP == "Schwartz") %>% group_by(donor_set, isSignificantPos, isSignificantNeg) %>% 
  filter(isSignificantPos | isSignificantNeg) %>% summarize(n = length(unique(ID)))

#############

#make top 10 tables
for(donor_set in unique(allGOResults$donor_set)) {
  results <- read_csv(here("results", expression_dir, "Allen_Schwartz_NiftiValue", donor_set, "GO_results.csv"))
  results %<>% mutate(directioned_rank = rank(-1*(sign(AUC-.5) * (max(rank) - rank))))
  
  results %>% arrange(directioned_rank) %>% head(10) %>% dplyr::select(Name=MainTitle, Genes = geneCount, ID, AUC, adj.P.Value) %>% 
    mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>% 
    dplyr::rename(`p-valueFDR` = adj.P.Value) %>%
    write_csv(here("results", expression_dir, "Allen_Schwartz_NiftiValue", donor_set, "GO_results.top10.csv"))
  
  results %>% arrange(-directioned_rank) %>% head(10) %>% dplyr::select(Name=MainTitle, Genes = geneCount, ID, AUC, adj.P.Value) %>% 
    mutate(AUC = signif(AUC, digits=3), adj.P.Value = signif(adj.P.Value, digits=3)) %>% 
    dplyr::rename(`p-valueFDR` = adj.P.Value) %>%
    write_csv(here("results", expression_dir, "Allen_Schwartz_NiftiValue", donor_set, "GO_results.bottom10.csv"))
}
  

#look at differences in GO group rankings
#sort with direction
allGOResults %<>% group_by(donor_set) %>% mutate(directioned_rank = rank(-1*(sign(AUC-.5) * (max(rank) - rank))))

if (use_rnaseq) {
  diff_rank <- allGOResults %>% select(ID, MainTitle, geneCount, donor_set, directioned_rank) %>% spread(donor_set, directioned_rank) %>% mutate(diff_rank = abs(donors_10021 - donors_9861))
  diff_rank %<>% arrange(-diff_rank)
  diff_rank %>% filter(donors_10021 < donors_9861) %>% arrange(-diff_rank)
  diff_rank %>% filter(donors_10021 > donors_9861) %>% arrange(-diff_rank)
  
} else {
  diff_rank <- allGOResults %>% select(ID, MainTitle, geneCount, donor_set, directioned_rank) %>% filter(donor_set %in% c("donors_12876_14380_15496_15697_9861", "donors_10021")) %>% spread(donor_set, directioned_rank) %>% mutate(diff_rank = abs(donors_10021 - donors_12876_14380_15496_15697_9861))
  diff_AUC <- allGOResults %>% select(ID, AUC, donor_set) %>% spread(donor_set, AUC)
  diff_padj <- allGOResults %>% select(ID, adj.P.Value, donor_set) %>% spread(donor_set, adj.P.Value) %>% rename(donors_10021.padj = donors_10021, donors_12876_14380_15496_15697_9861.padj = donors_12876_14380_15496_15697_9861)
  diff_rank <- inner_join(diff_rank, diff_AUC, suffix=c(".rank", ".AUC"), by="ID")
  diff_rank %<>% mutate(diff_AUC = abs(donors_10021.AUC -donors_12876_14380_15496_15697_9861.AUC))
  diff_rank <- inner_join(diff_rank, diff_padj, by="ID")
  diff_rank %<>% arrange(-diff_rank)
  diff_rank %>% filter(donors_10021.rank > donors_12876_14380_15496_15697_9861.rank) %>% arrange(-diff_rank)
  diff_rank %>% filter(donors_10021.rank < donors_12876_14380_15496_15697_9861.rank) %>% arrange(-diff_rank) 
  
  colnames(diff_rank) <- gsub("donors_9861_12876_14380_15496_15697", "remaining_five", colnames(diff_rank))
  
  #write it to both of the donor folders
  diff_rank %<>% select(Name=MainTitle, Genes = geneCount, ID, everything()) 
  for (donor_set in c("donors_12876_14380_15496_15697_9861", "donors_10021")) {
    diff_rank %>% write_csv(here("results", expression_dir, "Allen_Schwartz_NiftiValue", donor_set, "Supplement_rank_differences_comparison.csv"))
  }
}
