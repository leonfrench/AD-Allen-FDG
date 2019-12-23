library(here)
library(readr)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(tidyr)
library(magrittr)

exp_matrix <- 'microarray'
all6_geneStatistics <- read_csv(here("results", exp_matrix, "Allen_Schwartz_NiftiValue", "donors_9861_10021_12876_14380_15496_15697", "geneStatistics.csv"))
all6_geneStatistics %<>% mutate(isSignificant = (metaP.pos.adj < 0.025) | (metaP.neg.adj < 0.025))
sum(all6_geneStatistics$isSignificant)
range((all6_geneStatistics %>% filter(isSignificant))$average_log2_fold_change)
2**range((all6_geneStatistics %>% filter(isSignificant))$average_log2_fold_change)

donor_10021_geneStatistics <- read_csv(here("results", exp_matrix, "Allen_Schwartz_NiftiValue", "donors_10021", "geneStatistics.csv"))
donor_10021_geneStatistics %<>% mutate(isSignificant = (metaP.pos.adj < 0.025) | (metaP.neg.adj < 0.025))
sum(donor_10021_geneStatistics$isSignificant)
range((donor_10021_geneStatistics %>% filter(isSignificant))$average_log2_fold_change)
2**range((donor_10021_geneStatistics %>% filter(isSignificant))$average_log2_fold_change)


donor_12876_geneStatistics <- read_csv(here("results", exp_matrix, "Allen_Schwartz_NiftiValue", "donors_12876", "geneStatistics.csv"))
donor_12876_geneStatistics %<>% mutate(isSignificant = (metaP.pos.adj < 0.025) | (metaP.neg.adj < 0.025))
sum(donor_12876_geneStatistics$isSignificant)
range((donor_12876_geneStatistics %>% filter(isSignificant))$average_log2_fold_change)
2**range((donor_12876_geneStatistics %>% filter(isSignificant))$average_log2_fold_change)


#iterate all folders to get number of significant genes after correction
result_single_folder <- "Allen_Schwartz_NiftiValue"
sigGenes <- as_tibble()

for(run_single_folder in list.files(here("results", exp_matrix, result_single_folder), pattern = "donor", full.names = F)) {
  print(run_single_folder)
  geneStatistics <- read_csv(here("results", exp_matrix, result_single_folder, run_single_folder, "geneStatistics.csv"))
  geneStatistics %<>% mutate(isSignificant = (metaP.pos.adj < 0.025) | (metaP.neg.adj < 0.025))
  
  geneStatistics %<>% mutate(isSignificantPos = (metaP.pos.adj < 0.025))
  geneStatistics %<>% mutate(isSignificantNeg = (metaP.neg.adj < 0.025))
  
  sigGenes <- bind_rows(sigGenes,data.frame(donor_set = run_single_folder, num_sig= sum(geneStatistics$isSignificant), num_neg_sig= sum(geneStatistics$isSignificantNeg), num_ps_sig= sum(geneStatistics$isSignificantPos)))
}
sigGenes


