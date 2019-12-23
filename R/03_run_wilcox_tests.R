library(here)
library(readr)
library(metap)
library(dplyr)
library(tidyr)
library(magrittr)

#switch this flag to generate RNA sequence data
use_rnaseq <- FALSE

all_sample_annotations <- read_csv(here("/data/Combined_SampleAnnot_with_PET.csv"))

#restrict to in_cortex_excluding_piriform_hippocampus for this analysis
all_sample_annotations %<>% filter(in_cortex_excluding_piriform_hippocampus)

if (use_rnaseq) {
  full_expression_matrix <- read_csv(here("/data/R_processed_expression_data/rnaseq_averaged_qcNames.csv"))
  exp_matrix <- 'rnaseq'
  donor_groups <- c("10021", "9861")
} else {
  full_expression_matrix <- read_csv(here("/data/R_processed_expression_data/microarray_expression_probe_averaged_qcNames.csv"))
  exp_matrix <- 'microarray'
  donor_groups <- c(list(unique(all_sample_annotations$donor)), unique(all_sample_annotations$donor), list(setdiff(unique(all_sample_annotations$donor), "10021")))  
}


#for manually setting the donor groups use all
pet_map <- "Schwartz"
mni_map <- "Allen"
#groups for testing
#donor_group <- unique(all_sample_annotations$donor)
#donor_group <- c("10021")
#just the four for testing purposes 
donor_group <- c("14380", "15496", "15697", "9861")

for (donor_group in donor_groups) {
  donor_group <- sort(donor_group)
  print(paste("Donor group: ",paste("donors", paste(donor_group, collapse = "_"), sep="_")))

  pet_threshold <- 0.006

  #subset to donor group and available samples (only needed for RNAseq)
  sample_annotations <- all_sample_annotations %>% filter(donor %in% donor_group) %>% filter(well_id %in% full_expression_matrix$well_id)
  
  target_column <- paste(mni_map, pet_map, "NiftiValue", sep="_")
  print(paste("Using", nrow(sample_annotations), "data points"))
  sample_annotations %>% group_by(donor) %>% summarize(n=dplyr::n())
  
  #join to get expression
  donors_expression_matrix <- inner_join(full_expression_matrix, sample_annotations %<>% dplyr::select(well_id, target_column, donor))
  
  donors_expression_matrix %<>% mutate(inside = !!as.name(target_column) > pet_threshold) 
  
  #percent in/out of hypometabolism regions - should be moved to different code
  donors_expression_matrix %>% select(well_id, donor, inside) %>% distinct() %>% group_by(donor, inside) %>% summarize(n=dplyr::n()) %>% spread(inside,n) %>% mutate(percent = 100*`TRUE`/(`TRUE`+`FALSE`)) %>% arrange(percent)
  mean((donors_expression_matrix %>% select(well_id, donor, inside) %>% distinct() %>% group_by(donor, inside) %>% summarize(n=dplyr::n()) %>% spread(inside,n) %>% mutate(percent = 100*`TRUE`/(`TRUE`+`FALSE`)) %>% arrange(percent))$percent)
  
  wilcoxResultsTable <- donors_expression_matrix %>% group_by(gene_symbol, donor) %>%
    summarise(pvalue.pos = wilcox.test(expression~inside, alternative="less")$p.value,
              pvalue.neg = wilcox.test(expression~inside, alternative="greater")$p.value
    )
  wilcoxResultsTable %<>% mutate(direction = sign(pvalue.neg-pvalue.pos))
  
  message("Finding means and medians\n")
  averagesTable <- donors_expression_matrix %>% group_by(gene_symbol, donor, inside) %>%
    summarise(medianExpression=median(expression), meanExpression = mean(expression))
  medians <- averagesTable %>% select(-meanExpression) %>% spread(inside, medianExpression) %>% dplyr::rename(outside_median=`FALSE`, inside_median=`TRUE`)
  means <- averagesTable %>% select(-medianExpression) %>% spread(inside, meanExpression) %>% dplyr::rename(outside_mean=`FALSE`, inside_mean=`TRUE`)
  
  wilcoxResultsTable <- inner_join(means, wilcoxResultsTable)
  wilcoxResultsTable <- inner_join(medians, wilcoxResultsTable)
  
  wilcoxResultsTable %<>% mutate(fold_change = inside_mean/outside_mean) %>% mutate(log2_fold_change = log2(fold_change))
  
  if(length(donor_group) == 1 ) { # single brain
    (geneStatistics <- wilcoxResultsTable %>% group_by(gene_symbol) %>%
       summarise(medianDirection = median(direction),
                 directionSum = median(direction),
                 metaP.pos = pvalue.pos,
                 metaP.neg = pvalue.neg,
                 average_log2_fold_change = log2_fold_change))
  } else {
    (geneStatistics <- wilcoxResultsTable %>% group_by(gene_symbol) %>%
       summarise(medianDirection = median(direction),
                 directionSum = sum(direction),
                 metaP.pos = sumlog(pvalue.pos)$p,
                 metaP.neg = sumlog(pvalue.neg)$p,
                 average_log2_fold_change = mean(log2_fold_change)))
  }
  
  geneStatistics %<>% arrange(metaP.pos)
  
  geneStatistics %<>% filter(!is.na(average_log2_fold_change)) #used for RNASeq data which has genes with no variance in these samples
  
  geneStatistics %<>% mutate(metaP.pos.adj = p.adjust(metaP.pos, method="fdr"))
  geneStatistics %<>% mutate(metaP.neg.adj = p.adjust(metaP.neg, method="fdr"))
  
  #write out wilcoxResultsTable
  wilcoxResultsTable %<>% arrange(pvalue.pos)
  result_folder <- here("results", exp_matrix, target_column, paste("donors", paste(donor_group, collapse = "_"), sep="_"))
  
  dir.create(result_folder, showWarnings = F, recursive = T)
  
  wilcoxResultsTable %>% write_csv(file.path(result_folder, "wilcoxResultsTable.csv"))
  
  #sort it by whichever p-value is lower and add rank column
  geneStatistics %<>% arrange(metaP.pos)
  geneStatistics %<>% mutate(pos.rank = rank(metaP.pos), neg.rank = rank(-1*metaP.neg))
  geneStatistics %<>% mutate(rank = if_else(metaP.pos < metaP.neg, pos.rank, neg.rank))
  geneStatistics %<>% arrange(rank) %>% dplyr::select(-pos.rank, -neg.rank)
  
  #add in id's
  symbol_id_map <- read_csv(here("data/R_processed_expression_data/","symbol_to_IDs_table.csv"))
  geneStatistics <- left_join(geneStatistics, symbol_id_map) %>% select(gene_symbol, entrez_id, everything())
  
  #write out statistics
  geneStatistics %>% write_csv(file.path(result_folder, "geneStatistics.csv"))
}

