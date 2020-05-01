library(here)
library(readr)
library(dplyr)
library(magrittr)
source(here('R', 'AUCFunction.R'))
source(here("R", "GeneSetBuilders.R"))
getwd()

full_expression_matrix <-
  read_csv(
    here(
      "/data/R_processed_expression_data/microarray_expression_probe_averaged_qcNames.csv"
    )
  )
gene_universe <- unique(full_expression_matrix$gene_symbol)

prefixed_lists <- c("GO")
gene_sets <- loadFileSets(prefixed_lists)
srp_genes <-
  tbl_df(gene_sets$MODULES2GENES$SRPdependentTranslationalProteinTargetingMembrane)

# data downloaded from https://github.com/maggiecrow/DEprior
deprior <-
  read_delim(here('data/DEprior-master/', 'DE_Prior.txt'), delim = '\t')
deprior %<>% filter(Gene_Name %in% gene_universe)
# rerank gene_order because universe of genes was decreased
deprior %<>% mutate(Gene_Order = rank(deprior$Gene_Order))
deprior %<>% mutate(is_srp = Gene_Name %in% srp_genes$value)
# 62 of 87 srp genes are matched in deprior
length(srp_genes$value)
sum(deprior$is_srp)

# plot distribution of gene_order of SRP genes
deprior %>% 
  filter(is_srp) %>%
  select(Gene_Order) %>%
  ggplot(aes(x = .$Gene_Order)) +
  geom_histogram()

# reverse gene_order from de_prior to get the correct direction
auc <- auroc_analytic(deprior$Gene_Order, rev(deprior$is_srp))

#Alt methods
library(pROC)
head(deprior)
# note the direction
wilcox.test(Gene_Order ~ is_srp, deprior)
auc(deprior$is_srp, deprior$Gene_Order)
auc(deprior$is_srp, rev(deprior$Gene_Order))
