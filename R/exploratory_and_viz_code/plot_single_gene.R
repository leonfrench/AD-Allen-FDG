library(ggsignif)
library(here)
library(readr)
library(metap)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)

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

mni_map <- c("Allen")
pet_map <- c("Schwartz")


donor_group <- c("10021")
donor_group <- unique(all_sample_annotations$donor)

print(paste("Donor group: ", donor_group))


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

target_gene_symbol <- "RPL34"
#target_gene_symbol <- "VASP"
#target_gene_symbol <- "CNTN6"


donors_expression_matrix %<>% mutate(Region = if_else(inside, "Hypometabolism", "Normal")) %>% mutate(Region = factor(Region, levels = c("Normal", "Hypometabolism")))

ggplot(donors_expression_matrix %>% filter(gene_symbol == target_gene_symbol), aes(x = Region, y = expression)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .3) +
  facet_wrap(.~donor, scales="free") + xlab("Alzheimer's FDG-PET Classification") +
  geom_signif(comparisons = list(c("Hypometabolism", "Normal")), test = "wilcox.test") +
  theme_bw()
