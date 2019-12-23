library(here)
library(readr)
library(readxl)
library(stringr)
library(purrr)
library(homologene)
library(tmod)
source(here("R", "GeneSetBuilders.R"))
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:tidyr", unload = TRUE)
detach("package:dplyr", unload = TRUE)
library(dplyr)
library(tidyr)

universe <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% .$gene_symbol
target_genes <- read_csv(here("/data/gene_lists/GO.SRPdependentTranslationalProteinTargetingMembrane.txt"), col_names = F) %>% .$X1

universe_mouse <- unique(human2mouse(universe)$mouseGene)
target_genes_mouse <- unique(human2mouse(target_genes)$mouseGene)
mouse_data <- read_xlsx(here("data", "Keren-Shaul_et_al.","1-s2.0-S0092867417305780-mmc1.xlsx"))

mouse_data %<>% rename(gene_symbol = NAME) %>% filter(gene_symbol %in% universe_mouse)
universe_mouse <- intersect(mouse_data$gene_symbol, universe_mouse)

mouse_data %<>% mutate(is_target_gene = gene_symbol %in% target_genes_mouse)
mouse_data %>% filter(is_target_gene) %>%  summarise_if(is.numeric, mean, na.rm = TRUE)
mouse_data %>% filter(!is_target_gene) %>%  summarise_if(is.numeric, mean, na.rm = TRUE)

#shorten cell type cluster names
names(mouse_data) %<>% str_split(', ') %>% map(1)

#mark genes with a stepwise increase
mouse_data %<>% mutate(step_wise_increase = (Microglia1 <  Microglia2) & (Microglia2 <  Microglia3))
mouse_data %>% group_by(step_wise_increase) %>% count()
mouse_data %>% group_by(is_target_gene) %>% count()
mouse_data %>% group_by(is_target_gene, step_wise_increase) %>% count()

#test across all GO groups, using human GO annotations
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { 
} else {
  geneSetsGO <- loadGOSets(universe)
}

#convert GO groups to mouse - slow but less code than before
geneSetsGO_mouse <- geneSetsGO
convert_2_mouse_list <- function(human_genes) {
  mouse_genes <- unique(human2mouse(human_genes)$mouseGene)
  intersect(mouse_genes, universe_mouse)
}
start_time <- Sys.time()
geneSetsGO_mouse$MODULES2GENES <- lapply(geneSetsGO_mouse$MODULES2GENES, convert_2_mouse_list)
print(Sys.time() - start_time)
geneSetsGO_mouse <- makeTmod(modules = geneSetsGO_mouse$MODULES, modules2genes = geneSetsGO_mouse$MODULES2GENES)
#run hypergeometric test 
go_results <- as_tibble(tmodHGtest(fg = mouse_data %>% filter(step_wise_increase) %>% .$gene_symbol, bg=universe_mouse, 
                          mset = geneSetsGO_mouse, qval = 1.01, filter = TRUE))

dir.create(here('results', "mouse", "DAM_dataset"), showWarnings = FALSE)
write_csv(go_results, here('results', "mouse",  "DAM_dataset", "GO_results.csv"))

