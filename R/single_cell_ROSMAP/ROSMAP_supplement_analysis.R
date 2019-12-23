library(here)
library(readr)
library(cowplot)
source(here("R", "single_cell_ROSMAP",  "ROSMAP_utils.R"))

sheet_path <- here("data","Mathys_Jose-Davila-Velderrain_et_al.","41586_2019_1195_MOESM4_ESM.xlsx")
target_genes <- read_csv(here("/data/gene_lists/GO.SRPdependentTranslationalProteinTargetingMembrane.txt"), col_names = F) %>% .$X1

universe <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% .$gene_symbol

all_wilcox_results <- as_tibble()
all_types_slim <- as_tibble()
for (target_sheet in c("Ex", "Mic", "In", "Ast", "Oli", "Opc")) {
  celltype_slim <- transform_celltype_measures_to_tidy(sheet_path, target_sheet_celltype = target_sheet)
  # perform wilcox tests in specified celltype
  wilcox_results <- calculate_wilcox(tidy_celltype_data=celltype_slim, target_genes=target_genes)
  wilcox_results %<>% mutate(cell_type = target_sheet)
  all_wilcox_results <- bind_rows(all_wilcox_results, wilcox_results)
  
  celltype_slim %<>% mutate(cell_type = target_sheet) %>% mutate(is_target_gene = gene_symbol %in% target_genes) 
  all_types_slim <- bind_rows(all_types_slim, celltype_slim)
}

all_wilcox_results %<>% 
  arrange(pValue) %>% 
  mutate(p.adjust = p.adjust(pValue, method ="fdr")) %>% 
  ungroup() %>% 
  mutate(is_significant = p.adjust < 0.05)

all_wilcox_results %>% group_by(is_significant) %>% count()
all_wilcox_results %>% filter(is_significant) %>% group_by(contrast_name) %>% count()

all_wilcox_results %>% group_by(is_significant, auc > 0.5) %>% count()
mean(all_wilcox_results$auc)
all_wilcox_results %>% filter(cell_type == "Mic")
all_wilcox_results %>% filter(contrast_name == "no-pathology vs early-pathology differential expression")

# median p-value across the group
all_types_slim %>% 
  filter(is_target_gene) %>% 
  filter(cell_type == "Mic") %>%
  group_by(contrast_name) %>% 
  summarize(median_p=median(IndModel.adj.pvals), min_p=min(IndModel.adj.pvals), max_p=max(IndModel.adj.pvals))

###############################################################################################################################################
# Plotting
###############################################################################################################################################
# for AUC plots of microglia
gene_annotation <- as_tibble()
SRP_genes <- target_genes
gene_annotation <- bind_rows(data.frame(gene_symbol = SRP_genes, gene_group = "SRP membrane targeting genes", stringsAsFactors = F), gene_annotation)

microglia <- all_types_slim %>% filter(cell_type == "Mic")
microglia %<>% mutate(directioned_rank = directioned_rank * -1) #flip it for plotting

cell_type_data <- microglia
contrast_name_target <- "no-pathology vs early-pathology differential expression"

microglia_no_vrs_early_plots <- get_auc_raster_plots(microglia, 
                                                     "no-pathology vs early-pathology differential expression", 
                                                     "Overexpressed in early pathology", "Overexpressed in no pathology",
                                                     gene_annotation)

microglia_early_vrs_late_plots <- get_auc_raster_plots(microglia, 
                                                       "early-pathology vs late-pathology differential expression", 
                                                       "Overexpressed in late pathology", "Overexpressed in early pathology",
                                                       gene_annotation)

(four_plots <- plot_grid(microglia_no_vrs_early_plots$AUCPlot, microglia_early_vrs_late_plots$AUCPlot, microglia_no_vrs_early_plots$rasterPlot, microglia_early_vrs_late_plots$rasterPlot, rel_heights=c(1,0.4),align = "v", labels = c("a", "b", "c", "d")))
#save as 11x8 PDF
dir.create(here('results',  "ROSMAP_single_nuc"), showWarnings = FALSE)
ggsave(filename = here('results', "ROSMAP_single_nuc", 'SRP_mouse_microglia_ROCs.png'), dpi=300, width=11, height=8)
ggsave(filename = here('results', "ROSMAP_single_nuc", 'SRP_mouse_microglia_ROCs.pdf'), dpi=300, width=11, height=8, device="pdf")

###############################################################################################################################################
# check for two step increases 
###############################################################################################################################################
# cut it to genes that are in both contrasts
microglia <- all_types_slim %>% filter(cell_type == "Mic")
unique(microglia$contrast_name)
microglia %<>% filter(contrast_name != "no-pathology vs pathology differential expression")
# select genes that are in both contrasts
microglia %<>% filter(!is.na(IndModel.FC))
microglia %<>% inner_join(microglia %>% group_by(gene_symbol) %>% dplyr::count() %>% filter(n == 2))
microglia$gene_symbol %>% unique %>% length 

step_wise_increase_microglia <- microglia %>% group_by(gene_symbol) %>% 
  select(gene_symbol, IndModel.FC, contrast_name) %>% 
  pivot_wider(names_from = contrast_name, values_from = IndModel.FC) %>% 
  mutate(step_wise_increase = `early-pathology vs late-pathology differential expression` > 0 & `no-pathology vs early-pathology differential expression` > 0)

step_wise_increase_microglia %<>% mutate(is_target_gene = gene_symbol %in% target_genes)
step_wise_increase_microglia %>% group_by(is_target_gene) %>% count()
step_wise_increase_microglia %>% group_by(step_wise_increase) %>% count()
step_wise_increase_microglia %>% group_by(is_target_gene, step_wise_increase) %>% count()
step_wise_increase_microglia %>% filter(is_target_gene, step_wise_increase) %>% .$gene_symbol %>% sort

# run through gene ontology 
source(here("R", "GeneSetBuilders.R"))
detach("package:XML", unload = TRUE) 
detach("package:tidyr", unload = TRUE)
detach("package:dplyr", unload = TRUE)
library(dplyr)
library(tidyr)

if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { 
} else {
  geneSetsGO <- loadGOSets(universe)
}

#TMOD for test
universe_microglia <- intersect(universe, step_wise_increase_microglia$gene_symbol)
all_go_runs <- as_tibble(tmodHGtest(fg = step_wise_increase_microglia %>% filter(step_wise_increase) %>% .$gene_symbol, bg=universe_microglia, 
                                    mset = geneSetsGO, qval = 1.01, filter = TRUE))

all_go_runs %>% print()
write_csv(all_go_runs, here('results', "ROSMAP_single_nuc","microglia_increase_step_GO.csv"))

#Specificity check #2
#run tmod to get AUC for all groups
#micro_early_late <- all_types_slim %>% filter(cell_type == "Mic", contrast_name == )
for(contrast_name_target in c("no-pathology vs early-pathology differential expression", "early-pathology vs late-pathology differential expression")) {
  micro_early_late <- all_types_slim %>% filter(cell_type == "Mic", contrast_name == contrast_name_target)
  micro_early_late %>% filter(IndModel.adj.pvals < 0.05, is_target_gene)
  sum(micro_early_late$is_target_gene)
  micro_early_late %<>% arrange(directioned_rank)
  filterGenes <- TRUE
  result <- tbl_df(tmodUtest(c(micro_early_late$gene_symbol), mset=geneSetsGO, qval = 1.01, filter = filterGenes))
  result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
  result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsGO$MODULES2GENES[ID])), collapse = " "))
  result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
    summarize(MainTitle = dplyr::first(Title), ID=dplyr::first(ID), all_IDs=paste(ID, collapse=","), AUC = dplyr::first(AUC), P.Value= dplyr::first(P.Value), aspect= dplyr::first(aspect), otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
  result %<>% ungroup() %>% dplyr::select(-genes)
  
  #double p-value because tmod doesn't 
  result %<>% mutate(P.Value=P.Value*2) 
  
  #adjust again now that identical groups have been merged
  result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value)
  result %<>% mutate(rank = rank(P.Value))
  result %<>% dplyr::select(MainTitle, ID, geneCount = N1, AUC, P.Value, adj.P.Value, everything()) 
  write_csv(all_go_runs, here('results', "ROSMAP_single_nuc", paste0("microglia_AUC_step_GO_", contrast_name_target, ".csv")))
}
