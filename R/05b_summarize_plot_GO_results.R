library(ggplot2)
library(stringr)
library(here)
library(readr)
library(tidyr)
library(magrittr)

#load gene sets from other gene sets folder 
source(here("R", "GeneSetBuilders.R"))

base_folder <- "/results/microarray/Allen_Schwartz_NiftiValue/figures"
dir.create(here(base_folder), recursive = TRUE)

detach("package:dplyr", unload = TRUE)
library(dplyr)

prefixed_lists <- c("GO") 
gene_sets <- loadFileSets(prefixed_lists)


use_rnaseq <- FALSE
if (use_rnaseq) {
  expression_dir <- 'rnaseq'
} else {
  expression_dir <- 'microarray'
}

#store all donor data in one table
all_adj_results <- as_tibble()


for(result_single_folder in list.files(here("results", expression_dir), pattern = "_NiftiValue", full.names = F)) {
  print(result_single_folder)
  for(run_single_folder in list.files(here("results", expression_dir, result_single_folder), pattern = "donor", full.names = F)) {
    print(run_single_folder)
    
    run_gene_stats <- read_csv(here("results", expression_dir, result_single_folder, run_single_folder, "geneStatistics.csv"))
    run_gene_stats %<>% mutate(donor_set = run_single_folder, result_set = result_single_folder)
    run_gene_stats %<>% mutate(adj = (metaP.pos.adj <= 0.025) | (metaP.neg.adj <= 0.025))
    all_adj_results <- bind_rows(all_adj_results, run_gene_stats)
  } 
}

all_adj_results %<>% mutate(donor_count = str_count(donor_set, "_"))
all_adj_results_SingleD <- all_adj_results %>% filter(donor_count==1)

##Setting up dataframe to plot data
#add proper donar labels to the x axis
donor_info <- read_csv(here("/data/Donor_information.csv")) %>% select(donor_set = folder_name, id) %>% mutate(donor_set = as.character(donor_set))
all_adj_results_SingleD %<>% mutate(donor_set = gsub("donors_","", donor_set))
all_adj_results_DonorLabel <- inner_join(all_adj_results_SingleD, donor_info) %>% mutate(DonorLabel = paste(donor_set, id, sep="/"))


#add avgerage line: where direction of expression changes for all the genes  
mid_ranks <- all_adj_results_SingleD %>% group_by(donor_set) %>% arrange(abs(metaP.pos - metaP.neg)) %>% summarize(hline = dplyr::first(rank))
all_adj_results_DonorLabel_hline <- inner_join(all_adj_results_DonorLabel, mid_ranks)

sum(all_adj_results_SingleD$adj)

create_dot_plot <- function(gene_list, y_label) {
  #extract gene statistics for each GO group of interest by merging all data with only genes in the GO group to be plotted 
  gene_table <- inner_join(all_adj_results_DonorLabel_hline, gene_list, by=c("gene_symbol" = "value"))
  
  # add a horizontal line to summarise avg expression of SRP genes for each donor
  expression_summary <- gene_table %>% group_by(DonorLabel) %>% summarise(median_rank = median(rank))
  gene_table %<>% inner_join(expression_summary, by='DonorLabel')
  
  dot_plot <- ggplot(gene_table, aes(x=factor(DonorLabel), y=rank, fill=adj)) +
    geom_errorbar(aes(ymax=hline, ymin=hline), colour="grey", linetype=2) + # adds horizontal error bar for flip point from over-expressed to under-expressed
    geom_errorbar(aes(ymax=median_rank, ymin=median_rank), colour="black", linetype=1) + # adds horizontal error bar for mean expression rank of SRP genes
    scale_fill_manual(values=c("black", "red")) +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.8, binwidth = 200, stackgroups = F, binpositions = "bygroup", method = "histodot") +
    scale_x_discrete(name = "Donors")+ 
    scale_y_reverse( lim=c(15300,0), labels = c("Underexpression\n15000\n","10000","5000","\n0\nOverexpression"),
                     breaks = c(15000,10000,5000,0),
                     name = y_label, expand = c(0.01,0))+
    theme_bw(base_size = 14)+
    guides(fill = FALSE) +
    theme(panel.border = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"))
  
  dot_plot
}

#load GO SRP gene list to plot for each donor
srp_genes <- tbl_df(gene_sets$MODULES2GENES$SRPdependentTranslationalProteinTargetingMembrane)

srp_y_label <- "SRP-dependent cotranslational genes ranked by differential expression \n in hypometabolic regions associated with Alzheimer's disease"
srp_dot_plot <- create_dot_plot(srp_genes, srp_y_label)
ggsave(here(base_folder,"SRP-dependent_cotranslational_genes_ranked_differential_expression.pdf"), srp_dot_plot, dpi=300, width = 12*.9, height =  9.15*.85)

#supplmentary figure for mitochondrial ribosome GO term genes:
#load MR gene list to plot the genes  
mr_genes <- tbl_df(gene_sets$MODULES2GENES$MitochondrialRibosome)

mr_y_label <- "Mitochondrial ribosome genes ranked by differential expression \n in hypometabolic regions associated with Alzheimer's disease"
mr_dot_plot <- create_dot_plot(mr_genes, mr_y_label)
ggsave(here(base_folder,"Mitochondrial_Ribosome_genes_ranked_differential_expressiontest.pdf"), mr_dot_plot, dpi=300, width = 12*.9, height =  9.15*.85)


#load Darm microglia markers
prefixed_lists <- c("Darmanis")
gene_sets <- loadFileSets(prefixed_lists)

target_genes <- tbl_df(gene_sets$MODULES2GENES$Microglia)
microglia_y_label <- "Microglia marker genes ranked by differential expression \n in hypometabolic regions associated with Alzheimer's disease"
darmanis_dot_plot <- create_dot_plot(target_genes, microglia_y_label)
ggsave(here(base_folder,"Microglia_genes_ranked_differential_expressiontest.pdf"), darmanis_dot_plot, dpi=300, width = 12*.9, height =  9.15*.85)

target_genes <- tbl_df(gene_sets$MODULES2GENES$Neuron)
neuron_y_label <- "Neuron marker genes ranked by differential expression \n in hypometabolic regions associated with Alzheimer's disease"
darmanis_dot_plot <- create_dot_plot(target_genes, neuron_y_label)
ggsave(here(base_folder,"Neuron_genes_ranked_differential_expressiontest.pdf"), darmanis_dot_plot, dpi=300, width = 12*.9, height =  9.15*.85)
