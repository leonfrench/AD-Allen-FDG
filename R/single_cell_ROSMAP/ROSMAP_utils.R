library(here)
library(readr)
library(readxl)
library(dplyr)
library(magrittr)
library(purrr)
library(tidyr)
library(ggplot2)
library(plotROC)
source(here("R", "AUCFunction.R"))
# define target_sheet 
# target_sheet in c("Ex", "Mic", "In", "Ast", "Oli", "Opc")

sheet_path <- here("data","Mathys_Jose-Davila-Velderrain_et_al.","41586_2019_1195_MOESM4_ESM.xlsx")

# helper function to simplify cleaning
extract_3_blocks_of_expression <- function(sheet_path, target_sheet_celltype) {
  whole_sheet <- read_xlsx(sheet_path, sheet = target_sheet_celltype)
  block_names <-
    colnames(whole_sheet)[!grepl("[.][.][.]", colnames(whole_sheet))]
  
  #some code duplication
  first_block <-
    read_xlsx(sheet_path, sheet = target_sheet_celltype, range = "A2:I20000")
  first_block %<>% mutate(contrast_name = block_names[1]) %>% rename(gene_symbol = `...1`) %>% filter(!is.na(gene_symbol))
  first_block %<>% mutate(MixedModel.z = as.character(MixedModel.z),
                          MixedModel.p = as.character(MixedModel.p))
  
  second_block <-
    read_xlsx(sheet_path, sheet = target_sheet_celltype, range = "L2:T20000")
  second_block %<>% mutate(contrast_name = block_names[2]) %>% rename(gene_symbol = `...1`) %>% filter(!is.na(gene_symbol))
  second_block %<>% mutate(MixedModel.z = as.character(MixedModel.z),
                           MixedModel.p = as.character(MixedModel.p))
  
  third_block <-
    read_xlsx(sheet_path, sheet = target_sheet_celltype, range = "W2:AE20000")
  third_block %<>% mutate(contrast_name = block_names[3]) %>% rename(gene_symbol = `...1`) %>% filter(!is.na(gene_symbol))
  third_block %<>% mutate(MixedModel.z = as.character(MixedModel.z),
                          MixedModel.p = as.character(MixedModel.p))
  
  return(list(first_block, second_block, third_block))
}

extract_celltype_expression <- function(sheet_path, target_sheet_celltype) {
  whole_sheet <- read_xlsx(sheet_path, sheet = target_sheet_celltype)
  block_names <-
    colnames(whole_sheet)[!grepl("[.][.][.]", colnames(whole_sheet))]
  
  expression_blocks <- extract_3_blocks_of_expression(sheet_path = sheet_path, target_sheet_celltype = target_sheet_celltype)
  first_block <- expression_blocks[[1]]
  second_block <- expression_blocks[[2]]
  third_block <- expression_blocks[[3]]
  
  # avg the 2 no.pathology.mean values and 2 early.pathology.mean values
  mean_values <-
    inner_join(
      first_block %>% select(gene_symbol, no.pathology.mean),
      second_block %>% select(gene_symbol, no.pathology.mean, early.pathology.mean),
      by = 'gene_symbol'
    ) %>%
    mutate(no.pathology.mean = pmap_dbl(select(., no.pathology.mean.x, no.pathology.mean.y),
                                        function(...)
                                          mean(c(...)))) %>%
    select(-c(no.pathology.mean.x, no.pathology.mean.y)) %>%
    inner_join(third_block %>% select(gene_symbol, early.pathology.mean, late.pathology.mean),
               by = 'gene_symbol') %>%
    mutate(early.pathology.mean = pmap_dbl(select(., early.pathology.mean.x, early.pathology.mean.y),
                                           function(...)
                                             mean(c(...)))) %>%
    select(-c(early.pathology.mean.x, early.pathology.mean.y))
  
  mean_values# %>% select(gene_symbol, no.pathology.mean, early.pathology.mean, late.pathology.mean)
}

transform_celltype_measures_to_tidy<- function(sheet_path, target_sheet_celltype) { #, target_genes
  expression_blocks <- extract_3_blocks_of_expression(sheet_path = sheet_path, target_sheet_celltype = target_sheet_celltype)
  first_block <- expression_blocks[[1]]
  second_block <- expression_blocks[[2]]
  third_block <- expression_blocks[[3]]
  universe <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% .$gene_symbol
  
  one_cell_type <- bind_rows(first_block, second_block, third_block)
  tidy <-
    one_cell_type %>%
    filter(gene_symbol %in% universe) %>%
    group_by(contrast_name) %>%
    select(gene_symbol, IndModel.FC, contrast_name, IndModel.adj.pvals) %>%
    mutate(IndModel.adj.pvals.rank = rank(IndModel.adj.pvals))
  
  tidy %<>% group_by(contrast_name) %>% mutate(directioned_rank = (IndModel.adj.pvals.rank - max(IndModel.adj.pvals.rank)) * sign(IndModel.FC))
  tidy %<>% mutate(directioned_rank = rank(directioned_rank))
  #tidy %<>% mutate(is_target_gene = gene_symbol %in% target_genes)
  #tidy %>% group_by(contrast_name, is_target_gene) %>% summarize(n = dplyr::n(),
  #                                                               directioned_rank = median(directioned_rank, na.rm = T))
  tidy
}

calculate_wilcox <- function(tidy_celltype_data, target_genes) {
  tidy_celltype_data %<>% mutate(is_target_gene = gene_symbol %in% target_genes)  
  wilcox_results <- tidy_celltype_data %>% group_by(contrast_name) %>% summarize(
      auc = auroc_analytic(rank(-directioned_rank), as.numeric(is_target_gene)),
      pValue = wilcox.test(directioned_rank ~ is_target_gene , conf.int = F)$p.value
    )
  wilcox_results
}

get_auc_raster_plots <- function(cell_type_data, contrast_name_target, x_left_label, x_right_label, gene_annotation) {
  #convert to first letter upper case
  x_label <- paste(toupper(substr(contrast_name_target, 1, 1)), substr(contrast_name_target, 2, nchar(contrast_name_target)), sep="")
  microglia_contrast_specific <- cell_type_data %>% filter(contrast_name == contrast_name_target)
  microglia_contrast_specific %<>% select(gene_symbol, directioned_rank)
  
  microglia_contrast_specific <- left_join(microglia_contrast_specific, gene_annotation) %>% 
    spread(key=gene_group, value=gene_group) %>% 
    gather(gene_group, present, -gene_symbol, -directioned_rank) %>% 
    filter(gene_group != "<NA>")
  
  microglia_contrast_specific %<>% mutate(present = if_else(is.na(present), 0, 1)) 
  
  #print AUC scores
  microglia_contrast_specific %>% group_by(gene_group) %>% mutate(rank = rank(directioned_rank)) %>% summarize(auc = auroc_analytic(rank, present))
  
  (AUCPlot <- ggplot(microglia_contrast_specific, aes(d = present, m = directioned_rank, color=gene_group)) + ylab("") + 
      style_roc() + coord_cartesian(expand=F) +
      geom_roc(n.cuts=0, linetype = "solid") + 
      geom_abline(slope = 1, intercept = 0, colour = "grey90", size = 0.2) +
      labs(color='Gene Group')  + 
      theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) +
      #theme(legend.key.width = unit(2, "line"), legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) +
      theme(legend.position = "none") 
  )
  
  microglia_contrast_specific %<>% group_by(gene_group) %>% mutate(rank = rank(-directioned_rank)) #for the raster
  (rasterPlot <- ggplot(microglia_contrast_specific, aes(x = rank, y = present)) +
      geom_blank() + 
      geom_vline(data = filter(microglia_contrast_specific, present == 1), aes(xintercept=rank)) + #,color="black") + #, size=0.07) + 
      theme_bw()+coord_cartesian(expand=F) +
      scale_x_continuous(name = paste0( x_label," (",length(unique(microglia_contrast_specific$gene_symbol))," genes)"), breaks= c(min(microglia_contrast_specific$rank)+1000, max(microglia_contrast_specific$rank)-1300), labels = c(x_left_label, x_right_label)) + 
      theme(axis.title.y = element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank()) 
  )
  list(AUCPlot = AUCPlot, rasterPlot = rasterPlot)
}

