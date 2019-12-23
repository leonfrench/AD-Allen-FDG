library(huxtable)
library(here)
library(readr)
library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(homologene)
library(pROC)
source(here('R', 'AUCFunction.R'))


apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

convert_genes <- function(input_genes) {
  mouse_genes <- human2mouse(input_genes)
  return(unique(mouse_genes$mouseGene))
}


get_enrichment <- function(df, gene_list) {
  # df is matrix of ranks, with gene_symbol + celltypes as columns
  
  #forIndices <- as_tibble(rownames(linnarssonMatrixHumanReachable))
  forIndices <- as_tibble(df$gene_symbol)
  names(forIndices) <- 'gene_symbol'
  forIndices %<>% mutate(isTargetGene = gene_symbol %in% gene_list)
  
  targetIndices <- forIndices$isTargetGene
  print(paste0(
    sum(forIndices$isTargetGene),
    ' of ',
    length(gene_list),
    ' genes were found in the human reachable linnarson matrix'
  ))

  df %<>% select(-gene_symbol)
  AUROC <- map_df(df, auroc_analytic, targetIndices)
  MWU <- map_df(df, apply_MWU, targetIndices)
  
  enrichment <- gather(AUROC, key = 'celltype', value = 'AUC')
  pvals <- gather(MWU, key = 'celltype', value = 'pValue')
  enrichment %<>% inner_join(pvals, by = 'celltype') %>% mutate(rank = rank(desc(AUC))) %>% arrange(desc(AUC))
  #merge in celltype descriptions
  descriptions <-
    read_excel(here(
      'data',
      'Zeisel.TableS3.1-s2.0-S009286741830789X-mmc3.xlsx'
    )) %>% select(`Cluster name`, `Description`)
  enrichment %<>% left_join(descriptions, by = c('celltype' = 'Cluster name'))
  
  #cholinergic_celltypes <- c('DECHO1','DECHO2','ENT4','ENT5','ENT6','ENT7','ENT8','ENT9','HBCHO1','HBCHO2','MBCHO1','SYCHO1','SYCHO2','TECHO')
  cholinergic_celltypes <- enrichment %>% filter(grepl("holine", Description)) %>% .$celltype
  
  enrichment %<>% mutate(is_cholinergic = celltype %in% cholinergic_celltypes)
  enrichment %<>% mutate(adjusted_P = signif(p.adjust(pValue, method = 'fdr'), digits = 3))
}
# load full dataset:
# --> linnarssonMatrixMouse // linnarssonMatrixHumanReachable

load(here('data', 'Zeisel_etal_2018', 'l5_all.agg.tab.processed.RData'), verbose = TRUE)

# we have a matrix of gene ranks by celltype
head(linnarssonMatrixHumanReachable)

# process matrix to be operating in same gene_universe
df <- as_tibble(linnarssonMatrixHumanReachable)
celltypes <- names(df) # store celltype names in order to select them for re-ranking after filtering for gene_universe
df$gene_symbol <- rownames(linnarssonMatrixHumanReachable)
gene_universe <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% .$gene_symbol
mouse_gene_universe <- convert_genes(gene_universe)
# requires re-ranking the genes in the matrix after filtering some genes
df %<>% filter(gene_symbol %in% mouse_gene_universe) %>% modify_at(celltypes, rank)

# define gene list to test and convert symbols to mouse
# list loses 4 gene symbols in conversion from human to mouse
gene_list <- read_csv(here('data', 'gene_lists', 'GO.SRPdependentTranslationalProteinTargetingMembrane.txt'), col_names = 'gene_symbol')


mouse_gene_list <- convert_genes(gene_list$gene_symbol)
dim(gene_list)
length(mouse_gene_list)

srp_enrichment <- get_enrichment(df, mouse_gene_list)
chat_enrichment <- get_enrichment(df, c('Chat'))

dir.create(here("results/mouse/"), showWarnings = F, recursive = T)
supplement_table <- srp_enrichment %>% 
  select(Cluster_ID=celltype, Name=Description, AUC, pValue, adjusted_P, Rank=rank)
supplement_table %>% write_csv(here("results", "mouse", "Supplement table - All Zeisel clusters.csv"))

final_table <- srp_enrichment %>% 
  filter(is_cholinergic == TRUE) %>% 
  select(Cluster_ID=celltype, Name=Description, AUROC=AUC, pValue, adjusted_P, Rank=rank)

ht <- hux(final_table, add_colnames = TRUE, autoformat = TRUE) %>% 
  set_right_padding(1)           %>%
  set_left_padding(1)            %>% 
  set_bold(1, 1:6, TRUE)          %>% 
  set_bottom_border(1, 1:6, 1)    %>%
  set_right_border(everywhere, 1:5, 0.4) %>% 
  set_align(1:14, 3:6, 'right')      %>%
  #set_number_format(everywhere, 2)         %>% 
  set_caption('Cholinergic enrichment')

theme_plain(ht)
quick_html(ht, file = here("results", "mouse", "Table - Zeisel top clusters.html"), open=F)


# define different subsets of the data to facilitate plotting
merged_enrichment <-
  inner_join(
    srp_enrichment %>% select(celltype, AUC_SRP = AUC, is_cholinergic, Description),
    chat_enrichment %>% select(celltype, AUC_CHAT = AUC),
    by = 'celltype'
  )

sub1 <- merged_enrichment %>% filter(AUC_SRP > 0.75)
sub2 <- merged_enrichment %>% filter((0.75 > AUC_SRP) & (AUC_SRP > 0.625))
sub3 <- merged_enrichment %>% filter(AUC_SRP < 0.625)

srp_chat_enrichment <-
  ggplot(merged_enrichment, aes(
    y = AUC_SRP,
    x = AUC_CHAT,
    label = paste(Description, celltype)
  )) +
  geom_point(color = if_else(merged_enrichment$is_cholinergic, 'red', 'black')) +
  geom_text_repel(
    data = sub1 %>% filter(is_cholinergic),
    xlim = 1.05,
    hjust = 1,
    force = 9,
    #box.padding = 0.15,
    direction = 'both',
    size = 3,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub2 %>% filter(is_cholinergic),
    xlim = 1.05,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.1,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub3 %>% filter(is_cholinergic),
    xlim = 1.05,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.2,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  scale_x_continuous(breaks = c(0.25, 0.50, 0.75, 1.00),
                     limits = c(NA, 1.5)) +
  theme_bw()

ggsave(filename = here('results','mouse', 'zeisel-srp_chat_enrichment.png'),
       width = 12, height = 9,
       dpi=300)



### AUC of AUCs
# 
head(srp_enrichment)
wilcox.test(data=srp_enrichment, rank ~ is_cholinergic)
t.test(data=srp_enrichment, AUC ~ is_cholinergic) #to see the means
wilcox.test(data=srp_enrichment, AUC ~ is_cholinergic)

srp_enrichment %>% group_by(is_cholinergic) %>% summarize(medianAUC = median(AUC))

roc(srp_enrichment$is_cholinergic, srp_enrichment$AUC)
roc(is_cholinergic ~ AUC , srp_enrichment)


###################################################################
# plots showing distribution of enrichment scores for all celltypes
# highlighting cholinergic celltypes

# raster
ggplot(srp_enrichment, x=rank, y=is_cholinergic) +
  geom_vline(data = filter(srp_enrichment, is_cholinergic == TRUE), aes(xintercept=rank)) +
  xlim(range(srp_enrichment$rank)) +
  coord_cartesian(expand = F) +
  theme(strip.background = element_blank(), strip.placement = "inside") 

#distplot with rug
ggplot(srp_enrichment, aes(x=AUC)) +
  geom_line(stat = 'density') +
  geom_rug(data=filter(srp_enrichment, is_cholinergic==TRUE))
  
#displot with lines for cholinergic cells
ggplot(srp_enrichment, aes(x=AUC)) +
  geom_line(stat = 'density') +
  geom_vline(data=filter(srp_enrichment, is_cholinergic==TRUE), aes(xintercept=AUC))
