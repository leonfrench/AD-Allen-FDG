library(here)
library(readr)
library(readxl)
library(homologene)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

universe <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% .$gene_symbol
target_genes <- read_csv(here("/data/gene_lists/GO.SRPdependentTranslationalProteinTargetingMembrane.txt"), col_names = F) %>% .$X1
universe_mouse <- unique(human2mouse(universe)$mouseGene)
target_genes_mouse <- unique(human2mouse(target_genes)$mouseGene)

mouse_data <- read_xlsx(here("data", "Keren-Shaul_et_al.","1-s2.0-S0092867417305780-mmc1.xlsx"))
mouse_data %<>% rename(gene_symbol = NAME)
names(mouse_data) %<>% str_split(', ') %>% map(1) 
mouse_data %<>% filter(gene_symbol %in% universe_mouse)
mouse_data %<>% mutate(is_target_gene = gene_symbol %in% target_genes_mouse)


# mark genes with a stepwise increase
mouse_data %<>% mutate(step_wise_increase = (Microglia1 <  Microglia2) & (Microglia2 <  Microglia3))
mouse_data %>% group_by(is_target_gene, step_wise_increase) %>% count()

forMousePlot <- 
  mouse_data %>% 
  filter(is_target_gene == TRUE) %>% 
  select(gene_symbol, Microglia1, Microglia2, Microglia3) %>%
  rename(Normal = Microglia1, `Group II DAM` = Microglia2, `Group III DAM` = Microglia3)

forMousePlot %<>% 
  gather(key='celltype', value = 'expression', -gene_symbol) %>% 
  group_by(gene_symbol) %>% 
  mutate(zscore=scale(expression)) %>%
  select(-expression) %>% 
  spread(key='celltype', value='zscore') %>% 
  #mutate(diff3_1 = Microglia3 - Microglia1) %>% 
  #arrange(diff3_1) %>% 
  gather(key='celltype', value='z_expression', -gene_symbol) #-c(gene_symbol, diff3_1))

forMousePlot$celltype <- factor(forMousePlot$celltype, levels = c("Normal", "Group II DAM", "Group III DAM"))
forMousePlot %<>% rename(mouse_gene_symbol = gene_symbol, sample_origin = celltype)

###############################################################################################################################################
# Processing human data
###############################################################################################################################################
source(here('R', 'single_cell_ROSMAP', 'ROSMAP_utils.R'))
excel_path <- here("data","Mathys_Jose-Davila-Velderrain_et_al./41586_2019_1195_MOESM4_ESM.xlsx")
microglia_expression <- extract_celltype_expression(sheet_path = excel_path, target_sheet_celltype = 'Mic')
microglia_expression %<>%
  rename(`No Pathology` = no.pathology.mean, 
         `Early Pathology` = early.pathology.mean, 
         `Late Pathology` = late.pathology.mean)

forHumanPlot <- 
  microglia_expression %>% 
  filter(gene_symbol %in% target_genes) %>% 
  gather(key='pathology', value = 'expression', -gene_symbol) %>% 
  group_by(gene_symbol) %>% 
  mutate(z_expression = scale(expression)) %>% 
  select(-expression)

forHumanPlot$pathology <- factor(forHumanPlot$pathology, levels = c("No Pathology", "Early Pathology","Late Pathology"))
forHumanPlot %<>% rename(human_gene_symbol = gene_symbol, sample_origin = pathology)

# add humanGenes to mouse DF
# in order to bind rows, need to have the same columns in both DFs
forMousePlot %<>% mutate(human_gene_symbol = mouse2human(mouse_gene_symbol)$humanGene)
forMousePlot %>% filter_all(any_vars(is.na(.)))
# add mouseGene col to humanDF
mouse_analog_genes <- human2mouse(forHumanPlot$human_gene_symbol) %>% select(mouseGene, humanGene)
forHumanPlot %<>% left_join(mouse_analog_genes, by=c('human_gene_symbol' = 'humanGene')) %>% rename(mouse_gene_symbol = mouseGene)
forHumanPlot %>% filter_all(any_vars(is.na(.)))

# check dimensions match
dim(forMousePlot)
dim(forHumanPlot)

forHumanPlot %<>% mutate(species = 'human')
forMousePlot %<>% mutate(species = 'mouse')

forPlot <- bind_rows(forMousePlot, forHumanPlot) 
forPlot$sample_origin <- factor(forPlot$sample_origin, levels = c("Normal", "Group II DAM", "Group III DAM", "No Pathology", "Early Pathology","Late Pathology"))
forPlot$species <- factor(forPlot$species, levels = c('mouse', 'human'))
forPlot %<>% filter(!is.na(mouse_gene_symbol)) #%>% filter(!is.na(diff3_1))


###############################################################################################################################################
# working on reordering the heatmap
# drop NAs, hclust
no_mouse_genes <- forPlot %>% filter_all(any_vars(is.na(.)))

forClust <- 
  forPlot %>% 
  select(-species) %>%
  spread(key = sample_origin, value=z_expression) %>%
  drop_na() 

# helpful to use dataframe with rownames for hclust
forDistanceMatrix <- forClust %>% as.data.frame()
rownames(forDistanceMatrix) <- forClust$mouse_gene_symbol
forDistanceMatrix %<>% select(-mouse_gene_symbol, -human_gene_symbol)

#dist_matrix <- 
#  forClust %>% 
#  select(Normal, `Group II DAM`, `Group III DAM`, `No Pathology`, `Early Pathology`, `Late Pathology`) %>% #mouse_gene_symbol
#  dist()

hc <- hclust(dist(forDistanceMatrix))
plot(hc)

clustered_gene_order <- as_tibble(hc$labels[hc$order])
clustered_gene_order$rank_order <- 1:nrow(clustered_gene_order)


forClustHeatmap <-
  forClust %>%
  left_join(clustered_gene_order, by = c('mouse_gene_symbol' = 'value')) %>%
  #arrange(rank_order) %>%
  gather(key = sample_origin, value = z_expression, -mouse_gene_symbol, -human_gene_symbol, -rank_order) %>%
  mutate(species = case_when(
    sample_origin %in% c("Normal", "Group II DAM", "Group III DAM") ~ 'mouse',
    sample_origin %in% c("No Pathology", "Early Pathology", "Late Pathology") ~ 'human'
  ))

forClustHeatmap$sample_origin <- factor(forClustHeatmap$sample_origin, levels = c("Normal", "Group II DAM", "Group III DAM", "No Pathology", "Early Pathology","Late Pathology"))
forClustHeatmap$species <- factor(forClustHeatmap$species, levels = c('mouse', 'human'))

# some genes require duplicating to fix the match up of mouse-human homologs
forClustHeatmap %>% filter(mouse_gene_symbol == 'Gm10071')
forClustHeatmap %>% filter(human_gene_symbol == 'RPL13')
forClustHeatmap %>% filter(mouse_gene_symbol == 'Gm5428')
forClustHeatmap %>% filter(human_gene_symbol == 'RPL6')
forClustHeatmap %>% filter(mouse_gene_symbol == 'Srp54a') %>% mutate(human_gene_symbol = 'SRP54 ') 
#forHeat %>% filter(mouse_gene_symbol == 'Srp54b')

fixed_Gm5428 <- forClustHeatmap %>% filter(mouse_gene_symbol == 'Gm5428') %>% mutate(human_gene_symbol = 'RPL6 ') 
fixed_Gm10071 <- forClustHeatmap %>% filter(mouse_gene_symbol == 'Gm10071') %>% mutate(human_gene_symbol = 'RPL13 ') 
fixed_Srp54a <- forClustHeatmap %>% filter(mouse_gene_symbol == 'Srp54a') %>% mutate(human_gene_symbol = 'SRP54 ') 

genes_to_duplicate <- c('Gm10071', 'Gm5428', 'Srp54a')
# remove those that require duplication then bind in the fixed versions
forClustHeatmap %<>% filter(!mouse_gene_symbol %in% genes_to_duplicate)
forClustHeatmap %<>%  bind_rows(fixed_Gm10071, fixed_Gm5428, fixed_Srp54a)

# create 2 versions of the same plot, with mouse gene labels on left y-axis, and human labels on right y-axis
clustHeatMouse <- ggplot(forClustHeatmap, aes(x=sample_origin, y=reorder(mouse_gene_symbol, rank_order))) +
  facet_wrap(~species, scales='free_x') +
  geom_tile(aes(x=sample_origin, fill = z_expression), colour='black') +
  labs(x=NULL, fill = "Expression\n(z-score)") +
  scale_x_discrete(expand = c(0, 0)) + 
  theme(axis.text.y  = element_text(size=8)) +
  ylab('') +
  scale_fill_distiller(palette = 'RdBu') +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5), 
        plot.margin = grid::unit(c(t=0,r=0,b=2, l=0), "mm"),
        strip.text = element_blank())

clustHeatMouse
ggsave(filename = here('results', 'mouse', "DAM_dataset", 'clustered_heatmap_mouse.png'), dpi=300, width=4, height=12)


clustHeatHuman <- ggplot(forClustHeatmap, aes(x=sample_origin, y=reorder(human_gene_symbol, rank_order))) +
  facet_wrap(~species, scales='free_x') +#, space="free_x") + 
  geom_tile(aes(x=sample_origin, fill = z_expression), colour='black') +
  labs(x=NULL, fill = "Expression\n(z-score)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(position = 'right') +
  theme(axis.text.y  = element_text(size=8)) +
  ylab('') +
  scale_fill_distiller(palette = 'RdBu') +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5), 
        plot.margin = grid::unit(c(t=0,r=0,b=2, l=0), "mm"),
        strip.text = element_blank())

clustHeatHuman
ggsave(filename = here('results', 'mouse', "DAM_dataset", 'clustered_heatmap_human.png'), dpi=300, width=4, height=12)

