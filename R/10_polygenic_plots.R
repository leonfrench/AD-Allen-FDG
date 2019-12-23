library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(here)

#! This code requires running Allen_expression_for_region_specificty.ipynb in the root file first

enrichment_table <- read_csv(here('data', 'processed_HBA', 'srp_chat_enrichment_table.csv'))

#labels <- c('substantia innominata', 'basal nucleus of meynert', 'nucleus of the diagonal band', 
#            'red nucleus', 'nucleus of the anterior commissure', 'bed  nucleus of stria terminalis',
#            'olfactory tubercle', 'islands of calleja')

labels <- c('substantia innominata', 'basal nucleus of meynert', 'nucleus of the diagonal band', 
            'red nucleus', 'nucleus of the anterior commissure', 'bed  nucleus of stria terminalis',
            'olfactory tubercle', 'islands of calleja')


labels <- c('substantia innominata', 'basal nucleus of meynert', 'nucleus of the diagonal band', 'bed  nucleus of stria terminalis', 'olfactory tubercle','islands of calleja')
#labels <- c('substantia innominata',           'globus pallidus, external segment', 'globus pallidus, internal segment')
#labels <- c('substantia innominata','nucleus accumbens', "head of caudate nucleus")
#add in full names
enrichment_table
donor_info <- read_csv(here("/data/Donor_information.csv")) %>% select(donor = folder_name, id)
enrichment_table <- inner_join(enrichment_table, donor_info) %>% mutate(donor = paste(donor, id, sep="/"))

#setup the donor order
enrichment_table %<>% mutate(donor = factor(donor, levels= unique(c("9861/H0351.2001", "10021/H0351.2002", sort(unique(enrichment_table$donor))))))

ggplot(enrichment_table, aes(y=srp_AUC, x=CHAT, label=brain_structure)) +
  geom_point(data = enrichment_table %>% filter(!str_detect(brain_structure, paste(labels, collapse='|'))), color= 'darkgrey') +
  geom_point(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))), color= 'red') +
  #geom_point(color=if_else(enrichment_table$brain_structure %>% str_detect(paste(labels, collapse='|')), 'red', 'darkgrey')) +
  geom_text_repel(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))),
                  #force = 10,
                  #box.padding = 1.8,
                  size = 2.5,
                  #nudge_y = 0.6,
                  #nudge_x      = 0.15,
                  #nudge_x = -0.5,
                  direction='y',
                  segment.size  = 0.2,
                  hjust = 0,
                  #xlim  = c(10,NA),
                  segment.color = "red") +
  facet_wrap(~ donor) + 
  scale_y_continuous(limits = c(0, 1)) + 
  xlim(min(enrichment_table$CHAT), max(enrichment_table$CHAT)+10) +
  ylab('Regional specificity of ER translocation genes (AUC)') + 
  xlab('CHAT gene expression (log2 intensity)') + 
  theme_bw()

ggsave(filename = here('results', 'microarray', 'HBA-srp_chat_enrichment.png'), dpi=300)
#for PDF use 10x7

# correlations
cor.test(enrichment_table$CHAT, enrichment_table$srp_AUC, method = 's')

enrichment_table %>% 
  group_by(donor) %>% 
  summarise(correlation = cor(srp_AUC, CHAT, method='s'))

enrichment_table %>% 
  group_by(donor) %>% 
  summarise(meanSRP = mean(srp_AUC), meanCHAT = mean(CHAT))

#regions with both high colinergic and high SRP
means <- enrichment_table %>% 
  group_by(brain_structure) %>% 
  summarise(meanSRP = mean(srp_AUC), meanCHAT = mean(CHAT)) %>% arrange(-meanSRP)
means %>% mutate(SRPrank = rank(-meanSRP), CHATrank = rank(-meanCHAT)) %>% mutate(combined = SRPrank+CHATrank) %>% arrange(combined)
