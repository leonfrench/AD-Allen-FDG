library(readr)
library(dplyr)
library(corrr)
library(biomaRt)
library(ggplot2)
library(here)

srp_genes <- read_csv(here("data", "gene_lists", "GO.SRP.txt"), col_names = 'gene_symbol')

genestats9861 <- read_csv(here('/results/rnaseq/Allen_Schwartz_NiftiValue/donors_9861/geneStatistics.csv'))
genestats9861 %<>% mutate(donor = 9861)
genestats10021 <- read_csv(here('/results/rnaseq/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv'))
genestats10021 %<>% mutate(donor = 10021)

genestats <- bind_rows(genestats9861, genestats10021)

human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords <- getBM(attributes=c("hgnc_symbol", "start_position","end_position"), filters="hgnc_symbol", values=srp_genes, mart=human)
gene_coords$size=gene_coords$end_position - gene_coords$start_position

dim(genestats)
dim(gene_coords)
srp_size <- left_join(genestats, gene_coords, by=c("gene_symbol" = "hgnc_symbol"))

srp_size %>% 
  filter(!is.na(size)) %>% 
  ggplot(aes(x=rank, y=size, group=donor)) + 
  facet_wrap(. ~ donor) + 
  geom_point()


#srp_size %>% filter(!is.na(size)) %>% dplyr::select(rank, size, donor) %>% group_by(donor) %>% correlate() 
srp_size %>% filter(!is.na(size), donor==9861) %>% dplyr::select(rank, size) %>% correlate()
srp_size %>% filter(!is.na(size), donor==10021) %>% dplyr::select(rank, size) %>% correlate()
