library(ggplot2)
library(here)
library(readr)
library(stringr)
library(magrittr)
library(dplyr)
library(seqinr)
source(here("R", "GeneSetBuilders.R"))
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)


genecode_fasta <- read.fasta(file = here("data", "protein_gencode", "gencode.v32.pc_translations.fa.gz"))

name_table <- tibble(piped_name = names(genecode_fasta))
name_table %<>% dplyr::mutate(order = 1:n())

#grab gene_symbol from second last token
name_table %<>%
  mutate(
    splits = strsplit(piped_name, "[|]")
  ) %>% 
  rowwise() %>% 
  mutate(
    gene_symbol = splits[length(splits)-1]
  )

gene_universe <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv")) %>% .$gene_symbol

#load the GO sets we've been using throughout
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { 
} else {
  geneSetsGO <- loadGOSets(gene_universe)
}
filterGenes <- TRUE #work within the GO univerise


name_table %<>% mutate(in_universe = gene_symbol %in% gene_universe)
name_table %>% group_by(in_universe) %>% dplyr::summarise(count = length(gene_symbol))

name_table %<>% filter(in_universe)

name_table

#filter fasta file for genes in universe
genecode_fasta <- genecode_fasta[name_table$order]
length(genecode_fasta)

symbol_to_sequence <- data.frame(gene_symbol= name_table$gene_symbol, names=names(genecode_fasta), stringsAsFactors = F) %>% as_tibble()
symbol_to_sequence$protein_sequence = unlist(unname(lapply(genecode_fasta, c2s)))
symbol_to_sequence %<>% mutate(length = nchar(protein_sequence))

gene_stats <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv"))

amino_acid_alphabet <- c("a", "c", "d", "e", "f", "g", "h", "i", "k", "l", "m", "n", "p", "q", "r", "s", "t", "v", "w", "y")

amino_acid_first_target <- "r"
amino_acid_second_target <- "k"
all_pair_results <- as_tibble()
ptm <- Sys.time()
for(amino_acid_first_target in amino_acid_alphabet) {
  for(amino_acid_second_target in amino_acid_alphabet) {
    if (nrow (all_pair_results %>% filter(amino_acid_first == amino_acid_second_target & amino_acid_second == amino_acid_first_target)) == 1) {
      next
    }
    print(paste0(amino_acid_first_target, amino_acid_second_target))
    gene_stats_joined <- NULL
    symbol_to_sequence$target_pair_count <- NULL #clear for safety
    #count the amino acid hits
    symbol_to_sequence %<>% mutate(target_pair_count = str_count(protein_sequence, amino_acid_first_target) + str_count(protein_sequence, amino_acid_second_target))
    symbol_to_sequence %<>% mutate(percent_target_amino_acids = target_pair_count/length)

    #mean average percent across the many entries per gene
    symbol_to_sequence_averaged <- symbol_to_sequence %>% group_by(gene_symbol) %>% summarize(percent_target_amino_acids = mean(percent_target_amino_acids)) %>% arrange(-percent_target_amino_acids)
    
    symbol_to_sequence_averaged %<>% arrange(-percent_target_amino_acids)
    
    mean(symbol_to_sequence_averaged$percent_target_amino_acids)
    
    sortedGenes <- symbol_to_sequence_averaged$gene_symbol
    
    result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1.01, filter = filterGenes))
    
    #filter for positive AUC
    #result %<>% filter(AUC > 0.5)
    
    #merge same GO groups using the genes in the set - these are sets that have the same set of genes
    result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsGO$MODULES2GENES[ID])), collapse = " "))
    result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
      summarize(MainTitle = first(Title), ID=first(ID), all_IDs=paste(ID, collapse=","), 
                AUC = first(AUC), P.Value= first(P.Value), 
                otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
    result %<>% ungroup() %>% dplyr::select(-genes)
    
    #double p-value because tmod doesn't 
    result %<>% mutate(P.Value=P.Value*2) 
    
    #adjust again now that identical groups have been merged
    result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value)
    result %<>% mutate(rank = rank(P.Value))
    result %<>% dplyr::select(MainTitle, ID, geneCount = N1, AUC, P.Value, adj.P.Value, everything()) 
    
    if (amino_acid_first_target == "r" & amino_acid_second_target == "k") {
      write_csv(result, here("results", "protein", "Amino_acid_pair_enrichment_r_k.csv"))  
    }
    
    
    #get rank of SRP
    SRP_rank <- result %>% filter(ID == "GO:0006614") %>% head(1) %>% .$rank
    SRP_adj.P.Value <- result %>% filter(ID == "GO:0006614") %>% head(1) %>% .$adj.P.Value
    SRP_AUC <- result %>% filter(ID == "GO:0006614") %>% head(1) %>% .$AUC
    print(SRP_adj.P.Value)
    print(SRP_rank)
    print(SRP_AUC)
    
    symbol_to_sequence_averaged %>% group_by(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0006614`) %>% summarize(percent_target_amino_acids = mean(percent_target_amino_acids))
    
    single_row_result <- data.frame(amino_acid_first = amino_acid_first_target, amino_acid_second = amino_acid_second_target, SRP_AUC = SRP_AUC, SRP_rank = SRP_rank,SRP_adj.P.Value =SRP_adj.P.Value, stringsAsFactors = F)
    all_pair_results <- bind_rows(all_pair_results, single_row_result)
  }
}
Sys.time() - ptm

write_csv(all_pair_results, here("results", "protein", "Amino_acid_pair_enrichment_GO0006614.csv"))
all_pair_results %>% arrange(-SRP_AUC)

all_pair_results <- read_csv(here("results", "protein", "Amino_acid_pair_enrichment_GO0006614.csv"))
all_pair_results_dup <- all_pair_results
all_pair_results_dup %<>% filter(amino_acid_second != amino_acid_first)
all_pair_results_dup$SRP_AUC = NA
all_pair_results_dup %<>% rename(amino_acid_first = amino_acid_second, amino_acid_second = amino_acid_first)
all_pair_results <- bind_rows(all_pair_results, all_pair_results_dup)

ggplot(data = all_pair_results, aes(x = amino_acid_first, y = amino_acid_second)) +
  geom_tile(aes(fill = SRP_AUC), colour='white') +
  ylab("Amino acid") + xlab("Amino acid") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = rev(unique(all_pair_results$amino_acid_first))) +
  scale_fill_distiller(palette = 'RdBu', limits = c(0,1), na.value = "white", name ="ER translocation\ngene enrichment\n(AUC)" ) +
  theme(legend.position = c(0.85, 0.8))
  


