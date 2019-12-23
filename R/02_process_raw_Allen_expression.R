library(reshape2)
library(here)
library(AnnotationDbi)
library(annotate)
library(org.Hs.eg.db)
library(readxl)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readr)

folderPattern <- "normalized_microarray.*"
sampleFilename <- "SampleAnnot.csv"
probeFilename <- "Probes.csv"
expressionFilename <- "MicroarrayExpression.csv"
sourceExpression <- "allen_HBA"

allsampleAnnot = NULL
allExpression = NULL

#fetal brain has more/better gene symbol mappings
probeInfo <- read_csv(here("/data/probe_annotations/fetal_brain_rows_metadata.csv"))
probeInfo %<>% rename(probe_id = probeset_id, probe_name = probeset_name)

for (donorFolder in list.files(here("/data/raw/",sourceExpression,"/"), pattern = folderPattern)) {
  sampleAnnot <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/",sampleFilename))
  
  expressionMatrix <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder, "/",expressionFilename), col_names=F) 
  
  expressionMatrix %<>% rename(probe_id = X1)
  dim(expressionMatrix)
  
  sampleAnnot %<>% mutate(donorID = donorFolder)
  sampleAnnot %<>% mutate(uniqueID = well_id) %>% select(uniqueID, everything())
  
  colnames(expressionMatrix) <- c("probe_id", sampleAnnot$uniqueID)
  
  expressionMatrix <- inner_join(probeInfo %>% select(probe_id, probe_name), expressionMatrix) %>% select(-probe_id)
  
  #bind cols of expression matrix
  allExpression <- bind_cols(allExpression, expressionMatrix)
  
  
  #bind rows of sample annot
  allsampleAnnot <- bind_rows(allsampleAnnot, sampleAnnot)
}

sampleAnnot <- allsampleAnnot
expressionMatrix <- allExpression[, c("probe_name", sampleAnnot$uniqueID)]

print(paste("Number of unique regions:", length(unique(sampleAnnot$structure_name))))
print(paste("Number of samples:", nrow(sampleAnnot)))
print(paste("Number of donors:", length(unique(sampleAnnot$donorID))))


#create probe mapping file here
qc_table <- read_xlsx("./data/probe_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx",skip=1)
qc_table %<>% select(probe_name = "...1", gene_symbol = "...2", is_qc_pass = `Pass?`)

#fix numeric gene symbols from excel date conversions
qc_table <- inner_join( qc_table, probeInfo %>% select(probe_name, for_excel_problem = gene_symbol))
qc_table %<>% mutate(gene_symbol = if_else(!is.na(as.numeric(gene_symbol)), for_excel_problem, gene_symbol)) %>% select(-for_excel_problem)

#use gene symbol to get NCBI ID, then get updated symbol
symbolToID <- probeInfo %>% select(gene_symbol, entrez_id) %>% distinct()
symbolToID %<>% filter(!is.na(entrez_id)) 
symbolToID %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 

## output table for for python code
symbolToIDforPython <- probeInfo %>% select(probe_name, gene_symbol, entrez_id) %>% distinct()
symbolToIDforPython %<>% filter(!is.na(entrez_id)) 
symbolToIDforPython %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 
symbolToIDforPython %<>% mutate(legacySymbol = gene_symbol)
symbolToIDforPython %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% select(gene_symbol, entrez_id,  everything(), -new_symbol)
symbolToIDforPython %<>% select(-legacySymbol)

dir.create(here("data/probe_annotations/"), showWarnings = F, recursive = T)
write_csv(symbolToIDforPython, path=here('data', 'probe_annotations', 'updated_gene_symbol_annotations.csv'))

##

qc_table <- left_join(qc_table, symbolToID)
qc_table %<>% mutate(legacySymbol = gene_symbol) #save old symbol
qc_table %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% dplyr::select(gene_symbol, entrez_id,  everything(), -new_symbol)

#remove those not mapping
qc_table %<>% filter(!grepl("A_", gene_symbol)) %>% filter(!grepl("CUST_", gene_symbol)) 
qc_table %<>% filter(gene_symbol != "na")

print(paste("Probes after filtering for _ and numeric geneSymbols", nrow(qc_table)))
print(paste("Gene count",length(unique(qc_table$gene_symbol))))

#just probes that pass
qc_table %<>% filter(is_qc_pass == TRUE) 
print(paste("Probes after filtering for Miller qc_pass", nrow(qc_table)))
print(paste("Gene count",length(unique(qc_table$gene_symbol))))

#write out qc table for NCBI ID to gene symbol mapping
symbol_to_ids <- qc_table %>% select(gene_symbol, entrez_id) %>% distinct()
entrez_to_ensembl <- org.Hs.egENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(entrez_to_ensembl)
# Convert to a list
entrez_to_ensembl_list <- as.list(entrez_to_ensembl[mapped_genes])
symbol_to_ids %<>% rowwise() %>% mutate(ensembl_ids = paste(unlist(entrez_to_ensembl_list[as.character(entrez_id)]), collapse = "|"))
dir.create(here("data/R_processed_expression_data/"), showWarnings = F, recursive = T)
symbol_to_ids %>% write_csv(here("data/R_processed_expression_data/","symbol_to_IDs_table.csv"))

#melt to long form - can be slow
expressionMatrix_melted <- melt(expressionMatrix, variable.name = "well_id", value.name = "expression") %>% as_tibble() %>% mutate(well_id = as.character(well_id))
expressionMatrix <- NULL

#add in gene symbols
expressionMatrix_melted <- inner_join(qc_table %>% select(probe_name, gene_symbol), expressionMatrix_melted)


#mean average for a gene and well_id/sample
#warnings - takes a lot of memory/cpu, can be done on a mac with 16gb of RAM
expressionMatrix_melted

#write out probe level data for Marker Gene Profile
expressionMatrix_melted %>% write_csv(here("data/R_processed_expression_data/","microarray_expression_probe_level_qcNames.csv"))

expressionMatrix_melted %<>% select(-probe_name)
expressionMatrix_melted %<>% group_by(well_id, gene_symbol) 
expressionMatrix_melted %<>% summarize(expression = mean(expression))

expressionMatrix_melted %>% write_csv(here("data/R_processed_expression_data/","microarray_expression_probe_averaged_qcNames.csv"))
