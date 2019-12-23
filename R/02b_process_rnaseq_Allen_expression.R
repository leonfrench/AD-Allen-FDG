library(reshape2)
library(here)
library(AnnotationDbi)
library(annotate)
library(org.Hs.eg.db)
library(readxl)
detach("package:dplyr", unload = TRUE)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readr)

folderPattern <- "rnaseq*"
sampleFilename <- "SampleAnnot.csv"
ontologyFilename <- "Ontology.csv"
expressionFilename <- "RNAseqTPM.csv"
sourceExpression <- "allen_HBA"

allsampleAnnot = NULL
allExpression = NULL

#fetal brain has more/better gene symbol mappings
probeInfo <-
  read_csv(here("/data/probe_annotations/fetal_brain_rows_metadata.csv"))
probeInfo %<>% rename(probe_id = probeset_id, probe_name = probeset_name)

for (donorFolder in list.files(here("./data/raw/", sourceExpression, "/"), pattern = folderPattern)) {
  sampleAnnot <-
    read_csv(here(
      "./data/raw/",
      sourceExpression,
      donorFolder,
      sampleFilename
    ))
  ontology <-
    read_csv(here(
      "./data/raw/",
      sourceExpression,
      donorFolder,
      ontologyFilename
    ))
  ontology %<>% rename(structure_name = name)
  # merge in ontology to sampleAnnot here so that you can get the full structure_name
  sampleAnnot %<>% left_join(ontology %>% select(id, structure_name),
                             by = c('ontology_structure_id' = 'id'))
  
  expressionMatrix <- read_csv(
    here(
      "./data/raw/",
      sourceExpression,
      donorFolder,
      expressionFilename
    ),
    col_names = F
  )
  
  expressionMatrix %<>% rename(gene_symbol = X1)
  dim(expressionMatrix)
  
  sampleAnnot %<>% mutate(donorID = donorFolder)
  #should merge in ontology to sampleAnnot here so that you can get the full structure_name
  sampleAnnot %<>% mutate(uniqueID = well_id) %>% select(uniqueID, everything())
  
  colnames(expressionMatrix) <-
    c("gene_symbol", sampleAnnot$uniqueID)
  
  #bind cols of expression matrix
  allExpression <- bind_cols(allExpression, expressionMatrix)
  
  #bind rows of sample annot
  allsampleAnnot <- bind_rows(allsampleAnnot, sampleAnnot)
}

sampleAnnot <- allsampleAnnot
expressionMatrix <-
  allExpression[, c("gene_symbol", sampleAnnot$uniqueID)]

print(paste("Number of unique regions:", length(unique(
  sampleAnnot$structure_name
))))
print(paste("Number of samples:", nrow(sampleAnnot)))
print(paste("Number of donors:", length(unique(
  sampleAnnot$donorID
))))


#create probe mapping file here
qc_table <-
  read_xlsx(
    here(
      "./data/probe_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx"
    ),
    skip = 1
  )
qc_table %<>% select(probe_name = "...1",
                     gene_symbol = "...2",
                     is_qc_pass = `Pass?`)


#fix numeric gene symbols from excel date conversions
qc_table %>% filter(!is.na(as.numeric(gene_symbol)))

qc_table <-
  inner_join(qc_table,
             probeInfo %>% select(probe_name, for_excel_problem = gene_symbol))
qc_table %<>% mutate(gene_symbol = if_else(!is.na(as.numeric(gene_symbol)), for_excel_problem, gene_symbol)) %>% select(-for_excel_problem)

#use Genes.csv to get NCBI ID, then get updated symbol
symbolToID <- read_csv(here("data", "raw", "allen_HBA", donorFolder, "Genes.csv")) %>% select(gene_symbol, entrez_id) %>% distinct()
symbolToID %<>% filter(!is.na(entrez_id))
symbolToID %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data = 'org.Hs.eg'))

qc_table <- left_join(qc_table, symbolToID)
qc_table %<>% mutate(legacySymbol = gene_symbol) #save old symbol
qc_table %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% dplyr::select(gene_symbol, entrez_id,  everything(),-new_symbol)

#remove those not mapping
qc_table %<>% filter(gene_symbol != "na")

print(paste(
  "Probes after filtering for _ and numeric geneSymbols",
  nrow(qc_table)
))
print(paste("Gene count", length(unique(
  qc_table$gene_symbol
))))

#just probes that pass
qc_table %<>% filter(is_qc_pass == TRUE)
print(paste("Probes after filtering for Miller qc_pass", nrow(qc_table)))
print(paste("Gene count", length(unique(
  qc_table$gene_symbol
))))

#melt to long form - can be slow
expressionMatrix_melted <-
  melt(expressionMatrix,
       variable.name = "well_id",
       value.name = "expression") %>% as_tibble() %>% mutate(well_id = as.character(well_id))
expressionMatrix <- NULL

#update gene symbols
expressionMatrix_melted %<>% rename(legacySymbol = gene_symbol)
expressionMatrix_melted <- inner_join(expressionMatrix_melted, qc_table %>% select(legacySymbol, gene_symbol)) 
expressionMatrix_melted %<>% mutate(gene_symbol = if_else(is.na(gene_symbol), legacySymbol, gene_symbol)) %>% select(gene_symbol, well_id, expression)

#mean average for a gene and well_id/sample
length(unique(expressionMatrix_melted$gene_symbol))

expressionMatrix_melted %<>% group_by(well_id, gene_symbol)
expressionMatrix_melted %<>% summarize(expression = mean(expression))

dir.create(
  here("data/R_processed_expression_data/"),
  showWarnings = F,
  recursive = T
)

expressionMatrix_melted %>% write_csv(here(
  "data","R_processed_expression_data",
  "rnaseq_averaged_qcNames.csv"
))
