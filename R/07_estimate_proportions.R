library(reshape2)
library(readr)
library(dplyr)
library(magrittr)
library(here)
library(homologene)
library(markerGeneProfile)

pet_threshold <- 0.006
probeLevelExpression <- read_csv(here("data/R_processed_expression_data/","microarray_expression_probe_level_qcNames.csv"))
probeLevelColAnnotations <- read_csv(here("/data/Combined_SampleAnnot_with_PET.csv"))

probe_to_gene_map <- probeLevelExpression %>% select(probe_name, gene_symbol) %>% distinct()
all_donors_estimates_Darm = NULL #initialize to blank - does each donor individually
all_donors_estimates_Expresso = NULL
for (target_donor in unique(probeLevelColAnnotations$donor)) {
  print(target_donor)
  samplesForDonor <- probeLevelColAnnotations %>% filter(target_donor == donor)
  
  #filter cortex
  samplesForDonor %<>% filter(in_cortex_excluding_piriform_hippocampus == TRUE)
  
  #cast long form matrix
  target_donor_expression <- probeLevelExpression %>% filter(well_id %in% samplesForDonor$well_id)
  target_donor_expression %<>% select(-gene_symbol)
  target_donor_expression <- as_tibble(dcast(target_donor_expression, probe_name ~ well_id))
  
  
  #join back to get gene symbols
  expressionMatrix <- inner_join(probe_to_gene_map, target_donor_expression)
  expressionMatrix %<>% select( Probe = probe_name, Gene.Symbol = gene_symbol,everything())
  
  medExp = expressionMatrix %>% 
    sepExpr() %>% {.[[2]]} %>%
    unlist %>% median
  #medExp <- 0
  expressionMatrix = mostVariable(expressionMatrix, threshold = medExp, threshFun=median)
  
  #setup gene lists
  otherGeneListsFolder <- here("/data/gene_lists/")
  DarmanisLists <- list()
  for(geneListFilename in list.files(otherGeneListsFolder, pattern = "Darmanis.*txt", full.names = T)) {
    print(geneListFilename)
    genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
    genesOfInterest <- genesOfInterest$V1
    genesOfInterest <- intersect(genesOfInterest, expressionMatrix$Gene.Symbol)
    shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
    DarmanisLists[shortName] <- list(genesOfInterest)
    print(length(genesOfInterest))
  }
  
  ExpressoLists <- list()
  for(geneListFilename in list.files(otherGeneListsFolder, pattern = "NeuroExpresso.Cortex.*txt", full.names = T)) {
    print(geneListFilename)
    genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
    genesOfInterest <- genesOfInterest$V1
    genesOfInterest <- mouse2human(genesOfInterest)$humanGene
    genesOfInterest <- intersect(genesOfInterest, expressionMatrix$Gene.Symbol)
    shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
    ExpressoLists[shortName] <- list(genesOfInterest)
    print(length(genesOfInterest))
  }
  
  #ensure columns are in sorted order - by character sorting
  samplesForDonor %<>% mutate(well_id = as.character(well_id))
  expressionMatrix <- expressionMatrix[,sort(colnames(expressionMatrix))]
  expressionMatrix %<>% select(Gene.Symbol, Probe, everything())

  samplesForDonor %<>% arrange(well_id)
  identical(colnames(expressionMatrix %>% select(-Gene.Symbol, -Probe)), samplesForDonor$well_id)
  samplesForDonor %<>% mutate(inside = Allen_Schwartz_NiftiValue > pet_threshold) 
  
  
  #this is code from Ogan 
  estimationsDarm =  mgpEstimate(exprData=expressionMatrix,
                                 genes=DarmanisLists,
                                 geneColName='Gene.Symbol',
                                 outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                                 geneTransform = function(x){ x }, # this is the default option for geneTransform
                                 groups=samplesForDonor$inside, #if there are experimental groups provide them here. if not desired set to NULL
                                 seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                                 removeMinority = TRUE) # removes genes if they are the minority in terms of rotation sign from estimation process
  
  
  estimationsDarm <- as.data.frame(estimationsDarm$estimates)
  estimationsDarm$well_id <- rownames(estimationsDarm)
  estimationsDarm %<>% as_tibble() %>% select(well_id, everything())
  all_donors_estimates_Darm <- bind_rows(estimationsDarm, all_donors_estimates_Darm)

  if (length(ExpressoLists) != 0) {
    estimationsExpresso =  mgpEstimate(exprData=expressionMatrix,
                                       genes=ExpressoLists,
                                       geneColName='Gene.Symbol',
                                       outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                                       geneTransform = function(x){ x }, # this is the default option for geneTransform
                                       groups=samplesForDonor$inside, #if there are experimental groups provide them here. if not desired set to NULL
                                       seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                                       removeMinority = TRUE) # removes genes if they are the minority in terms of rotation sign from estimation process
    estimationsExpresso <- estimationsExpresso$estimates
    estimationsExpresso[names(estimationsExpresso[lapply(estimationsExpresso, length) == 0])] <- NULL
    estimationsExpresso <- as.data.frame(estimationsExpresso)
    estimationsExpresso$well_id <- rownames(estimationsExpresso)
    estimationsExpresso %<>% as_tibble() %>% select(well_id, everything())
    all_donors_estimates_Expresso <- bind_rows(estimationsExpresso, all_donors_estimates_Expresso)
  }
}

dir.create(here("./data/Marker_Gene_Profile/"), showWarnings = FALSE)
all_donors_estimates_Darm %>% write_csv(here("./data/Marker_Gene_Profile/", "Darmanis_markers_cortex_only.csv"))
all_donors_estimates_Expresso %>% write_csv(here("./data/Marker_Gene_Profile/", "Neuroexpresso_markers_cortex_only.csv"))
