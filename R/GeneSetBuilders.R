library(readr)
library(tmod)
library(annotate)
library(GO.db)
library(org.Hs.eg.db)
library(here)

otherGeneListsFolder <- here("data", "gene_lists")

loadGOSets <- function(geneBackground) {
  go_object <- as.list(org.Hs.egGO2ALLEGS)
  
  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')
  
  #build GO sets for tmod -slow
  tmodNames <- data.frame()
  modules2genes <- list()
  goGroupName <- names(go_object)[1]
  showMethods(Term)
  
  goCount <- length(go_object)
  count <- 1
  for(goGroupName in names(go_object)) {
    if (count %% 1000 == 0) print(paste(count, "of", goCount))
    count <- count + 1
    
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
    
    genesymbols <- intersect(genesymbols, geneBackground)
    if (!(length(genesymbols) >= 10 & length(genesymbols) <= 200)) next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}

loadFileSets <- function(prefix, isMouse = FALSE) {
  tmodNames <- data.frame()
  modules2genes <- list()
  
  for(geneListFilename in list.files(otherGeneListsFolder, pattern = paste0(prefix,".*txt"), full.names = T)) {
    print(geneListFilename)
    genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
    shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
    shortName <- gsub(paste0(prefix,"."),"",shortName)
    
    genesOfInterest$term <- shortName
    
    if (isMouse) {
      modules2genes[shortName] <- list(mouse2human(genesOfInterest$V1)$humanGene)
    } else {
      modules2genes[shortName] <- list(genesOfInterest$V1)
    }
    
    tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
  }
  geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
  geneSets
}



