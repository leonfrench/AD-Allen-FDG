########################
#Venn Diagram
#install.packages('venneuler')

library(dplyr)
library(here)
library(venneuler)
library(ggrepel)
detach("package:AnnotationDbi", unload = TRUE)
detach("package:org.Hs.eg.db", unload = TRUE)
source(here("R", "GeneSetBuilders.R"))


#input is a tmod result table and the number of top groups - don't send the whole result
getEulerDiagram <- function(goResult, geneSets) {
  combined <- NULL
  for(group in goResult$ID) {
    groupName <- as.character(filter(geneSets$MODULES, ID == group)$Title[1])
    combined <- rbind(combined, data.frame(elements= unlist(geneSets$MODULES2GENES[group]), sets=groupName))
  }
  
  v <- venneuler(combined)
  plot(v)
  
  (vForGGplot <- tbl_df(data.frame(diameter=v$diameters, v$centers, color=v$colors, goName=v$labels, stringsAsFactors = F)))
  
  #credited to user5061 and baptiste via stackoverflow.com
  circularise <- function(d, n=360){
    angle <- seq(-pi, pi, length = n)
    make_circle <- function(x,y,r,goName){
      data.frame(x=x+r*cos(angle), y=y+r*sin(angle), goName)
    }
    lmat <- mapply(make_circle, goName = d[,"goName"], 
                   x = d[,"x"], y=d[,"y"], r=d[,"diameter"]/2, SIMPLIFY = FALSE)
    do.call(rbind, lmat)
  }
  
  circles <- circularise(vForGGplot)
  vForGGplot <- as.data.frame(vForGGplot)
  diagram <- ggplot(vForGGplot) + geom_blank(aes(x, y)) + 
    geom_polygon(aes(x,y, group=goName, fill=goName), data=circles, alpha=0.4) +
    coord_fixed() + theme_void(base_size = 16) + theme(legend.position="none") +
    geom_label_repel(aes(x, y, fill=goName, label = goName), color="black", segment.alpha=0, box.padding = unit(0.2, "lines")) 
  diagram
}


allGenes <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv"))
allGenes <- allGenes$gene_symbol
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  geneSetsGO <- loadGOSets(allGenes)
}
geneSetsToUse <- geneSetsGO

length(intersect(geneSetsGO$MODULES2GENES$`GO:0006614`,geneSetsGO$MODULES2GENES$`GO:0006613`))
setdiff(geneSetsGO$MODULES2GENES$`GO:0006613`,geneSetsGO$MODULES2GENES$`GO:0006614`)
length(intersect(geneSetsGO$MODULES2GENES$`GO:0048167`,geneSetsGO$MODULES2GENES$`GO:0099504`))
#SRP and non-sense decay
length(intersect(geneSetsGO$MODULES2GENES$`GO:0000184`,geneSetsGO$MODULES2GENES$`GO:0006614`))
setdiff(geneSetsGO$MODULES2GENES$`GO:0000184`,geneSetsGO$MODULES2GENES$`GO:0006614`)



geneSetsGO$MODULES2GENES$`GO:0044450`
#learning and synaptic vesicle cycle intersect
length(intersect(geneSetsGO$MODULES2GENES$`GO:0099504`,geneSetsGO$MODULES2GENES$`GO:0007612`))



#results <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/GO_results.csv"))
results <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_12876_14380_15496_15697_9861/GO_results.csv"))
results %<>% mutate(directioned_rank = rank(-1*(sign(AUC-.5) * (max(rank) - rank))))
results %<>% arrange(directioned_rank)

getEulerDiagram(results %>% head(10), geneSetsToUse)

getEulerDiagram(results %>% tail(10), geneSetsToUse)


geneStatistics <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_10021/geneStatistics.csv"))
geneStatistics %>% filter(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0044450`) %>% tail()
geneStatistics %>% filter(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0030897`) %>% print()
#decay
geneStatistics %>% filter(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0000184`) %>% head(20) %>% as.data.frame()
geneStatistics %>% filter(gene_symbol %in% setdiff(geneSetsGO$MODULES2GENES$`GO:0000184`,geneSetsGO$MODULES2GENES$`GO:0006614`))

#ribosome in all 5
geneStatistics <- read_csv(here("/results/microarray/Allen_Schwartz_NiftiValue/donors_12876_14380_15496_15697_9861/geneStatistics.csv"))
geneStatistics %>% filter(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0003735`) %>% tail(20) %>% print(20)


geneSetsToUse <- geneSetsGO
targetGroupID <- dplyr::filter(tbl_df(geneSetsToUse$MODULES), Title == "SRP-dependent cotranslational protein targeting to membrane")$ID
