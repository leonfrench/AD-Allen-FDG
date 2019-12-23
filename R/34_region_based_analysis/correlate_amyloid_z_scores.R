library(reshape2)
library(magrittr)
library(readr)
library(annotate)
library(GO.db)
library(org.Hs.eg.db)
library(tmod)
library(here)
library(dplyr)

freesurfer_six_regions <- read_tsv(here("/data/Grothe_et_al./DKRegionStatistics.tsv"))
freesurfer_six_regions %<>% rename(region=X1)

freesurfer_six <- read_tsv(here("/data/Grothe_et_al./AllenHBA_DK_ExpressionMatrix.tsv"))
freesurfer_six %<>% rename(Gene=X1)
freesurfer_six %<>% select(-`Average donor correlation to median`)
freesurfer_six <- as_tibble(melt(freesurfer_six, variable.name = "region"))

amyloid <- read_csv(here("/data/Grothe_et_al./Amyloid dep.estimates.edits.csv")) %>% select(region, Amyloid_deposition_z)
amyloid %>% arrange(-Amyloid_deposition_z)
freesurfer_six_regions %>% inner_join(amyloid)

freesurfer_six %<>% inner_join(amyloid)

freesurfer_six_correlations <- freesurfer_six %>% group_by(Gene) %>% summarize(r = cor(value, Amyloid_deposition_z, m='s'))

#agrees with Grothe et al. 
freesurfer_six_correlations %>% filter(Gene == "MAPT")
freesurfer_six_correlations %>% filter(Gene == "APP")

freesurfer_six_correlations %>% arrange(r)

freesurfer_six_correlations %>% filter(Gene == "RPS28")
freesurfer_six_correlations %>% filter(Gene == "RPS27A")
freesurfer_six_correlations %>% filter(Gene == "RPS25")
freesurfer_six_correlations %>% filter(Gene == "RPL34")

ggplot(freesurfer_six %>% filter(Gene == "RPS28"), aes(value, Amyloid_deposition_z)) + geom_point()

#look in single brain 10021
freesurfer_10021 <- read_tsv(here("/data/Grothe_et_al./AllenHBA_DK_ExpessionMatrix.10021.tsv"))
freesurfer_10021 %<>% rename(Gene=X1)
freesurfer_10021 %<>% select(-`Average donor correlation to median`)
freesurfer_10021 <- as_tibble(melt(freesurfer_10021, variable.name = "region"))
freesurfer_10021 %<>% inner_join(amyloid)

ggplot(freesurfer_10021 %>% filter(Gene == "RPS28"), aes(value, Amyloid_deposition_z)) + geom_point()
freesurfer_10021_correlations <- freesurfer_10021 %>% group_by(Gene) %>% summarize(r = cor(value, Amyloid_deposition_z, m='s', use="complete.obs"))

freesurfer_10021_correlations %>% filter(Gene == "MAPT")
freesurfer_10021_correlations %>% filter(Gene == "APP")
freesurfer_10021_correlations %>% filter(Gene == "RPS28")
freesurfer_10021_correlations %>% filter(Gene == "RPL34")
freesurfer_10021_correlations %>% filter(Gene == "RPS27A")
freesurfer_10021_correlations %>% filter(Gene == "RPS25")
freesurfer_10021_correlations %>% filter(Gene == "CD33")

Grothe_top_genes <- read_csv(here("/data/Grothe_et_al./Grothe.SupplementTable1.page1.genes.txt"), col_names = F) %>% rename(Gene = X1)
freesurfer_10021_correlations_top <- inner_join(freesurfer_10021_correlations, Grothe_top_genes)
freesurfer_six_correlations_top <- inner_join(freesurfer_six_correlations, Grothe_top_genes)

median(freesurfer_10021_correlations_top$r)
median(freesurfer_six_correlations_top$r)

joined <- inner_join(freesurfer_six_correlations_top, freesurfer_10021_correlations_top, by="Gene", suffix = c("_six", "_brain_10021")) %>% mutate(diff = abs(r_six-r_brain_10021))

joined %>% arrange(-diff)
joined %>% arrange(r_six)
wilcox.test(joined$r_six, joined$r_brain_10021)

median(freesurfer_10021_correlations_top$r)
median(freesurfer_six_correlations_top$r)

#load GO group for SRP, run AUC analyses


modules2genes <- list()
go_object <- as.list(org.Hs.egGO2ALLEGS)
GO_0006614_symbols <- unname(getSYMBOL(unique(unlist(go_object$`GO:0006614`)), data='org.Hs.eg'))

tmodNames <- data.frame(ID="GO_0006614_symbols", Title = "GO_0006614_symbols")
modules2genes["GO_0006614_symbols"] <- list(GO_0006614_symbols)
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)

sortedGenes <- freesurfer_six_correlations %>% arrange(-r) %>% .$Gene
tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)
all_six_AUC <- tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)["GO_0006614_symbols", "AUC"]

sortedGenes <- freesurfer_10021_correlations %>% arrange(-r) %>% .$Gene
tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)
intersect(sortedGenes, geneSets$MODULES2GENES$GO_0006614_symbols)

median(freesurfer_six_correlations %>% filter(Gene %in% geneSets$MODULES2GENES$GO_0006614_symbols) %>% .$r)
median(freesurfer_10021_correlations %>% filter(Gene %in% geneSets$MODULES2GENES$GO_0006614_symbols) %>% .$r)

#iterate single brains, duplicated code
for (file in list.files(here("/data/Grothe_et_al./"), pattern="AllenHBA_DK_ExpessionMatrix.*.tsv", full.names = TRUE)) {
  freesurfer_single <- read_tsv(file)
  freesurfer_single %<>% rename(Gene=X1)
  freesurfer_single %<>% select(-`Average donor correlation to median`)
  freesurfer_single <- as_tibble(melt(freesurfer_single, variable.name = "region"))
  freesurfer_single %<>% inner_join(amyloid)
  freesurfer_single_correlations <- freesurfer_single %>% group_by(Gene) %>% summarize(r = cor(value, Amyloid_deposition_z, m='s', use="complete.obs"))
  
  sortedGenes <- freesurfer_single_correlations %>% arrange(-r) %>% .$Gene
  print(file)
  #these p-values need to be doubled
  print(tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F) %>% mutate(P.Value = P.Value*2))
  print(tmodUtest(sortedGenes, mset=geneSets, qval = 1, filter = F)["GO_0006614_symbols", "AUC"] - all_six_AUC)
}
