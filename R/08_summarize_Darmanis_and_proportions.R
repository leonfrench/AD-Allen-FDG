library(reshape2)
library(readr)
library(dplyr)
library(magrittr)
library(here)
library(ggplot2)

#########################################
#just marker gene enrichment first
donor_10021 <- read_csv(here("results/microarray/Allen_Schwartz_NiftiValue/donors_10021/","Darmanis_custom_genelists_results.csv"))
five_brains <- read_csv(here("results/microarray/Allen_Schwartz_NiftiValue/donors_9861_10021_12876_14380_15496_15697/","Darmanis_custom_genelists_results.csv"))
donor_10021 %<>% select(-ID, -U, -P.Value)
five_brains %<>% select(-ID, -U, -P.Value, -N1)
combined <- inner_join(five_brains, donor_10021, by=c("Title"), suffix=c(".five_brains", ".10021"))
combined %<>% mutate(AUC_difference = AUC.10021-AUC.five_brains)
combined %<>% select(Type = Title, 'Genes' = N1, everything())

combined %<>% mutate(AUC.10021 = signif(AUC.10021, digits=2), 
                     AUC.five_brains = signif(AUC.five_brains, digits=2),
                     adj.P.Val.five_brains = signif(adj.P.Val.five_brains, digits=2),
                     adj.P.Val.10021 = signif(adj.P.Val.10021, digits=2),
                     AUC_difference = signif(AUC_difference, digits=2)) %>% arrange(-abs(AUC_difference)) 
colnames(combined) <- gsub("adj.P.Val.", "pFDR ", colnames(combined))
colnames(combined) <- gsub("[._]", " ", colnames(combined))
combined %>% write_csv("results/microarray/Allen_Schwartz_NiftiValue/donors_10021/Darmanis_custom_genelists_combined_table.csv")

#########################################
#linear models on estimated proportions
all_donors_estimates <- read_csv(here("./data/Marker_Gene_Profile/", "Darmanis_markers_cortex_only.csv"))

probeLevelColAnnotations <- read_csv(here("/data/Combined_SampleAnnot_with_PET.csv"))
probeLevelColAnnotations %<>% mutate(well_id = as.character(well_id)) %>% select(well_id, donor, structure_name, Allen_Schwartz_NiftiValue)
pet_threshold <- 0.006
probeLevelColAnnotations %<>% mutate(inside = Allen_Schwartz_NiftiValue > pet_threshold) 

all_donors_estimates %<>% mutate(well_id = as.character(well_id))
all_donors_estimates %<>% melt() %>% as_tibble()
all_donors_estimates <- inner_join(probeLevelColAnnotations, all_donors_estimates)
all_donors_estimates %<>% mutate(variable = gsub("Darmanis.", "", variable))

unique(all_donors_estimates$variable)
all_donors_estimates %>% filter(variable == "Microglia") %>% group_by(donor, inside) %>% summarize(n=dplyr::n(), meanMicro = mean(value))
all_donors_estimates %<>% mutate(donor = as.factor(donor))

all_donors_estimates %<>% mutate(is10021 = donor == "10021")

all_donors_estimates %<>% filter(!is.na(value))
for(celltype in unique(all_donors_estimates$variable)) {
  print(celltype)
  summary(lm(value ~ inside * is10021 + donor, data = all_donors_estimates %>% filter(variable == celltype)))
  z <- summary(lm(value ~ inside * is10021 + donor, data = all_donors_estimates %>% filter(variable == celltype)))
  print(z$coefficients["insideTRUE:is10021TRUE", "Pr(>|t|)"])
}

summary(lm(value ~ inside * is10021 + donor, data = all_donors_estimates %>% filter(variable == "Microglia")))
summary(lm(value ~ inside * is10021 + donor, data = all_donors_estimates %>% filter(variable == "Neuron")))

#plot
donors_10021_estimates <- all_donors_estimates #%>% filter(donor == "10021")
donors_10021_estimates %<>% group_by(variable, inside, donor) %>% summarize(average_prop = mean(value))
ggplot(data = donors_10021_estimates, aes(x=variable, y=average_prop, fill=inside)) + 
  geom_bar(stat="identity", position=position_dodge()) + facet_wrap( ~ donor)

