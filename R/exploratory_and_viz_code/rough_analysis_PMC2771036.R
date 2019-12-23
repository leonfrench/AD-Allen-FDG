library(here)
library(readr)
library(tidyr)
library(magrittr)
library(readxl)

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2771036/#S15
#Additional file 15:  
#from GSE16134

combined <- as_tibble()
for(bacteria_correlations in list.files(here("data", "Papapanou_et_al.", "1471-2180-9-221-S15"), pattern = "xls", full.names = T)) {
  print(bacteria_correlations)
  table <- read_xls(bacteria_correlations)
  table %<>% mutate(rank = dplyr::row_number())
  table %<>% filter(ID == "GO:0006613") #cotranslational protein targeting to membrane
  table %<>% mutate(file = bacteria_correlations)
  combined <- bind_rows(combined , table)
}
combined %>% select(`p-value`, file, rank) %>% arrange(`p-value`)

#'health-associated' bacterial species are A. naeslundii and V. parvula


