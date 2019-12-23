library(here)
library(dplyr)
#obtained with WebPlotDigitizer
figure2_data_APP <- read_csv(here("/data/Grothe_et_al./Amyloid dep vrs APP.Grothe et al.from Figure 2.csv"))
figure2_data_MAPT <- read_csv(here("/data/Grothe_et_al./Amyloid dep vrs MAPT.Grothe et al.from Figure 2.csv"))

freesurfer_six_regions <- read_tsv(here("/data/Grothe_et_al./DKRegionStatistics.tsv"))
freesurfer_six_regions %<>% rename(region=X1)

freesurfer_six <- read_tsv(here("/data/Grothe_et_al./AllenHBA_DK_ExpressionMatrix.tsv"))
freesurfer_six %<>% rename(Gene=X1)
freesurfer_six %<>% select(-`Average donor correlation to median`)
freesurfer_six <- as_tibble(melt(freesurfer_six, variable.name = "region"))

freesurfer_six_MAPT <- freesurfer_six %>% filter(Gene == "MAPT")
freesurfer_six_APP <- freesurfer_six %>% filter(Gene == "APP")
freesurfer_six_MAPT <- inner_join(freesurfer_six_regions %>% filter(Hemisphere == "lh") %>% select(region, BrainsWithSamples), freesurfer_six_MAPT)
freesurfer_six_APP <- inner_join(freesurfer_six_regions %>% filter(Hemisphere == "lh") %>% select(region, BrainsWithSamples), freesurfer_six_APP)

freesurfer_six_APP %<>% arrange(value) %>% rename(APP_freesurfer = value)
freesurfer_six_MAPT %<>% arrange(value) %>% rename(MAPT_freesurfer = value)

#sort from low to high then merge both
figure2_data_APP %<>% arrange(APP_expression)
merged_APP <- bind_cols(freesurfer_six_APP, figure2_data_APP) %>% arrange(Amyloid_deposition_z)

#sort from low to high then merge both
figure2_data_MAPT %<>% arrange(MAPT_expression)
merged_MAPT <- bind_cols(freesurfer_six_MAPT, figure2_data_MAPT) %>% arrange(Amyloid_deposition_z)

#some disagreement
bind_cols(merged_MAPT %>% select(region_MAPT = region, Amyloid_deposition_z, MAPT_expression), merged_APP %>% select(region_APP = region, Amyloid_deposition_z, APP_expression))
#inner join and average to get a single amyloid value
merged_full <- inner_join(merged_MAPT %>% select(region, Amyloid_deposition_z_MAPT = Amyloid_deposition_z, MAPT_expression), merged_APP %>% select(region, Amyloid_deposition_z_APP = Amyloid_deposition_z, APP_expression), by=c("region"))
cor.test(merged_full$Amyloid_deposition_z_MAPT, merged_full$Amyloid_deposition_z_APP)
plot(merged_full$Amyloid_deposition_z_MAPT, merged_full$Amyloid_deposition_z_APP)
merged_full %<>% mutate(diff = abs(Amyloid_deposition_z_MAPT - Amyloid_deposition_z_APP)) %>% arrange(-diff)
merged_full %>% write_csv(here("/data/Grothe_et_al./Amyloid dep.estimates.csv"))

#manual edits to resolve diagreements according to gene exp differences - the data from the gene with the largest difference was selected for amyloid data (for disagreeing pairs)
merged_full <- read_csv(here("/data/Grothe_et_al./Amyloid dep.estimates.edits.csv")) %>% select(region, Amyloid_deposition_z)
merged_full %>% arrange(Amyloid_deposition_z) 


