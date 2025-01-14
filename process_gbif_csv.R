library(arrow)
library(dplyr)
library(readr)

butterfly_path <- file.path("D:", "data", "gbif_observation_org_butterflies", "gbif_butterfly_observation-org", "occurences.csv")
file.exists(butterfly_path)
d <- arrow::open_delim_dataset(butterfly_path, delim='\t')
d %>% filter(year==2023, countryCode=="NL", order=="Lepidoptera") %>% collect() -> df
dim(df)
length(unique(df$species))
save_path <- file.path("D:", "data", "gbif_observation_org_butterflies", "gbif_butterfly_observation-org", "gbif_subset_netherlands_lepidoptera.csv")
write_delim(df, file=save_path, delim="\t")

