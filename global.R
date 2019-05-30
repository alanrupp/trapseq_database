library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

# - Read in data information ---------------------------------------------------
data_files <- list.files(path = "data", pattern = ".csv")
pulldowns <- data.frame("pulldown" = str_extract(data_files, ".+(?=\\.csv)")) %>%
  separate(pulldown, into = c("subregion", "pop"), sep = "_")
  
# determine if it's either hypothalamus or hindbrain
areas_hypo <- c('Hypothalamus', 'ARH', 'VMH')
areas_hb <- c('Hindbrain', 'NTS')
pulldowns <- mutate(pulldowns,
                    region = case_when(
                      subregion %in% areas_hypo ~ "Hypothalamus",
                      subregion %in% areas_hb ~ "Hindbrain"
                    ))

# turn into a list for user selection
pulldown_list <- 
  map(unique(pulldowns$region), ~ paste(
    filter(pulldowns, region == .x)$pop,
    filter(pulldowns, region == .x)$subregion))
names(pulldown_list) <- unique(pulldowns$region)



# all genes
genes <- read_csv("gene_biotypes.csv")
gene_list <- unique(genes$Gene)

biotypes <-
  genes %>%
  group_by(Biotype) %>%
  tally(.) %>%
  arrange(desc(n)) %>%
  .$Biotype
