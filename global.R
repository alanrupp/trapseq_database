library(dplyr)
library(readr)

# - Read in data information ---------------------------------------------------
data_files <- list.files(path = "data", pattern = ".csv")

# make list of either hypo or hindbrain
areas_hypo <- c('Hypothalamus','Arc','VMH')
areas_hb <- c('Hindbrain','NTS')

regions <- data.frame(
  "subregion" = str_extract(data_files, "(.+)(?=_)")
) %>%
  mutate(region = case_when(
    subregion %in% areas_hypo ~ "Hypothalamus",
    subregion %in% areas_hb ~ "Hindbrain"
  ))


# Make options list of brain areas and genes based on data
pulldown_list <- unique(pops)
regions <- data.frame(subregion = unique(areas))


regions <- regions %>%
  mutate(region = case_when(
    subregion %in% areas_hypo ~ "Hypothalamus",
    subregion %in% areas_hb ~ "Hindbrain"
  ))

region_list <- 
  list(Hypothalamus = filter(regions, region == 'Hypothalamus')$subregion,
       Midbrain = filter(regions, region == 'Midbrain')$subregion,
       Hindbrain = filter(regions, region == 'Hindbrain')$subregion)

# all genes
genes <- read_csv("gene_biotypes.csv")

gene_list <- unique(genes$Gene)

biotypes <-
  genes %>%
  group_by(Biotype) %>%
  tally(.) %>%
  arrange(desc(n)) %>%
  select(-n)

biotypes_list <- biotypes$Biotype
