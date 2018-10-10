library(dplyr)

# - Read in data --------------------------------------------------------------
data_files <- list.files(path = "Enrich_data", 
                         pattern = ".csv")
all_data <- lapply(paste0("Enrich_data/", data_files), 
                   readr::read_csv, 
                   col_types = readr::cols())
pops <- gsub("(.+)_(.+)\\.csv", "\\2", data_files)
areas <- gsub("(.+)_(.+)\\.csv", "\\1", data_files)

datasets <- paste(pops, areas)
names(all_data) <- datasets


# Make options list of brain areas and genes based on data
pulldown_list <- unique(pops)
regions <- data.frame(subregion = unique(areas))
areas_hypo <- c('Hypothalamus','Arc','VMH')
areas_hb <- c('Hindbrain','NTS')

regions <- regions %>%
  mutate(region = case_when(
    subregion %in% areas_hypo ~ "Hypothalamus",
    subregion %in% areas_hb ~ "Hindbrain"
  ))

region_list <- 
  list(Hypothalamus = dplyr::filter(regions, region == 'Hypothalamus')$subregion,
       Midbrain = dplyr::filter(regions, region == 'Midbrain')$subregion,
       Hindbrain = dplyr::filter(regions, region == 'Hindbrain')$subregion)

# all genes
genes <- readr::read_csv("gene_biotypes.csv")

gene_list <- unique(genes$Gene)

biotypes <-
  genes %>%
  group_by(Biotype) %>%
  tally(.) %>%
  dplyr::arrange(desc(n)) %>%
  select(-n)

biotypes_list <- biotypes$Biotype
