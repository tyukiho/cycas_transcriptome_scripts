> library(dplyr)
> library(tidyr)
> library(UpSetR)
> library(tibble)
> df_species1 <- read.csv("data/species1_OG_GO.csv") #Importing the GO term and orthogroup　dataset for each gene
  df_species2 <- read.csv("data/species2_OG_GO.csv")
> df_species1_long <- df_species1 %>% #Data formatting, GO is one per line
  pivot_longer(
  cols = starts_with("GO"),   
  names_to = "GO_col",        
  values_to = "GO"            
  ) %>%
  filter(!is.na(GO) & GO != "") 
  df_species2_long <- df_species2 %>%
  pivot_longer(
  cols = starts_with("GO"),   
  names_to = "GO_col",        
  values_to = "GO"            
  ) %>%
  filter(!is.na(GO) & GO != "") 
> id_species1 <- readLines("data/tissue-specific_geneid_species1.txt")  # tissue-specifically expressed gene ID
  id_species2 <- readLines("data/tissue-specific_geneid_species2.txt")
> go_species1 <- readLines("data/fertilization-cluster_GO_species1.txt")  # GO terms included in the fertilization cluster
  go_species2 <- readLines("data/fertilization-cluster_GO_species1.txt")
> specific_species1_id <- df_species1_long %>% #
  filter(GeneID %in% id_species1)
  specific_species1_id_go <- specific_species1_id %>%
  filter(GO %in% go_species1)
  species1_result <- specific_species1_id_go %>%
  select(GeneID, eggNOG_OGs) %>%
  distinct()
  #Extraction of orthogroups of genes corresponding to GO in the fertilization cluster among tissue-specific genes
> specific_species2_id <- df_species2_long %>%
  filter(GeneID %in% id_species2)
  specific_species2_id_go <- specific_species2_id %>%
  filter(GO %in% go_species2)
  species2_result <- specific_species2_id_go %>%
  select(GeneID, eggNOG_OGs) %>%
  distinct()
> OG_list_species1 <- unique(species1_result$eggNOG_OGs)
  OG_list_species2 <- unique(species2_result$eggNOG_OGs)
> OG_list_for_upset <- list(Species1 = OG_list_species1, Species2 = OG_list_species2)
> inputupset <- fromList(OG_list_for_upset) #Data formatting for UpSetR
  ogs <- unique(unlist(OG_list_for_upset))
  rownames(inputupset) <- ogs
> upset(inputupset, sets = c("Species2", "Species1") , order.by = "freq", keep.order = TRUE) #output of UpSet diagram 
> write.csv(inputupset, "upset_output.csv", row.names = TRUE) #output of UpSet diagram breakdown