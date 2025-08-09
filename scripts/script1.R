> species_og <- read.csv("data/species_OGlist.csv") #data of all genes and corresponding orthogroups
> specific_geneid <- readLines("data/specific_geneid.txt") #GeneID list of specifically expressed gene
> all_geneid <- c(species_og$GeneID)
> og_list <- split(species_og$GeneID, species_og$eggNOG_OGs)
> results_species <- data.frame()
  for (og in names(og_list)) {
  og_genes_species <- og_list[[og]]
  a <- sum(og_genes_species %in% specific_geneid)
  b <- sum(og_genes_species %in% setdiff(all_geneid, specific_geneid))
  c <- sum(setdiff(specific_geneid, og_genes_species) %in% all_geneid)
  d <- length(setdiff(all_geneid, union(og_genes_species, specific_geneid)))
  fisher_mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(fisher_mat, alternative = "greater") #Fisher's exact test
  results_species <- rbind(results_species, data.frame(
  OG = og,
  target_count = a,
  background_count = b,
  pvalue = ft$p.value
  ))
  }
> results_species$FDR <- p.adjust(results_species$pvalue, method = "BH")
> results_species_sig <- results_species[order(results_species$FDR), ] #sorting by FDR
> write.csv(results_species_sig, "output.csv")
