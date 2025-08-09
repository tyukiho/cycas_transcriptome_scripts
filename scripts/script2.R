> library(GO.db)
> library(AnnotationDbi)
> library(ggplot2)
> library(ggrepel)
> species1_go = readLines("data/topGOlist_species1.txt") #topGO lists for each species
  species2_go = readLines("data/topGOlist_species2.txt")
> all_go <- unique(c(species1_go, species2_go)) #unique GO list
> get_ancestors <- function(go_id) {
    res <- AnnotationDbi::mget(go_id, GOBPANCESTOR, ifnotfound = NA)[[1]]
    if (is.na(res[1])) return(go_id) else return(unique(c(go_id, res)))
  } #get ancestors GOs
> go_ancestors <- lapply(all_go, get_ancestors) #ancestral GO assignment
  names(go_ancestors_go) <- all_go
> jaccardSim <- function(a, b) length(intersect(a, b)) / length(union(a, b)) #jaccard similarity matrix
  n <- length(all_go)
  sim_mat <- matrix(0, n, n)
  rownames(sim_mat) <- colnames(sim_mat) <- all_go
  for (i in seq_len(n)) {
    for (j in i:n) {
        s <- jaccardSim(go_ancestors[[i]], go_ancestors[[j]])
        sim_mat[i, j] <- s
        sim_mat[j, i] <- s
    }
  }
> dist_mat <- 1 - sim_mat #conversion to distant matrix
  dist_mat[is.na(dist_mat)] <- 1
> mds_coords <- cmdscale(as.dist(dist_mat), k = 2) #two-dimentional mapping
> make_df <- function(gos, species) {
    df <- as.data.frame(mds_coords[gos, ])
    df$GO <- gos
    df$Species <- species
    return(df)
  }
  df_species1 <- make_df(species1_go, "species1")
  df_species2 <- make_df(species2_go, "species2")
  plot_go <- rbind(df_species1, df_species2)
> go_counts <- table(plot_go$GO) #detecting overlapping GO between species
  plot_go$GO_Freq <- as.integer(go_counts[plot_go$GO])
  plot_go2 <- plot_go %>%
    group_by(GO) %>%
    summarise(
        Combo = paste(sort(unique(Species)), collapse = "-"),
        GO_Freq = first(GO_Freq),
        V1 = mean(V1), 
        V2 = mean(V2)
    ) %>%
    ungroup()
> custom_colors <- c(
    "species1" = "grey80",
    "species2" = "grey80",
    "species1-species2" = "grey50",
    ) #example
> ggplot(plot_go2, aes(x = V1, y = V2)) +
    geom_point(aes(color = Combo, size = factor(GO_Freq)), alpha = 0.8) +
    scale_size_manual(values = c("1" = 2, "2" = 5)) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    labs(
        title = "GO semantic map",
        subtitle = "Grouped by species combinations (Combo)",
        x = "MDS1", y = "MDS2", color = "Shared by"
    )
