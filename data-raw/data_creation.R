## code to prepare datasets in data/
# All raw files were generated in the code associated with the paper
library(here)

#----mod_enrich.rda----
load(here("data", "shiny_enrichment.rda"))
enrichment_results$TermID
# Mapman: 1:29 | GO: 30:1383
enrichment_results$Category <- "GO"
enrichment_results$Category[1:29] <- "MapMan"
mod_enrich <- enrichment_results
mod_enrich <- mod_enrich[, c(1,5,7,8)]
usethis::use_data(mod_enrich, compress="xz", overwrite=TRUE)

#----scaled_degree.rda----
load(here("data", "shiny_degree.rda"))
degree_mod <- merge(degree, genes_modules, by.x="row.names", by.y="Genes")
degree_l <- split(degree_mod, degree_mod$Modules)
degree_l2 <- lapply(degree_l, function(x) {
    y <- cbind(x, Max=max(x$kWithin))
    y$scaled <- y$kWithin / y$Max
    y <- y[, c("Row.names", "scaled")]
    colnames(y) <- c("Gene", "Scaled")
    return(y)
})
scaled_degree <- Reduce(rbind, degree_l2)
scaled_degree$Scaled <- round(scaled_degree$Scaled, 2)
usethis::use_data(scaled_degree, compress="xz", overwrite=TRUE)

#----plotdata.rda----
load(here("data", "shiny_edgelists.rda"))

# Find optimal r cutoff
find_optimal_r <- function(edges, cutoff = seq(0.4, 0.9, by=0.1)) {
    list_edges <- replicate(length(cutoff), edges, simplify = FALSE)
    list_degree <- BiocParallel::bplapply(seq_along(list_edges), function(x) {
        y <- list_edges[[x]][list_edges[[x]]$Weight >= cutoff[x], ]
        graph <- igraph::graph_from_data_frame(y)
        deg <- igraph::degree(graph, mode="all")
        return(deg)
    }, BPPARAM = BiocParallel::SerialParam())
    sft <- 0.4
    try(
        sft <- unlist(lapply(list_degree, function(x) {
            return(WGCNA::scaleFreeFitIndex(x)$Rsquared.SFT)
        })),
        silent = TRUE
    )
    best_index <- which.max(sft)
    suppressWarnings(try(min_val <- min(sft[sft > 0.8]), silent=TRUE))
    if(min_val == "Inf") { min_val <- max(sft) }
    try( best_index <- which(sft == min_val), silent=TRUE)
    best <- cutoff[best_index]
    return(best)
}

optimal_r <- lapply(seq_along(edges), function(x) {
    message("Working on module ", names(edges)[x])
    y <- find_optimal_r(edges = edges[[x]])
    return(y)
})

r_df <- data.frame(
    Module = names(edges), optimal_r = unlist(optimal_r)
)

# Filter edge list
filtered_edges <- lapply(seq_along(edges), function(x) {
    y <- edges[[x]][edges[[x]]$Weight >= r_df[x, 2], ]
    return(y)
})
filtered_edges_df <- Reduce(rbind, lapply(seq_along(filtered_edges), function(x) {
    y <- cbind(filtered_edges[[x]], Module = r_df[x, ])
    return(y)
}))
filtered_edges_df <- filtered_edges_df[, c(1:4)]
names(filtered_edges_df) <- c("Gene1", "Gene2", "Weight", "Module")
edges <- filtered_edges_df

# Create pre-processed data for network visualization
load(here("data", "genes_modules.rda"))
load(here("data", "scaled_degree.rda"))

create_plot_data <- function(edges, module) {
    fedges <- edges[edges$Module == module, ]
    geneIDs <- unique(c(as.character(fedges[,1]), 
                        as.character(fedges[,2])))
    col <- "module"
    nod_at <- data.frame(Gene = geneIDs)
    nod_at <- merge(nod_at, scaled_degree)
    graph <- igraph::graph_from_data_frame(d = fedges[, 1:3],
                                           vertices = nod_at, 
                                           directed = FALSE)
    n <- ggnetwork::ggnetwork(graph, layout = igraph::with_kk(), arrow.gap = 0)
    n$Scaled10 <- n$Scaled * 10
    n$Scaled4 <- n$Scaled * 4
    return(n)
}

modules <- unique(genes_modules$Modules)
n_list <- lapply(modules, function(x) {
    y <- create_plot_data(edges, x)
    return(y)
})
n_list_df <- Reduce(rbind, n_list)
n_list_df <- merge(n_list_df, genes_modules, by.x="name", by.y="Genes")
plotdata <- n_list_df

# Check minimum r for each module
plotdata %>%
    group_by(Modules) %>%
    summarise(min = min(Weight, na.rm=TRUE)) %>%
    print(n=Inf)
    
usethis::use_data(plotdata, compress="xz", overwrite=TRUE)


