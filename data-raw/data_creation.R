## code to prepare datasets in data/
# All raw files were generated in the code associated with the paper
library(here)

#----mod_enrich.rda----
load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_insect_enrichment.rda")
enrichment_insect <- enrichment_results
enrichment_insect$taxon <- "insect"

load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_nematode_enrichment.rda")
enrichment_nematode <- enrichment_results
enrichment_nematode$taxon <- "nematode"

mod_enrich <- rbind(enrichment_insect, enrichment_nematode)
usethis::use_data(mod_enrich, compress="xz", overwrite=TRUE)

#----genes_modules.rda----------------------------------------------------------
load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_gmodules.rda")
genes_modules_insect <- genes_modules
genes_modules_insect$taxon <- "insect"

load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_nematode_gmodules.rda")
genes_modules_nematode <- genes_modules
genes_modules_nematode$taxon <- "nematode"

genes_modules <- rbind(genes_modules_insect, genes_modules_nematode)
usethis::use_data(genes_modules, compress = "xz", overwrite = TRUE)

#----scaled_degree.rda----
scale_degree_by_module <- function(degree_df = NULL, genes_modules = NULL) {
    degree_mod <- merge(degree_df, genes_modules, 
                        by.x = "row.names", by.y="Genes")
    
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
    return(scaled_degree)
}

load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_degree.rda")
scaled_degree_insect <- scale_degree_by_module(
    degree, genes_modules[genes_modules$taxon == "insect", ]
)
scaled_degree_insect$taxon <- "insect"

load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_nematode_degree.rda")
scaled_degree_nematode <- scale_degree_by_module(
    degree, genes_modules[genes_modules$taxon == "nematode", ]
)
scaled_degree_nematode$taxon <- "nematode"

scaled_degree <- rbind(scaled_degree_insect, scaled_degree_nematode)
usethis::use_data(scaled_degree, compress="xz", overwrite=TRUE)


#----plotdata.rda---------------------------------------------------------------
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

create_plot_data <- function(edges, module, scaled_degree) {
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

edgelist2plotdata <- function(edges = NULL, genes_modules = NULL,
                              scaled_degree = NULL) {
    if(methods::is(edges, "list")) {
        optimal_r <- lapply(seq_along(edges), function(x) {
            message("Working on module ", names(edges)[x])
            y <- find_optimal_r(edges = edges[[x]])
            return(y)
        })
        
        r_df <- data.frame(
            Module = names(edges), optimal_r = unlist(optimal_r)
        )
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
    } else {
        filtered_edges_df <- edges
    }

    modules <- unique(genes_modules$Modules)
    n_list <- lapply(modules, function(x) {
        y <- create_plot_data(filtered_edges_df, x, scaled_degree)
        return(y)
    })
    n_list_df <- Reduce(rbind, n_list)
    plotdata <- merge(n_list_df, genes_modules, by.x="name", by.y="Genes")
    return(plotdata)
}


load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_edgelists.rda")
load(here("data", "genes_modules.rda"))
load(here("data", "scaled_degree.rda"))
genes_modules_insect <- genes_modules[genes_modules$taxon == "insect", 1:2]
scaled_degree_insect <- scaled_degree[scaled_degree$taxon == "insect", 1:2]
plotdata_insect <- edgelist2plotdata(
    edges, genes_modules_insect, scaled_degree_insect
)

load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_nematode_edgelists.rda")

genes_modules_nematode <- genes_modules[genes_modules$taxon == "nematode", 1:2]
scaled_degree_nematode <- scaled_degree[scaled_degree$taxon == "nematode", 1:2]
edges <- merge(edgelist, genes_modules_nematode, by.x="Node1", by.y="Genes")
rm(edgelist)
names(edges) <- c("Gene1", "Gene2", "Weight", "Modules")
edges <- merge(edges, genes_modules_nematode, by.x = "Gene2", by.y="Genes")
edges <- edges[edges$Modules.x == edges$Modules.y, ]
edges <- edges[, 1:4]
names(edges)[4] <- "Modules"
edges <- split(edges, edges$Modules)
edges2 <- lapply(edges, function(x) {
    y <- x
    y$Modules <- NULL
    return(y)
})
edges <- edges2
rm(edges2)
    
plotdata_nematode <- edgelist2plotdata(
    edges, genes_modules_nematode, scaled_degree_nematode
)

plotdata <- list(insect = plotdata_insect,
                 nematode = plotdata_nematode)
plotdata_reduced <- lapply(plotdata, function(x) {
    y <- x[x$Weight > 0.8, ]
    y$Weight <- NULL
    return(y)
})
usethis::use_data(plotdata_reduced, compress="xz", overwrite=TRUE) # then renamed to plotdata.rda

# This was saved and then moved to a different directory for size issues
usethis::use_data(plotdata, compress="xz", overwrite=TRUE)



#----Hubs.rda-------------------------------------------------------------------
load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_nematode_hubs.rda")
hubs_nematode <- hubs
load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_hubs.rda")
hubs_insect <- hubs$Gene

hubs <- list(insect = hubs_insect, nematode = hubs_nematode)
usethis::use_data(hubs, compress="xz", overwrite=TRUE)

