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


#----Hubs.rda-------------------------------------------------------------------
load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_nematode_hubs.rda")
hubs_nematode <- hubs
load("~/Dropbox/Working_from_home/GWAS_GCN_pest/products/result_files/shiny_hubs.rda")
hubs_insect <- hubs$Gene

hubs <- list(insect = hubs_insect, nematode = hubs_nematode)
usethis::use_data(hubs, compress="xz", overwrite=TRUE)
