
#' Data frame of genes and their scaled degree in the network
#'
#' @name scaled_degree
#' @format A data frame with gene IDs in column 1 and scaled degree in column 2.
#' @examples 
#' data(scaled_degree)
#' @usage data(scaled_degree)
"scaled_degree"


#' Data frame of genes and modules
#' 
#' @name genes_modules
#' @format A data frame with gene IDs in column 1 and modules in column 2.
#' @examples 
#' data(genes_modules)
#' @usage data(genes_modules)
"genes_modules"


#' Hub genes
#' 
#' @name hubs
#' @format A data frame with hub gene IDs in column 1, and module and degree
#' in columns 2 and 3.
#' @examples 
#' data(hubs)
#' @usage data(hubs)
"hubs"


#' Network plotting data
#' 
#' Pre-computed data frame of network plotting data.
#' 
#' @name plotdata
#' @format A data frame with the x and y coordinates of nodes and edges between
#' them for plotting.
#' @examples 
#' data(plodata)
#' @usage data(plotdata)
"plotdata"


#' Module enrichment results
#' 
#' @name mod_enrich
#' @format A data frame with module enrichment results for MapMan bins, 
#' GO terms and Interpro domains, as well as associated adjusted P-values.
#' @examples 
#' data(mod_enrich)
#' @usage data(mod_enrich)
"mod_enrich"
