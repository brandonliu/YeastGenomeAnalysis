# Performs wilcox test comparing histone enrichment genome-wide against genes with significant fold-change and p-values

# I generated two versions of the data, one with p-values adjusted and another with p-values not adjusted

source("/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/chromatin_remodeler_genes.R")


histoneData_wilcoxonTest <- function() {
  Chromatin_Remod = lapply(CHROMREMOD, function(x) {
    p.adjust(mapply(function(x, y) {
      wilcox.test(x, y)$p.value  
    }, data[get_significant_gene_interactions(x),], data, SIMPLIFY = TRUE))
  })
  names(Chromatin_Remod) = CHROMREMOD
  return (Chromatin_Remod)
}



# p adjust each of these?

# Sapply with the list of genes that we have been given
# append each list to a data frame


# Multiple hypothesis testing
# look up function: p.adjust

# Function to get the genome-wide data for 