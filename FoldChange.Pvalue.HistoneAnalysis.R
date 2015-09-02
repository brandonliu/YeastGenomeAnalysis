# 1. Download the ChromatinDB data
# 2. Use function for significant p-values and fold-changes to generate a subset of data
# 3. Generate boxplots for each histone marker and look at the distributions
# 4. Then compare those boxplots to boxplots of the subsets of significant p-value and fold-change data to see if there is significant variation

# Additional Notes:
# -----------------
# Function: Specify gene name, return back a list of significant genes
# Optional: Include p-value of interest, include positive fold change threshold, negative fold change threshold
load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/deleteome_raw_data.RData')

# <<<<<<<<< NEED TO MAKE THE FOLLOWING FUNCTION CALLS >>>>>>>>>

# 2. Need to call fixColnames to reformat the columns
# -------------------------------------------------
# colnames(deleteome$FC) = sapply(colnames(deleteome$FC), function(x) fixColnames(x))
# colnames(deleteome$p.value) = sapply(colnames(deleteome$p.value, function(x) fixColnames(x)))
# 
# 3. To fix the replicates:
# -------------------------
# make.unique(colnames(deleteome$p.value))
# make.unique(colnames(deleteome$FC))


# Convert all characters to upper case

colnames(deleteome$FC) = toupper(colnames(deleteome$FC))
#colnames(deleteome$p.value) = toupper(colnames(deleteome$p.value))
colnames(deleteome$p.value) = colnames(deleteome$FC)

# The deleteome data set includes a deleteome data object of two field (Fold Change (FC) and p-value (p.value))

library(org.Sc.sgd.db)

print(colnames(deleteome$FC))
print(colnames(deleteome$p.value))

get_significant_gene_interactions <- function(gene, p.value = 0.05, pos.FC = log2(1.7), neg.FC = log2(1/1.7)) {
  # First check if gene exists: user needs to input in capital letters
  if (!(gene %in%colnames(deleteome$FC)) | !(gene %in%colnames(deleteome$p.value))) stop("Gene not available in the data set. Please input gene names in uppercase format")
  print('Finding your significant gene interactions...')  
  
  # Generates a matrix of significant fold-changes based on a specified threshold
  sig_fold_changes = as.matrix(deleteome$FC[colnames(deleteome$FC) == gene])
  rownames(sig_fold_changes) = rownames(deleteome$FC)
  sig_fold_changes = subset(sig_fold_changes, sig_fold_changes >= pos.FC | sig_fold_changes <= neg.FC)
  
  # Generates a matrix of significant p-values according to a specified p-value threshold
  sig_p.values = as.matrix(deleteome$p.value[colnames(deleteome$p.value) == gene])
  rownames(sig_p.values) = rownames(deleteome$p.value)
  sig_p.values = subset(sig_p.values, sig_p.values <= p.value)
  
  
  #This optional phase returns a list of the two objects: significant fold changes, significant p-values
  # results = list("Significant Fold Changes" = sig_fold_changes, "Significant P-values" = sig_p.values)
  
  
  # The gene list returns the list of names of the genes with both these qualifications
  geneList = intersect(rownames(sig_p.values), rownames(sig_fold_changes))
  return(geneList)  
}






# Trim trailing white space in row names:
library(reshape2)
library(ggplot2)
library(IRanges)



# Or use load to load the data from somewhere else..... (.RData file)

data = read.csv(file = "/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/ChromatinDBPromoterData.csv", sep = "\t", row.names = 1)


trim.trailing <- function (x) sub("\\s+$", "", x)
rownames(data) = trim.trailing(rownames(data))

# All the row names have an extra space at the end of the name for some reason

# Function: compareHistoneEnrichment
# -------------------------------

# From the Holstege lab deleteome data (which is essentially a data set containing values for gene expression levels after single gene deletions), we take genes with significant fold change and p-value. We then take that subset of genes and create a boxplot of histone enrichment for the aggregate of those genes (histone data from ChromatinDB), then overlay that for the genome-wide averages for histone enrichment, thus creating a spread of boxplots for each histone marker.

compareHistoneEnrichment <- function(gene) {
  significantGenes = get_significant_gene_interactions(gene)

  histoneSubset = data[significantGenes,]

  print('Generating boxplots')
  
  ggplot(melt(data), aes(x = variable, y = value)) + geom_boxplot(aes(fill = variable),position = "identity", alpha = 0.5) + geom_boxplot(data = melt(histoneSubset), aes(x = variable, y = value, fill = "red"), position = "identity", main = gene) + theme(plot.title = element_text(size = 20, face = "bold", vjust = 1), axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Enrichment Score") + xlab("Histone") + ggtitle(paste(gene, "Histone Enrichment Comparison"))

}
