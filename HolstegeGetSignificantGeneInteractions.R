# Function: Specify gene name, return back a list of significant genes
# Optional: Include p-value of interest, include positive fold change threshold, negative fold change threshold

# 1. need to load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/deleteome_raw_data.RData')
# 2. Need to call fixColnames to reformat the columns
  # colnames(deleteome$FC) = sapply(colnames(deleteome$FC), function(x) fixColnames(x))
  # colnames(deleteome$p.value) = sapply(colnames(deleteome$p.value, function(x) fixColnames(x)))
# 
# 3. make.unique(colnames(deleteome$p.value)) to fix the replicates
  # make.unique(colnames(deleteome$FC))

# The deleteome data set includes a deleteome data object of two field (Fold Change (FC) and p-value (p.value))

library(org.Sc.sgd.db)

colnames(deleteome$FC) = toupper(colnames(deleteome$FC))
colnames(deleteome$p.value) = toupper(colnames(deleteome$p.value))


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

  
  #This optional phase returns a list of two objects: significant fold changes, significant p-values
#         results = list("Significant Fold Changes" = sig_fold_changes, "Significant P-values" = sig_p.values)

  geneList = intersect(rownames(sig_p.values), rownames(sig_fold_changes))
  return(geneList)  
}