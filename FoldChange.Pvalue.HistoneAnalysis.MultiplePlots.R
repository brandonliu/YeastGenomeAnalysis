
# This script performs correlation analysis for histone modification data and gene expression data. 
# Unlike the other version of this script, this one also performs additional analysis that plots
# multiple correlations side-by-side, so you can compare averages of histone modifications for 
# genes that are up/down regulated in the presence of a single gene deletion.


# INSTRUCTIONS FOR USE
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


convertToStandard <- function(systematicNames) {
  return(select(org.Sc.sgd.db,keys=systematicNames,keytype='ORF',columns='GENENAME')$GENENAME)
}

convertToSystematic <- function(standardNames) {
  return(select(org.Sc.sgd.db,keys=standardNames,keytype='GENENAME',columns='ORF')$ORF)
}

get_significant_gene_interactions <- function(gene, p.value = 0.05, pos.FC = log2(1.7), neg.FC = log2(1/1.7)) {
  # First check if gene exists: user needs to input in capital letters
  if (!(gene %in%colnames(deleteome$FC)) | !(gene %in%colnames(deleteome$p.value))) stop("Gene not available in the data set. Please input gene names in uppercase format")
  print('Finding your significant gene interactions...')  
  
  # Generates a matrix of significant fold-changes based on a specified threshold
  sig_fold_changes = as.matrix(deleteome$FC[colnames(deleteome$FC) == gene])
  rownames(sig_fold_changes) = rownames(deleteome$FC)
  sig_pos_fold_changes = subset(sig_fold_changes, sig_fold_changes >= pos.FC)
  sig_neg_fold_changes = subset(sig_fold_changes, sig_fold_changes <= neg.FC)
  total_sig_fold_changes = rbind(sig_pos_fold_changes, sig_neg_fold_changes)
  #sig_fold_changes = subset(sig_fold_changes, sig_fold_changes >= pos.FC | sig_fold_changes <= neg.FC)
  
  # Generates a matrix of significant p-values according to a specified p-value threshold
  sig_p.values = as.matrix(deleteome$p.value[colnames(deleteome$p.value) == gene])
  rownames(sig_p.values) = rownames(deleteome$p.value)
  sig_p.values = subset(sig_p.values, sig_p.values <= p.value)
  
  
  #This optional phase returns a list of the two objects: significant fold changes, significant p-values
  # results = list("Significant Fold Changes" = sig_fold_changes, "Significant P-values" = sig_p.values)
  
  
  # The gene list returns the list of names of the genes with both these qualifications
  totalGeneList = intersect(rownames(sig_p.values), rownames(total_sig_fold_changes))
  posFoldChange = intersect(rownames(sig_p.values), rownames(sig_pos_fold_changes))
  negFoldChange = intersect(rownames(sig_p.values), rownames(sig_neg_fold_changes))
  pValues = rownames(sig_p.values)
  Result = list(totalGeneList = totalGeneList, posFoldChange = posFoldChange, negFoldChange = negFoldChange, pValues = pValues)
  print("Returning results: 'totalGeneList', 'posFoldChange', 'negFoldChange', 'pValues")
  return(Result)
}






# Trim trailing white space in row names:
library(reshape2)
library(ggplot2)
library(IRanges)
library(plyr)
library(gridExtra)

# Or use load to load the data from somewhere else..... (.RData file)

data = read.csv(file = "/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/ChromatinDBPromoterData.csv", sep = "\t", row.names = 1)


trim.trailing <- function (x) sub("\\s+$", "", x)
rownames(data) = trim.trailing(rownames(data))


# All the row names have an extra space at the end of the name for some reason

# Function: compareHistoneEnrichment
# -------------------------------

# From the Holstege lab deleteome data (which is essentially a data set containing values for gene expression levels after single gene deletions), we take genes with significant fold change and p-value. We then take that subset of genes and create a boxplot of histone enrichment for the aggregate of those genes (histone data from ChromatinDB), then overlay that for the genome-wide averages for histone enrichment, thus creating a spread of boxplots for each histone marker.



compareHistoneEnrichment <- function(gene, ...) {
  ResultGenes = get_significant_gene_interactions(gene, ...)

  totalGenes = ResultGenes$totalGeneList
  posGenes = ResultGenes$posFoldChange
  negGenes = ResultGenes$negFoldChange
  total_histone_subset = data[totalGenes,]
  pos_histone_subset = data[posGenes,]
  neg_histone_subset = data[negGenes,]
  
  data$Subset = "genome"
  total_histone_subset$Subset = gene
  pos_histone_subset$Subset = gene
  neg_histone_subset$Subset = gene
  
  total_df = data.frame(rbind(data, total_histone_subset))
  pos_df = data.frame(rbind(data, pos_histone_subset))
  neg_df = data.frame(rbind(data, neg_histone_subset))
  
  
  print('Generating boxplots')
  print ('This will print side-by-side boxplots')
  total_plot <- ggplot(melt(total_df), aes(x = variable, y = value, fill = Subset)) + geom_boxplot() + ggtitle(paste("Total genes regulated by ", gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab('Histone') + ylab('Enrichment Score(log2)')
  pos_plot <- ggplot(melt(pos_df), aes(x = variable, y = value, fill = Subset)) + geom_boxplot() + ggtitle(paste("Genes up-regulated by ", gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab('Histone') + ylab('Enrichment Score(log2)')
  neg_plot <- ggplot(melt(neg_df), aes(x = variable, y = value, fill = Subset)) + geom_boxplot() + ggtitle(paste("Genes down-regulated by ", gene)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab('Histone') + ylab('Enrichment Score(log2)')
  grid.arrange(total_plot, pos_plot, neg_plot, nrow = 1, ncol = 3)
}

# Things to do:
# Change the framing so it isn't so hard to read
# Vertical x-axis labels
# Change the names of the x and y labels
# Fix the legend name instead of just "Subset"
# Print out names of genes that correspond to each graph
