# Task: get correlations of Arp5 occupancy and histone modifications on Arp5-regulated genes.  I believe Devin has two Arp5 occupancy datasets.  Would it be possible for you to do these correlations?  As a comparison it would also be great to do correlations with other INO80 subunits, like Ino80, Ies6, Arp8 etc.

# The levels of subunits present might be regulating gene expression?
# Actual numbers are proportional to how many reads mapped their from the sequencing

load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/Arp5_occu/MNase_ChIPexo_Arp5_Ino80.RData')
#contains two objects: 'MNase_Arp5_Ino80', 'ChIPexo_Arp5_Ino80'

source('~/Documents/Stanford/Data Analysis/RelationalDatabaseProject/FoldChange.Pvalue.HistoneAnalysis.MultiplePlots.R')


Histone_Data = read.csv(file = "/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/ChromatinDBPromoterData.csv", sep = "\t", row.names = 1)

# Function:
# Inputs:
# -----------------
# 1) Identify genes that are regulated (differentially expressed) by X input gene. (using get_significant_gene_interactions) --> list of gene names
# 2) Correlate histone modifications and X gene occupancy between different differentially expressed genes (e.g. take the set with histone modifications, take the set of gene occupancy numbers, correlate them? -- for each histone?)
# 3) Look at all INO80 subunits

# Function: Correlates histone modification with ARP5/INO80 occupancy based on a subset of genes that are differentially expressed when said gene is knocked out
# Prints out a data frame for all the correlations between values
NucOccup_Histone_Corr_Total <- function(OccupyingGene, DiffGenes = get_significant_gene_interactions('ARP5')$totalGeneList) {
  # Subset of genes for ARP5 deletion
  print(DiffGenes)
  mna_subset = MNase_Arp5_Ino80[DiffGenes, OccupyingGene]
  exo_subset = ChIPexo_Arp5_Ino80[DiffGenes, OccupyingGene]
  m = matrix(0, ncol = 2, nrow = 22)
  result = data.frame(m)
  names(result) = c('MNase Correlation', 'ChIPExo Correlation')
  result[,1] = as.character(result[,1])
  print('Processing...')
  for (i in 1: length(colnames(Histone_Data))) {
    hist_subset = Histone_Data[DiffGenes, colnames(Histone_Data)[i]]
    MNaseCor = cor(hist_subset, mna_subset, use = "complete.obs")
    ChIPExoCor = cor(hist_subset, exo_subset, use = "complete.obs")
    entry = data.frame(MNaseCor, ChIPExoCor)
    names(entry) = c('MNase Correlation', 'ChIPExo Correlation')
    result[i,] = entry
  }
  rownames(result) = colnames(Histone_Data)
  print(result)
  
  
}


# correlate total linear regression of nucleosome occupancy with the linear regression of histone modification

# specifically H2AZ?