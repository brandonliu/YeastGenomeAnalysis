# Generates correlations based on queries inputted into data file


load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/deleteome_raw_data.RData')

# <<<<<<<<< NEED TO MAKE THE FOLLOWING FUNCTION CALLS >>>>>>>>>
# Note: some of the data is not formatted well, so it is necessary to modify it
# 2. Need to call fixColnames to reformat the columns
# -------------------------------------------------
colnames(deleteome$FC) = sapply(colnames(deleteome$FC), function(x) fixColnames(x))
colnames(deleteome$p.value) = sapply(colnames(deleteome$p.value), function(x) fixColnames(x))
# 
# 3. To fix the replicates:
# -------------------------
make.unique(colnames(deleteome$p.value))
make.unique(colnames(deleteome$FC))


# Convert all characters to upper case, reformats the loaded data appropriately
colnames(deleteome$FC) = toupper(colnames(deleteome$FC))
colnames(deleteome$p.value) = colnames(deleteome$FC)

library(org.Sc.sgd.db)

# Converts systematic yeast names to standard names
convertToStandard <- function(systematicNames) {
  return(select(org.Sc.sgd.db,keys=systematicNames,keytype='ORF',columns='GENENAME')$GENENAME)
}

# Converts yeast standard names to systematic names
convertToSystematic <- function(standardNames) {
  return(select(org.Sc.sgd.db,keys=standardNames,keytype='GENENAME',columns='ORF')$ORF)
}

# Performs correlations based on data from the Holstege lab deleteome
HolstegeCorrelations <- function(firstGene, geneSet) {
  print(paste('Printing correlations against: ', firstGene))
  firstGeneName = convertToSystematic(standardNames = firstGene)
  for (i in 1:length(geneSet)) {
    geneName = convertToSystematic(standardNames = geneSet[i])
    value = cor(deleteome$FC[firstGeneName,], deleteome$FC[geneName,])
    print(paste('Correlation: ', firstGene, '-', geneSet[i], ' ', value))
  }  
}
