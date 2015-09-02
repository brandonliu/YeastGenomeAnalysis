
library(org.Sc.sgd.db)

# Function: getCommonNames
# Takes a vector of system gene names and converts it to a vector of readable common names
getCommonNames <- function(systemNames) {
  x <- org.Sc.sgdGENENAME
  # Get the gene names that are mapped to an ORF identifier
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  xx = as.data.frame(x[mapped_genes])
  indeces = match(systemNames, xx$systematic_name)
  print('If no match was found in the org.Sc.sgd.db package, then the systematic name is printed')
  result = xx[indeces,]
  missingIndeces = which(is.na(result))
  result[is.na(result)] = systemNames[missingIndeces]
  return (result)
}