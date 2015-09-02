
# Generates subunits based on an appropriate filename and given a set of gene names

subsetCorrelationMatrix <- function(filename, EMAPCondition, genes) {
#   data = read.csv(file = filename)
  load(filename)
#   rownames(data) = data[,1]
#   data = data[,-1]
#   print(dim(data))
  y = melt(EMAPCondition[genes, genes])
  plot(hclust(y))
#   filename = paste("NewSubunits",
#   write.csv(y, file = "Subunits_UntreatedQueryCor")
  
}