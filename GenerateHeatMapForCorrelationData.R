# Generate heatmap

library(reshape)
library(ggplot2)
library(ggdendro)

# Input a file name from one of the correlation matrices
createHeatMap <- function(filename) {
  data = read.csv(file = filename)
  rownames(data) = data[,1]
  y = melt(data)
  names(y) <- c("Gene1", "Gene2", "value")
  print(y)
  p <- ggplot(y, aes(x = Gene1, y = Gene2))
  p + geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "darkblue") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

  # Need to cluster them though
  
}

# Make a heatmap using the heatmap function
# data = read.csv(file = "/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/EMAP2015CorrelationFiles/UntreatedQueries.csv")
# result = produceDistMatrix(data)
# UntreatedQueryCor = result
# heatmap(as.matrix(UntreatedQueryCor))
