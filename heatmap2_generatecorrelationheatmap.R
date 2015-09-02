# Generate heat map using heatmap.2 function and choose your own color palette
# Automatically generates the dendrogram for you

library(gplots)
library(reshape2)
library(RColorBrewer)

#produceDistMatrix function
source('~/Documents/Stanford/Data Analysis/Test Scripts/DistMatrixPractice.R')

createHeatMap3 <- function(filename) {
  data = read.csv(file = filename)
  result<- produceDistMatrix(data)
  UntreatedQueryCor = result
  print(UntreatedQueryCor)
  # can add labels, etc.

  
  # Color breaks
  #quantile.range <- quantile(UntreatedQueryCor, probs = seq(0, 1, 0.01))
  #palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
  
  
  color.palette <- colorRampPalette(c("#3182bd",'#deebf7', '#ffeda0', '#feb24c',  '#f03b20'))#(length(palette.breaks) - 1)
#   par(mar = c(10, 2, 2, 10))
  heatmap.2(as.matrix(UntreatedQueryCor), col = color.palette, main = "Untreated Query Correlation Heatmap",, key.title = 'Legend', keysize = 1.0, offsetCol = -0.5, margins = c(9, 9))
  
}

# Function subsets the gene matrix based on the names that you input
# Takes in a vector

INO80_Subunits = c('ARP5_DELETION2','ARP5_DELETION2CLONE4','ARP5_DELETION3CLONE2','ARP5_DELETION3CLONE3','ARP6_DELETION','ARP8_DELETIONCLONE1','ARP8_DELETIONCLONE4','IES1_DELETION','IES2_DELETION','IES3_DELETION','IES4_DELETION','IES5_DELETION','IES6_DELETION1CLONE1','IES6_DELETION2CLONE2','IES6_DELETION3','IES6_DELETION4','IES6_DELETION6','IES6_K145|166A','NHP10_DELETION','RVB1_DAMP')

SWR1_genes = c('ARP6_DELETION', 'HTZ1_DAMP', 'HTZ1_DELETION', 'HTZ1_DELETION2BATCH', 'SWC3_DELETION', 'SWC5_DELETION','SWC7_DELETION','SWR1_DELETIONCLONE1', 'SWR1_DELETIONCLONE2', 'VPS71_DELETION', 'VPS72_DELETION', 'VPS72_DELETION2', 'VPS72_YL1C')

subsetCorrelationMatrix <- function(filename, genes) {
  data = read.csv(file = filename)
  rownames(data) = data[,1]
  data = data[,-1]
  print(dim(data))
  y = melt(data[genes, genes])
  write.csv(y, file = "SWR1Subunits_UntreatedQueryCor")
  
}