# Builds upon existing heatmap script (ggplot2) in order to generate heat map diagrams that come with the dendrogram format.

source('~/Documents/Stanford/Data Analysis/Test Scripts/DistMatrixPractice.R')

library(reshape2)
library(ggplot2)
library(ggdendro)
library(grid)
library(proxy)


#http://stackoverflow.com/questions/6673162/reproducing-lattice-dendrogram-graph-with-ggplot2

# data(mtcars)
# x <- as.matrix(scale(mtcars))
# dd.col <- as.dendrogram(hclust(dist(x)))
# col.ord <- order.dendrogram(dd.col)
# 
# dd.row <- as.dendrogram(hclust(dist(t(x))))
# row.ord <- order.dendrogram(dd.row)
# 
# xx <- scale(mtcars)[col.ord, row.ord]
# xx_names <- attr(xx, "dimnames")
# df <- as.data.frame(xx)
# colnames(df) <- xx_names[[2]]
# df$car <- xx_names[[1]]
# df$car <- with(df, factor(car, levels=car, ordered=TRUE))
# 
# mdf <- melt(df, id.vars="car")


# Input a file name from one of the correlation matrices
createHeatMap2 <- function(filename) {
  data = read.csv(file = filename)
  data[is.na(data)] <- 0
  rownames(data) = data[,1]
  data = data[,-1]
  dend.col <- as.dendrogram(hclust(dist(data)))
  dend.row <- as.dendrogram(hclust(dist(t(data))))
  col.order <- order.dendrogram(dend.col)
  row.order <- order.dendrogram(dend.row)
  
  data1 <- scale(data)[col.order, row.order]
  
  y = melt(data1)
  
  #To generate the corresponding dendrograms

  ddata_x <- dendro_data(dend.row)
  ddata_y <- dendro_data(dend.col)
  
  ### Set up a blank theme
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  names(y) <- c("Gene1", "Gene2", "Correlation")
  p1 <- ggplot(y, aes(x = Gene1, y = Gene2)) + geom_tile(aes(fill = Correlation)) + scale_fill_gradient2(low = "blue", high = "red") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('Untreated Correlation Heatmap') + theme(plot.title=element_text(face="bold", size = 36), legend.justification = c(0, 1), legend.position = c(1,1), plot.margin = unit(c(3, 3, 3, 3),"cm")) + theme(legend.background = element_rect(fill="grey93", size=0.5, linetype="solid", colour ="grey93"))
  plot(p1)
  
#   # Dendrogram 1
#   p2 <- ggplot(segment(ddata_x)) + 
#     geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
#     theme_none + theme(axis.title.x=element_blank())
#   
#   # Dendrogram 2
#   p3 <- ggplot(segment(ddata_y)) + 
#     geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
#     coord_flip() + theme_none
# 
#   grid.newpage()
#   print(p1, vp=viewport(0.8, 0.8, x=0.4, y=0.4))
#   print(p2, vp=viewport(0.52, 0.2, x=0.45, y=0.9))
#   print(p3, vp=viewport(0.2, 0.8, x=0.9, y=0.4))
  
  
  
  # Need to cluster them though
  
}