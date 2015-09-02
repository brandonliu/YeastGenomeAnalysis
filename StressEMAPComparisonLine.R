# Differential E-MAP plot

library(ggplot2)

load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/EMAPObjects.RData')

plotDifference <- function(gene, Condition1, Condition2, ..., Subset = FALSE) {
  if (gene %in% colnames(Condition1)) {
    Diff = Condition1[,gene] - Condition2[,gene]
    result = data.frame(rownames(Condition1),Diff, row.names= NULL)
    names(result) <- c('Gene', 'Difference')
    ggplot(result, aes(x = Gene, y = Difference, group = 1)) + geom_line() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else if (gene %in% rownames(Condition1)) {
    Diff = Condition1[gene,] - Condition2[gene,]
    dataNames = data.frame(t(colnames(Condition1)))
    result = data.frame(dataNames, as.data.frame(Diff), row.names = NULL)
    print(dim(result))
    names(result) <- c('Gene', 'Difference')
    ggplot(result, aes(x = Gene, y = Difference, group = 1)) + geom_line() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
    stop (paste('Could not find '), gene)
  }
}