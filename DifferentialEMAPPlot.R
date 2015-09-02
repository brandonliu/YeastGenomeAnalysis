# Differential E-MAP plot

library(ggplot2)

load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/RelationalDatabaseProject/EMAPObjects.RData')

plotDifferential <- function(gene, Condition1, Condition2, ..., Subset = FALSE) {
  if (gene %in% colnames(Condition1)) {
    result = data.frame(Condition1[,1], Condition1[,gene], Condition2[,gene])
    names(result) <- c('Gene', 'SScore', 'StressScore')
    ggplot(result, aes(x = Gene, y = SScore, group = 1)) + geom_line() + geom_line(data = result, aes(x = Gene, y = StressScore, group = 1), color = "red") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else if (gene %in% Condition1[,1]) {
    result = data.frame(colnames(Condition1), Condition1[gene,], Condition2[gene,])
    names(result) <- c('Gene', 'SScore', 'StressScore')
    ggplot(result, aes(x = Gene, y = SScore, group = 1)) + geom_line() + geom_line(data = result, aes(x = Gene, y = StressScore, group = 1),color = "red") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
    stop (paste('Could not find '), gene)
  }
}