# Generate a cross-correlation plot from two different E-MAPs (i.e. x-axis represents correlation scores for untreated and y-axis represents stress test for the other one).
# Query a single gene, and two conditions and it should return a differential plot.
# Print out a list of genes whose differential is above a certain threshhold?

library(ggplot2)
load('EMAPCorrelationObjects.RData')

cross_correlation_differential_plot <- function(gene, DifferentialX, DifferentialY, title = "Untreated vs Stress Condition", Labels = FALSE, threshold = 0.3) {
    if (dim(DifferentialX) != dim(DifferentialY)) {
      stop('Incorrect comparison .... Please choose tables of equal size.')
    }
    if (!(gene %in% names(DifferentialX)) | !(gene %in% names(DifferentialY))) {
      stop('Gene not found in one or more tables. Please include a full gene name.')
    }
    Correlations = data.frame(DifferentialX[,1], DifferentialX[,gene], DifferentialY[,gene])
    names(Correlations) <- c('Gene', 'Untreated', 'Stress_Condition')
    p <- ggplot(Correlations, aes(x = Untreated, y = Stress_Condition)) + geom_point() + ggtitle(paste(gene, " ", title))
    if (Labels == TRUE) {
      p <- p + geom_text(aes(label = Gene))
    }
    significant_genes = subset(Correlations, abs(Correlations[,2] - Correlations[,3]) > threshold)
    plot(p)
 print(significant_genes[,1])
}
