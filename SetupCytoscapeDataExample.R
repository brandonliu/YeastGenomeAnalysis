# Using the INO80_UntreatedQueryMatrix, which is a subset of the UntreatedQueryCor matrix that only includes INO80 subunit data
# There is a script that I can use to subset data from these files by specifying a vector of gene deletion names.


# Function takes in a square correlation matrix. Must have correct row and column names
cytoscapeTable <- function(matrix) {
  
  data = cbind(rownames(matrix), matrix)
  
  # Need to first make the first column the names that way you can properly melt the data
  
  names(data)[1] = "Gene"
  
  y = melt(data, id.vars = "Gene", measure.vars = names(data)[2:21])
  
  
  # Need to threshold the data by a certain amount otherwise cytoscape will just include all the edges
  
  temp = subset(y, abs(y[,3]) >= 0.25)
  temp = temp[order(temp[,1]),]
  
  return(temp)
  #write.csv(y, file = "MeltedINO80SubunitCorrelationData")
  
}