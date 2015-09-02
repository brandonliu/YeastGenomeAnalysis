# Generate the top interactions for query / test gene inputs
# for a given EMAP condition.
# Introduces a cutoff of -2.5 for aggravating (negative) genetic interactions
# and a cutoff of 2.0 for alleviating (positive) genetic interactions

# Responses are written into CSV files of the format "Queryname.SignificantInteractions" 

# Requires that the user have approrpiately formatted EMAP data
# for object EMAPCondition. You need to run GenerateDataStructuresForEMAPData.R
# to generate all the appropriate files first

load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/2015EMAP_Rap_EtOH_NaCl_DataSets/EMAPObjects.RData')

generateTopQueryInteractions <- function(..., EMAPCondition) {
  queries = list(...)
  numQueries = length(queries)
  for (i in 1:numQueries) {
    # Temp holds a row that corresponds to the current query of interest
    temp = data.frame(EMAPCondition[EMAPCondition[,1] == queries[i],])
    print(temp)
    # Transpose the horizontal matrix
    temp1 = t(temp[,-1])
    temp1 = na.omit(temp1)
    
    print(temp1)
    
    # Subsets the data to only include significant interactions
    
    rownames(newtemp) = rownames(temp1)
      temp2 = subset(newtemp, newtemp <= -2.5 | newtemp >= 2.0)
      
      # Reformat the data frame to include the gene name as a column
      temp2 = data.frame(rownames(temp2), temp2)
      names(temp2) = c(as.character(queries[i]), "SScore")
      # Orders the data in descending order
      
      temp2 = temp2[order(temp2$SScore, decreasing = TRUE),]
      rownames(temp2) = NULL
      
      print(temp2)
    
    
  }
}

# Generates significant interactions for inputted test genes
# for a specified EMAP Stress condition
generateTopTestInteractions <- function(..., EMAPCondition) {
  tests = list(...)
  numTests = length(tests)
  for (i in 1:numTests) {
    #temp = data.frame(EMAPCondition[names(EMAPCondition) == tests[i]], row.names = EMAPCondition[,1])
    
    
    temp = data.frame(EMAPCondition[grepl(pattern = tests[i], x = names(EMAPCondition)),])
    
    if (nrow(temp) == 0) {
      stop('Could not find requested gene. Please type full gene name.')
    } 
    if(length(temp[,1]) > 1) {
      print('We have found the following results...')
      print(temp[,1])
      print('Consider rerunning function with a specific query')
    }
  
    for (j in 1:ncol(temp)) {
      # Subsets the data to only include significant interactions
      newtemp = as.(t(temp[j,]))
      rownames(newtemp) = colnames(temp)
      print(newtemp)
      newtemp = newtemp[,-1]
      newtemp = subset(newtemp, newtemp <= -2.5 | newtemp >= 2.0)
      print(dim(newtemp))
      # Reformat the data frame to include the gene name as a column
      temp3 = data.frame(rownames(newtemp), newtemp)
      names(temp3) = c(as.character(temp[j,1]), "SScore")
      
      # Orders the data in descending order
      
      #temp2 = temp2[order(temp2$SScore, decreasing = TRUE),]
      
      
      rownames(temp3) = NULL
      
    }
    
    
#     temp = na.omit(temp)
#     temp2= subset(temp, temp <= -2.5 | temp >= 2.0)
#     
#     # Generate a reformatted data frame to include a column for query names
#     temp2 = data.frame(rownames(temp2), temp2)
#     names(temp2) = c(tests[i], "SScore")  
#     
#     # Orders the data in descending order
#     temp2 = temp2[order(temp2$SScore, decreasing = TRUE),]
#     rownames(temp2) = NULL
#     print(temp2)    
  }   
}