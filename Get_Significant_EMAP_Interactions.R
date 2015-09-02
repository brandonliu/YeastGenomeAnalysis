# Generate the top interactions for query / test gene inputs
# for a given EMAP condition.
# Introduces a cutoff of -2.5 for aggravating (negative) genetic interactions
# and a cutoff of 2.0 for alleviating (positive) genetic interactions


# Requires that the user have approrpiately formatted EMAP data
# for object EMAPCondition. You need to run GenerateDataStructuresForEMAPData.R
# to generate all the appropriate files first


# The user inputs a gene name and the program will return results for all entries containing that string.

load('/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/2015EMAP_Rap_EtOH_NaCl_DataSets/EMAPObjects.RData')

generateTopQueryInteractions <- function(..., EMAPCondition) {
  queries = list(...)
  numQueries = length(queries)
  for (i in 1:numQueries) {
    # Temp holds a row that corresponds to the current query of interest
    #temp = data.frame(EMAPCondition[EMAPCondition[,1] == queries[i],])
    temp = data.frame(EMAPCondition[grepl(pattern = queries[i], x = EMAPCondition[,1]),])
    
    if (nrow(temp) == 0) {
      stop('Could not find requested gene.')
    }
    print('We have found the following results...')
    print(as.character(temp[,1]))
    print('Showing you the first one now...')
    print('Processing....')
    
    # Transpose the horizontal matrix
    temp1 = t(temp[,-1])
    temp1 = na.omit(temp1)
      # Subsets the data to only include significant interactions
      newtemp = matrix(temp1[,1])
      rownames(newtemp) = rownames(temp1)
      temp2 = subset(newtemp, newtemp <= -2.5 | newtemp >= 2.0)
      
      # Reformat the data frame to include the gene name as a column
      temp2 = data.frame(rownames(temp2), temp2)
      names(temp2) = c(as.character(temp[1,1]), "SScore")
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
      temp = EMAPCondition[,-1]
      
      # Method essentially follows the compositon of the generateTopQueryInteractions method except we transpose the dataframe so that queries are on the x-axis and tests are on the y-axis
      temp = t(temp)
      
      # Subset the data frame to find tests that contain the search string
      temp = subset(temp,grepl(pattern = tests[i], x = as.matrix(rownames(temp))),)
      
      if (nrow(temp) == 0) {
        stop('Could not find requested gene.')
      }
      print('We have found the following results....')
      print(rownames(temp))
      print('Showing you the first one now.....')
      print('Processing...')

      colnames(temp) = (EMAPCondition[,1])
      temp1 = t(temp)

      # Remove NA values
      temp1 = temp1[!is.na(temp),]
                
      
      
      # Subset based on significance value
      SScore = subset(temp1,temp1 <= -2.5 | temp1 >= 2.0)
        
      temp3 = data.frame(SScore, row.names = names(SScore))
      print(temp3)
  }
}