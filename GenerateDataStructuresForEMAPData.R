library(reshape2)
# --------------------------------------------
# User inputs filename
filename = "/Users/BrandonLiu/Documents/Stanford/Data\ Analysis/2015_Rap_NaCl_EtOH\ CSV.csv"
# generateData <- function(filename) {

#}

temp = as.data.frame(read.csv(filename))
numDataFrames = ncol(temp) - 6 # 6 non value columns

# --------------------------------------------
# Subsets out duplicates, melts data before it is converted to wide format
# Horizontal : test names, Vertical : query names
# The Data frame is then assigned to an appropriately named variable
  # e.g. EtOH, Rapamycin, 
for (i in 1: numDataFrames) {
  temp1 = data.frame(temp[2], temp[3], temp[i + 6])
  temp2 = subset(temp1, duplicated(temp1) == FALSE) # Eliminates duplicates (basically the NAs that can't be read)
  temp3 = melt(temp2, id.vars = c("testName", "queryName"))
  temp4 <- dcast(temp3, queryName ~ testName, value.vars = "value")
  assign(toString(names(temp1)[3]), temp4)
}
