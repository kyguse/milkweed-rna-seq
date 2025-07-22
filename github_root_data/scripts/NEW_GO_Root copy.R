#####NEW GO Analysis on the Root data taking out the NA or absent values
###starting with the de_filtered_clean data.frame
go_data<-de_filtered_clean
go_data <- go_data %>%
  filter(go_terms_short != "." & go_terms_short != "" & !is.na(go_terms_short))

###Trying again
# Remove rows where go_terms_short is NA or "."
go_data <- go_data %>%
  filter(!is.na(go_terms_short) & go_terms_short != "NA")


# Check number of rows with valid GO terms
cat("Number of genes with GO annotations:", nrow(go_data), "\n")

# Load libraries
library(dplyr)
library(ggplot2)
library(forcats)

