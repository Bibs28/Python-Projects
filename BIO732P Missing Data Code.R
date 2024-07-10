## Missing Data Challenge
## A. Moreno, F. Carocci, N. Malik, U. Afshar, 2024.

#Working directory to access fo reproducibility
mypath <- "/Users/francescacarocci/Desktop/BIO732P/"

#Verify that the path is set correctly
print(mypath)

#Install packages
install.packages("naniar")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("caret")
install.packages("VIM")
install.packages("tibble")
install.packages("impute")
install.packages("visdat")
install.packages("missForest")
install.packages("BiocManager")
install.packages("ggprofiler2")
install.packages("methods")
install.packages("tools")
install.packages("ggvenn")
install.packages("RColorBrewer")
install.packages("gridExtra")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
BiocManager::install("DESeq2", force = TRUE, dependencies = TRUE)
BiocManager::install('EnhancedVolcano') #package for pretty volcano plots

#Load packages
library(naniar)
library(ggplot2)
library(caret)
library(VIM)
library(dplyr)
library(tibble)
library(visdat)
library(missForest)
library(BiocManager)
library(impute)
library(DESeq2)
library(gprofiler2)
library(EnhancedVolcano)
library(ggvenn)
library(RColorBrewer)
library(gridExtra)

#Construct the full file path using paste0
countdata5 <- paste0(mypath, "countdata5.txt")
countdata10 <- paste0(mypath, "countdata10.txt")
countdata30 <- paste0(mypath, "countdata30.txt")

#Load datasets using the full file path
countdata5 <- read.table(countdata5, header = TRUE, row.names = 1)
countdata10 <- read.table(countdata10, header = TRUE, row.names = 1)
countdata30 <- read.table(countdata30, header = TRUE, row.names = 1)

#Check datasets dimensions
dim(countdata5)
dim(countdata10)
dim(countdata30)

#Transpose datasets and convert to data frames
countdata5 <- as.data.frame(t(countdata5))
countdata10 <- as.data.frame(t(countdata10))
countdata30 <- as.data.frame(t(countdata30))

#Check datasets dimensions
dim(countdata5)
dim(countdata10)
dim(countdata30)

#View datasets
View(countdata5)
View(countdata10)
View(countdata30)

#Function to remove columns with only zero values
remove_zero_columns <- function(data) {
  
  #Identify all the rows where all the values are zero
  non_zero_columns <- colSums(data != 0, na.rm = TRUE) > 0
  
  #Subset the dataset to only keep columns with non-zero values
  cleaned_data <- data[, non_zero_columns]
  
  #Remove columns with only zeros
  num_removed <- sum(!non_zero_columns)
  print(paste("Number of columns removed with only zero values:", num_removed))
  
  return(cleaned_data)
}

#Remove rows with only zero values from datasets
countdata5_clean <- remove_zero_columns(countdata5)
countdata10_clean <- remove_zero_columns(countdata10)
countdata30_clean <- remove_zero_columns(countdata30)

#Check datasets dimensions
dim(countdata5_clean)
dim(countdata10_clean)
dim(countdata30_clean)

#Visualisations of missing data
#Countdata5
countdata5viss <- vis_miss(countdata5_clean[, 1:50], warn_large_data = FALSE) +
  ggtitle("Missingness Visualisation of Data in countdata5 (First 50 genes)") +  #Add title to the plot
  labs(y = "Mice", x = "Genes") +  #Label the axes
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  #Rotate x-axis text for better readability
    plot.title = element_text(size = 18),  #Set the size of the plot title
    axis.title.x = element_text(size = 16),  #Set the size of the x-axis title
    axis.title.y = element_text(size = 16),  #Set the size of the y-axis title
    legend.text = element_text(size = 12),  #Set the size of the legend text
    legend.title = element_text(size = 14)  #Set the size of the legend title
  )

#Print the visualisation to the console
print(countdata5viss)

#Construct the full file path using paste0
output_path_5 <- paste0(mypath, "countdata5_missingness.png")

#Save the visualisation as a high-quality image
ggsave(filename = output_path_5, plot = countdata5viss, device = "png", dpi = 300, width = 10, height = 8)

#Countdata10
countdata10viss <- vis_miss(countdata10_clean[, 1:50], warn_large_data = FALSE) +
  ggtitle("Missingness Visualisation of Data in countdata10 (First 50 genes)") +  #Add title to the plot
  labs(y = "Mice", x = "Genes") +  #Label the axes
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  #Rotate x-axis text for better readability
    plot.title = element_text(size = 18),  #Set the size of the plot title
    axis.title.x = element_text(size = 16),  #Set the size of the x-axis title
    axis.title.y = element_text(size = 16),  #Set the size of the y-axis title
    legend.text = element_text(size = 12),  #Set the size of the legend text
    legend.title = element_text(size = 14)  #Set the size of the legend title
  )

#Print the visualisation to the console
print(countdata10viss)

#Construct the full file path using paste0
output_path_10 <- paste0(mypath, "countdata10_missingness.png")

#Save the visualisation as a high-quality image
ggsave(filename = output_path_10, plot = countdata10viss, device = "png", dpi = 300, width = 10, height = 8)

#Countdata30
countdata30viss <- vis_miss(countdata30_clean[, 1:50], warn_large_data = FALSE) +
  ggtitle("Missingness Visualisation of Data in countdata30 (First 50 genes)") +  #Add title to the plot
  labs(y = "Mice", x = "Genes") +  #Label the axes
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  #Rotate x-axis text for better readability
    plot.title = element_text(size = 18),  #Set the size of the plot title
    axis.title.x = element_text(size = 16),  #Set the size of the x-axis title
    axis.title.y = element_text(size = 16),  #Set the size of the y-axis title
    legend.text = element_text(size = 12),  #Set the size of the legend text
    legend.title = element_text(size = 14)  #Set the size of the legend title
  )

#Print the visualisation to the console
print(countdata30viss)

#Construct the full file path using paste0
output_path_30 <- paste0(mypath, "countdata30_missingness.png")

#Save the visualisation as a high-quality image
ggsave(filename = output_path_30, plot = countdata30viss, device = "png", dpi = 300, width = 10, height = 8)

#Imputations functions
#Mean imputation
mean_imputation <- function(data) {
  data %>% mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
}

#Hot-deck imputation
hot_deck_imputation <- function(data) {
  
  #Create a copy of the data to store imputed values
  imputed_data <- data
  
  #Iterate over each column in the data frame
  for (j in 1:ncol(data)) {
    
    #Identify the missing values in the current column
    missing <- is.na(data[, j])
    
    #Replace missing values with the mean of the non-missing values in the current column
    imputed_data[missing, j] <- mean(data[!missing, j], na.rm = TRUE)
  }
  
  return(imputed_data)
}

#KNN imputation
knn_imputation <- function(data, k = 10) {
  
  #Convert data to matrix
  data_matrix <- as.matrix(data)
  
  #Imputation
  knn_imputed <- impute.knn(data_matrix, k = k)
  
  #Imputed data
  imputed_data <- knn_imputed$data
  
  #Convert matrix to dataframe
  imputed_df <- as.data.frame(imputed_data)
  
  return(imputed_df)
}

#Random Forest imputation
RF_imputation <- function(data, subsets = 1000, maxiter = 5, ntree = 50) {
  
  #Set the seed for reproducibility
  set.seed(123)
  
  #Determine the number of subsets
  num_subsets <- ceiling(ncol(data) / subsets)
  
  #Initialise a list to store imputed subsets
  imputed_subsets <- list()
  
  #Impute each subset of variables
  for (i in 1:num_subsets) {
    cat("Running imputation on subset", i, "of", num_subsets, "\n")
    
    #Determine the columns for the current subset
    start_col <- (i - 1) * subsets + 1
    end_col <- min(i * subsets, ncol(data))
    data_subset <- data[, start_col:end_col]
    
    #Perform imputation on the subset of variables
    imputed_subset <- missForest(data_subset, maxiter = maxiter, ntree = ntree)$ximp
    imputed_subsets[[i]] <- imputed_subset
  }
  
  #Combine all imputed subsets into one dataset
  combined_imputed_data <- do.call(cbind, imputed_subsets)
  
  #Return the combined imputed data
  return(combined_imputed_data)
}

#Countdata5
#Mean imputation
countdata5_imputedMean <- mean_imputation(countdata5_clean)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata5_imputedMean))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after mean imputation:")
  print(NA_counts)
} else {
  print("No NAs found after mean imputation.")
}

#Transpose the imputed dataset
countdata5_imputedMean_transposed <- t(countdata5_imputedMean)

#Construct the full file path using paste0
output_file <- paste0(mypath, "countdata5_imputedMean.txt")

#Save transposed data as a text file
write.table(countdata5_imputedMean_transposed, file = output_file, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Hot-deck imputation
countdata5_imputedHotDeck <- hot_deck_imputation(countdata5_clean)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata5_imputedHotDeck))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after hot-deck imputation:")
  print(NA_counts)
} else {
  print("No NAs found after hot-deck imputation.")
}

#Transpose the imputed dataset
countdata5_imputedHotDeck_transposed <- t(countdata5_imputedHotDeck)

#Construct the full file path using paste0
output_file_hotdeck <- paste0(mypath, "countdata5_imputedHotDeck.txt")

#Save transposed data as a text file
write.table(countdata5_imputedHotDeck_transposed, file = output_file_hotdeck, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#KNN imputation
countdata5_imputedKNN <- knn_imputation(countdata5_clean, k = 10)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata5_imputedKNN))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after KNN imputation:")
  print(NA_counts)
} else {
  print("No NAs found after KNN imputation.")
}

#Transpose the imputed dataset
countdata5_imputedKNN_transposed <- t(countdata5_imputedKNN)

#Construct the full file path using paste0
output_file_knn <- paste0(mypath, "countdata5_imputedKNN.txt")

#Save transposed data as a text file
write.table(countdata5_imputedKNN_transposed, file = output_file_knn, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Random Forest imputation
countdata5_imputedRF <- RF_imputation(countdata5_clean, subsets = 1000)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata5_imputedRF))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after Random Forest imputation:")
  print(NA_counts)
} else {
  print("No NAs found after Random Forest imputation.")
}

#Transpose datasets and convert to data frames
countdata5_imputedRF_transposed <- as.data.frame(t(countdata5_imputedRF))

#Construct the full file path using paste0
output_file_rf <- paste0(mypath, "countdata5_imputedRF.txt")

#Save transposed data as a text file
write.table(countdata5_imputedRF_transposed, file = output_file_rf, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Countdata10
#Mean imputation
countdata10_imputedMean <- mean_imputation(countdata10_clean)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata10_imputedMean))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after mean imputation:")
  print(NA_counts)
} else {
  print("No NAs found after mean imputation.")
}

#Transpose the imputed dataset
countdata10_imputedMean_transposed <- t(countdata10_imputedMean)

#Construct the full file path using paste0
output_file_mean <- paste0(mypath, "countdata10_imputedMean.txt")

#Save transposed data as a text file
write.table(countdata10_imputedMean_transposed, file = output_file_mean, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Hot-deck imputation
countdata10_imputedHotDeck <- hot_deck_imputation(countdata10_clean)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata10_imputedHotDeck))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after hot-deck imputation:")
  print(NA_counts)
} else {
  print("No NAs found after hot-deck imputation.")
}

#Transpose the imputed dataset
countdata10_imputedHotDeck_transposed <- t(countdata10_imputedHotDeck)

#Construct the full file path using paste0
output_file_hotdeck <- paste0(mypath, "countdata10_imputedHotDeck.txt")

#Save transposed data as a text file
write.table(countdata10_imputedHotDeck_transposed, file = output_file_hotdeck, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#KNN imputation
countdata10_imputedKNN <- knn_imputation(countdata10_clean, k = 10)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata10_imputedKNN))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after KNN imputation:")
  print(NA_counts)
} else {
  print("No NAs found after KNN imputation.")
}

#Transpose the imputed dataset
countdata10_imputedKNN_transposed <- t(countdata10_imputedKNN)

#Construct the full file path using paste0
output_file_knn <- paste0(mypath, "countdata10_imputedKNN.txt")

#Save transposed data as a text file
write.table(countdata10_imputedKNN_transposed, file = output_file_knn, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Random Forest imputation
countdata10_imputedRF <- RF_imputation(countdata10_clean, subsets = 1000)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata10_imputedRF))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after Random Forest imputation:")
  print(NA_counts)
} else {
  print("No NAs found after Random Forest imputation.")
}

#Transpose datasets and convert to data frames
countdata10_imputedRF_transposed <- as.data.frame(t(countdata10_imputedRF))

#Construct the full file path using paste0
output_file_rf <- paste0(mypath, "countdata10_imputedRF.txt")

#Save transposed data as a text file
write.table(countdata10_imputedRF_transposed, file = output_file_rf, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Countdata30
#Mean imputation
countdata30_imputedMean <- mean_imputation(countdata30_clean)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata30_imputedMean))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after mean imputation:")
  print(NA_counts)
} else {
  print("No NAs found after mean imputation.")
}

#Transpose the imputed dataset
countdata30_imputedMean_transposed <- t(countdata30_imputedMean)

#Construct the full file path using paste0
output_file_mean <- paste0(mypath, "countdata30_imputedMean.txt")

#Save transposed data as a text file
write.table(countdata30_imputedMean_transposed, file = output_file_mean, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Hot-deck imputation
countdata30_imputedHotDeck <- hot_deck_imputation(countdata30_clean)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata30_imputedHotDeck))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after hot-deck imputation:")
  print(NA_counts)
} else {
  print("No NAs found after hot-deck imputation.")
}

#Transpose the imputed dataset
countdata30_imputedHotDeck_transposed <- t(countdata30_imputedHotDeck)

# Construct the full file path using paste0
output_file_hotdeck <- paste0(mypath, "countdata30_imputedHotDeck.txt")

# Save transposed data as a text file
write.table(countdata30_imputedHotDeck_transposed, file = output_file_hotdeck, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#KNN imputation
countdata30_imputedKNN <- knn_imputation(countdata30_clean, k = 10)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata30_imputedKNN))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after KNN imputation:")
  print(NA_counts)
} else {
  print("No NAs found after KNN imputation.")
}

#Transpose the imputed dataset
countdata30_imputedKNN_transposed <- t(countdata30_imputedKNN)

#Construct the full file path using paste0
output_file_knn <- paste0(mypath, "countdata30_imputedKNN.txt")

#Save transposed data as a text file
write.table(countdata30_imputedKNN_transposed, file = output_file_knn, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Random Forest imputation
countdata30_imputedRF <- RF_imputation(countdata30_clean, subsets = 1000)

#Check there are no NAs after imputation
NA_counts <- colSums(is.na(countdata30_imputedRF))

#Print the number of NAs if any are present
if (any(NA_counts > 0)) {
  print("Number of NAs after Random Forest imputation:")
  print(NA_counts)
} else {
  print("No NAs found after Random Forest imputation.")
}

#Transpose datasets and convert to data frames
countdata30_imputedRF_transposed <- as.data.frame(t(countdata30_imputedRF))

#Construct the full file path using paste0
output_file_rf <- paste0(mypath, "countdata30_imputedRF.txt")

#Save transposed data as a text file
write.table(countdata30_imputedRF_transposed, file = output_file_rf, 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)

#Define Mean Absolute Error and Root Mean Square Error as numeric value functions
mae <- function(true_values, predicted_values) {
  
  #Calculate the absolute differences between true and predicted values
  absolute_errors <- abs(true_values - predicted_values)
  
  #Calculate the mean of the absolute errors
  mean_absolute_error <- mean(absolute_errors)
  return(mean_absolute_error)
}

rmse <- function(true_values, predicted_values) {
  
  #Calculate the squared differences between true and predicted values
  squared_errors <- (true_values - predicted_values)^2
  
  #Calculate the mean of the squared errors
  mean_squared_error <- mean(squared_errors)
  
  #Calculate the square root of the mean squared error
  root_mean_square_error <- sqrt(mean_squared_error)
  return(root_mean_square_error)
}

#Evaluation function of Mean imputation
evaluate_imputationMean <- function(dataset_path, na_fraction = 0.3, seed_split = 123, seed_na = 456) {
  
  #Load the dataset from the specified file path
  dataset <- read.table(dataset_path, header = TRUE, row.names = 1)
  
  #Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed_split) #Set seed for reproducibility of train-test split
  total_samples <- nrow(dataset) #Get the total number of samples
  train_size <- floor(0.7 * total_samples) #Calculate the number of training samples (70% of total)
  
  train_indices <- sample(seq_len(total_samples), size = train_size) #Randomly sample indices for training set
  train_data <- dataset[train_indices, ] #Create training set
  test_data <- dataset[-train_indices, ] #Create testing set with the remaining data
  
  #Introduce artificial NAs randomly in the testing set
  set.seed(seed_na) #Set seed for reproducibility of NA introduction
  test_data_with_na <- test_data #Create a copy of the testing set to introduce NAs
  num_values <- prod(dim(test_data_with_na)) #Calculate the total number of values in the testing set
  num_na <- round(num_values * na_fraction) #Calculate the number of NAs to introduce based on the specified fraction
  na_indices <- arrayInd(sample(num_values, num_na), dim(test_data_with_na)) #Randomly select indices to introduce NAs
  test_data_with_na[na_indices] <- NA #Introduce NAs at the selected indices
  
  #Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.") #Stop execution if no NAs are introduced
  }
  
  #Re-impute the missing values in the testing set using mean imputation
  test_data_reimputed <- test_data_with_na %>%
    mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) #Impute missing values with column means
  
  #Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.") #Stop execution if any NAs are still present after imputation
  }
  
  #Extract true and imputed values for comparison
  true_values <- test_data[na_indices] #Extract true values from the original testing set
  imputed_values <- test_data_reimputed[na_indices] #Extract imputed values from the re-imputed testing set
  
  #Ensure there are no NAs in the extracted true and imputed values
  if (any(is.na(true_values)) || any(is.na(imputed_values))) {
    stop("NAs present in the true or imputed values for comparison.") #Stop execution if any NAs are present in the extracted values
  }
  
  #Calculate MAE and RMSE between the imputed values and the true values in the test set
  mae_value <- mae(true_values, imputed_values) #Calculate Mean Absolute Error
  rmse_value <- rmse(true_values, imputed_values) #Calculate Root Mean Squared Error
  
  #Return the results as a list
  return(list(MAE = mae_value, RMSE = rmse_value))
}

#Apply evaluation function of Mean imputation
countdata5_imputedMean_evaluation <- evaluate_imputationMean(paste0(mypath, "countdata5_imputedMean.txt"))
cat("Mean Absolute Error (MAE):", countdata5_imputedMean_evaluation$MAE, "\n")
cat("Root Mean Squared Error (RMSE):", countdata5_imputedMean_evaluation$RMSE, "\n")

countdata10_imputedMean_evaluation <- evaluate_imputationMean(paste0(mypath,"countdata10_imputedMean.txt"))
cat("Mean Absolute Error (MAE):", countdata10_imputedMean_evaluation$MAE, "\n")
cat("Root Mean Squared Error (RMSE):", countdata10_imputedMean_evaluation$RMSE, "\n")

countdata30_imputedMean_evaluation <- evaluate_imputationMean(paste0(mypath,"countdata30_imputedMean.txt"))
cat("Mean Absolute Error (MAE):", countdata30_imputedMean_evaluation$MAE, "\n")
cat("Root Mean Squared Error (RMSE):", countdata30_imputedMean_evaluation$RMSE, "\n")


#Evaluation function of Hot-Deck imputation
evaluate_imputationHotDeck <- function(dataset_path, na_fraction = 0.3, seed_split = 123, seed_na = 456) {
  
  #Load the dataset from the specified file path
  dataset <- read.table(dataset_path, header = TRUE, row.names = 1)
  
  #Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed_split) #Set seed for reproducibility of train-test split
  total_samples <- nrow(dataset) #Get the total number of samples
  train_size <- floor(0.7 * total_samples) #Calculate the number of training samples (70% of total)
  
  train_indices <- sample(seq_len(total_samples), size = train_size) #Randomly sample indices for training set
  train_data <- dataset[train_indices, ] #Create training set
  test_data <- dataset[-train_indices, ] #Create testing set with the remaining data
  
  #Introduce artificial NAs randomly in the testing set
  set.seed(seed_na) #Set seed for reproducibility of NA introduction
  test_data_with_na <- test_data #Create a copy of the testing set to introduce NAs
  num_values <- prod(dim(test_data_with_na)) #Calculate the total number of values in the testing set
  num_na <- round(num_values * na_fraction) #Calculate the number of NAs to introduce based on the specified fraction
  na_indices <- arrayInd(sample(num_values, num_na), dim(test_data_with_na)) #Randomly select indices to introduce NAs
  test_data_with_na[na_indices] <- NA #Introduce NAs at the selected indices
  
  #Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.") #Stop execution if no NAs are introduced
  }
  
  #Re-impute the missing values in the testing set using hot-deck imputation
  test_data_reimputed <- test_data_with_na #Create a copy of the testing set with NAs for re-imputation
  for (j in 1:ncol(test_data_with_na)) {
    missing <- is.na(test_data_with_na[, j]) #Identify missing values in the current column
    test_data_reimputed[missing, j] <- mean(test_data_with_na[!missing, j], na.rm = TRUE) #Impute missing values with column means
  }
  
  #Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.") #Stop execution if any NAs are still present after re-imputation
  }
  
  #Extract true and imputed values for comparison
  true_values <- test_data[na_indices] #Extract true values from the original testing set
  imputed_values <- test_data_reimputed[na_indices] #Extract imputed values from the re-imputed testing set
  
  #Ensure there are no NAs in the extracted true and imputed values
  if (any(is.na(true_values)) || any(is.na(imputed_values))) {
    stop("NAs present in the true or imputed values for comparison.") #Stop execution if any NAs are present in the extracted values
  }
  
  #Calculate MAE and RMSE between the imputed values and the true values in the test set
  mae_value <- mae(true_values, imputed_values) #Calculate Mean Absolute Error
  rmse_value <- rmse(true_values, imputed_values) #Calculate Root Mean Squared Error
  
  # Return the results as a list
  return(list(MAE = mae_value, RMSE = rmse_value))
}

#Apply evaluation function of Hot Deck imputation
countdata5_imputedHotDeck_evaluation <- evaluate_imputationHotDeck(paste0(mypath,"countdata5_imputedKNN.txt"))
cat("Hot Deck Imputation - Mean Absolute Error (MAE):", countdata5_imputedHotDeck_evaluation$MAE, "\n")
cat("Hot Deck Imputation - Root Mean Squared Error (RMSE):", countdata5_imputedHotDeck_evaluation$RMSE, "\n")

countdata10_imputedHotDeck_evaluation <- evaluate_imputationHotDeck(paste0(mypath,"countdata10_imputedKNN.txt"))
cat("Hot Deck Imputation - Mean Absolute Error (MAE):", countdata10_imputedHotDeck_evaluation$MAE, "\n")
cat("Hot Deck Imputation - Root Mean Squared Error (RMSE):", countdata10_imputedHotDeck_evaluation$RMSE, "\n")

countdata30_imputedHotDeck_evaluation <- evaluate_imputationHotDeck(paste0(mypath,"countdata5_imputedKNN.txt"))
cat("Hot Deck Imputation - Mean Absolute Error (MAE):", countdata5_imputedHotDeck_evaluation$MAE, "\n")
cat("Hot Deck Imputation - Root Mean Squared Error (RMSE):", countdata5_imputedHotDeck_evaluation$RMSE, "\n")

#Evaluation function of KNN imputation
evaluate_imputationKNN <- function(dataset_path, na_fraction = 0.3, k = 10, seed_split = 123, seed_na = 456) {
  
  #Load the dataset from the specified file path
  dataset <- read.table(dataset_path, header = TRUE, row.names = 1)
  
  #Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed_split) #Set seed for reproducibility of train-test split
  total_samples <- nrow(dataset) #Get the total number of samples
  train_size <- floor(0.7 * total_samples) #Calculate the number of training samples (70% of total)
  
  train_indices <- sample(seq_len(total_samples), size = train_size) #Randomly sample indices for training set
  train_data <- dataset[train_indices, ] #Create training set
  test_data <- dataset[-train_indices, ] #Create testing set with the remaining data
  
  #Introduce artificial NAs randomly in the testing set
  set.seed(seed_na) #Set seed for reproducibility of NA introduction
  test_data_with_na <- test_data #Create a copy of the testing set to introduce NAs
  num_values <- prod(dim(test_data_with_na)) #Calculate the total number of values in the testing set
  num_na <- round(num_values * na_fraction) #Calculate the number of NAs to introduce based on the specified fraction
  na_indices <- arrayInd(sample(num_values, num_na), dim(test_data_with_na)) #Randomly select indices to introduce NAs
  test_data_with_na[na_indices] <- NA #Introduce NAs at the selected indices
  
  #Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.") #Stop execution if no NAs are introduced
  }
  
  #Re-impute the missing values in the testing set using kNN
  #Suppress warnings and capture output
  suppressWarnings({
    capture.output({
      test_data_reimputed <- knn_imputation(test_data_with_na, k = k) #Apply kNN imputation
    }, file = nullfile())
  })
  
  #Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.") #Stop execution if any NAs are still present after re-imputation
  }
  
  #Extract true and imputed values for comparison
  true_values <- test_data[na_indices] #Extract true values from the original testing set
  imputed_values <- test_data_reimputed[na_indices] #Extract imputed values from the re-imputed testing set
  
  #Ensure there are no NAs in the extracted true and imputed values
  if (any(is.na(true_values)) || any(is.na(imputed_values))) {
    stop("NAs present in the true or imputed values for comparison.") #Stop execution if any NAs are present in the extracted values
  }
  
  #Calculate MAE and RMSE between the imputed values and the true values in the test set
  mae_value <- mae(true_values, imputed_values) #Calculate Mean Absolute Error
  rmse_value <- rmse(true_values, imputed_values) #Calculate Root Mean Squared Error
  
  #Return the results as a list
  return(list(MAE = mae_value, RMSE = rmse_value))
}

#Apply evaluation function of KNN imputation
countdata5_imputedKNN_evaluation <- evaluate_imputationKNN(paste0(mypath,"countdata5_imputedKNN.txt"))
cat("Mean Absolute Error (MAE):", countdata5_imputedKNN_evaluation$MAE, "\n")
cat("Root Mean Squared Error (RMSE):", countdata5_imputedKNN_evaluation$RMSE, "\n")

countdata10_imputedKNN_evaluation <- evaluate_imputationKNN(paste0(mypath,"countdata10_imputedKNN.txt"))
cat("Mean Absolute Error (MAE):", countdata10_imputedKNN_evaluation$MAE, "\n")
cat("Root Mean Squared Error (RMSE):", countdata10_imputedKNN_evaluation$RMSE, "\n")

countdata30_imputedKNN_evaluation <- evaluate_imputationKNN(paste0(mypath,"countdata30_imputedKNN.txt"))
cat("Mean Absolute Error (MAE):", countdata30_imputedKNN_evaluation$MAE, "\n")
cat("Root Mean Squared Error (RMSE):", countdata30_imputedKNN_evaluation$RMSE, "\n")

#Evaluation function to perform Random Forest imputation with verbose output
RF_imputation_verbose <- function(data, maxiter, ntree) {
  for (i in seq_len(maxiter)) {
    set.seed(123 + i) # For reproducibility
    imputed_data <- missForest(data, maxiter = i, ntree = ntree, verbose = FALSE)
    cat(paste("Iteration:", i, "out of", maxiter, "completed\n"))
    print(head(imputed_data$ximp, 10)) # Print the first 10 rows of the imputed data for each iteration
  }
  return(imputed_data$ximp)
}

#Random Forest evaluation function
evaluate_imputationRF <- function(dataset_path, na_fraction = 0.3, maxiter = 5, ntree = 50, seed_split = 123, seed_na = 456) {
  
  #Load the dataset from the specified file path
  dataset <- read.table(dataset_path, header = TRUE, row.names = 1)
  
  #Split the dataset into training (70%) and testing (30%) sets
  set.seed(seed_split) #Set seed for reproducibility of train-test split
  total_samples <- nrow(dataset) #Get the total number of samples
  train_size <- floor(0.7 * total_samples) #Calculate the number of training samples (70% of total)
  
  train_indices <- sample(seq_len(total_samples), size = train_size) #Randomly sample indices for training set
  train_data <- dataset[train_indices, ] #Create training set
  test_data <- dataset[-train_indices, ] #Create testing set with the remaining data
  
  #Introduce artificial NAs randomly in the testing set
  set.seed(seed_na) #Set seed for reproducibility of NA introduction
  test_data_with_na <- test_data #Create a copy of the testing set to introduce NAs
  num_values <- prod(dim(test_data_with_na)) #Calculate the total number of values in the testing set
  num_na <- round(num_values * na_fraction) #Calculate the number of NAs to introduce based on the specified fraction
  na_indices <- arrayInd(sample(num_values, num_na), dim(test_data_with_na)) #Randomly select indices to introduce NAs
  test_data_with_na[na_indices] <- NA #Introduce NAs at the selected indices
  
  #Check if there are NAs introduced
  if (sum(is.na(test_data_with_na)) == 0) {
    stop("No NAs introduced in the test data.") # Stop execution if no NAs are introduced
  }
  
  #Perform imputation and print iteration progress
  print("Starting Random Forest imputation...")
  test_data_reimputed <- RF_imputation_verbose(test_data_with_na, maxiter = maxiter, ntree = ntree)
  
  #Ensure there are no NAs after imputation
  if (sum(is.na(test_data_reimputed)) > 0) {
    stop("NAs still present after re-imputation.") #Stop execution if any NAs are still present after re-imputation
  }
  
  #Extract true and imputed values for comparison
  true_values <- test_data[na_indices] #Extract true values from the original testing set
  imputed_values <- test_data_reimputed[na_indices] #Extract imputed values from the re-imputed testing set
  
  #Ensure there are no NAs in the extracted true and imputed values
  if (any(is.na(true_values)) || any(is.na(imputed_values))) {
    stop("NAs present in the true or imputed values for comparison.") #Stop execution if any NAs are present in the extracted values
  }
  
  #Calculate MAE and RMSE between the imputed values and the true values in the test set
  mae_value <- mae(true_values, imputed_values) #Calculate Mean Absolute Error
  rmse_value <- rmse(true_values, imputed_values) #Calculate Root Mean Squared Error
  
  #Return the results as a list
  return(list(MAE = mae_value, RMSE = rmse_value))
}

#Apply evaluation function of Random Forest imputation for Countdata5
file_path_5 <- paste0(mypath, "countdata5_imputedRF.txt")
countdata5_imputedRF_evaluation <- evaluate_imputationRF(
  dataset_path = file_path_5, 
  maxiter = 5, 
  ntree = 50
)
cat("Random Forest Imputation - Mean Absolute Error (MAE) for countdata5:", countdata5_imputedRF_evaluation$MAE, "\n")
cat("Random Forest Imputation - Root Mean Squared Error (RMSE) for countdata5:", countdata5_imputedRF_evaluation$RMSE, "\n")

#Countdata10
file_path_10 <- paste0(mypath, "countdata10_imputedRF.txt")
countdata10_imputedRF_evaluation <- evaluate_imputationRF(
  dataset_path = file_path_10, 
  maxiter = 5, 
  ntree = 50
)
cat("Random Forest Imputation - Mean Absolute Error (MAE) for countdata10:", countdata10_imputedRF_evaluation$MAE, "\n")
cat("Random Forest Imputation - Root Mean Squared Error (RMSE) for countdata10:", countdata10_imputedRF_evaluation$RMSE, "\n")

#Countdata30
file_path_30 <- paste0(mypath, "countdata30_imputedRF.txt")
countdata30_imputedRF_evaluation <- evaluate_imputationRF(
  dataset_path = file_path_30, 
  maxiter = 5, 
  ntree = 50
)
cat("Random Forest Imputation - Mean Absolute Error (MAE) for countdata30:", countdata30_imputedRF_evaluation$MAE, "\n")
cat("Random Forest Imputation - Root Mean Squared Error (RMSE) for countdata30:", countdata30_imputedRF_evaluation$RMSE, "\n")

#Visualisations of MAE and RMSE
#Create a data frame to store the evaluation results
results <- data.frame(
  Dataset = rep(c("countdata5", "countdata10", "countdata30"), each = 4),  #Repeats dataset names
  Method = rep(c("Mean", "Hot-Deck", "KNN", "Random Forest"), times = 3),  #Repeats method names
  MAE = c(
    countdata5_imputedMean_evaluation$MAE,  #MAE for countdata5 with Mean imputation
    countdata5_imputedHotDeck_evaluation$MAE,  #MAE for countdata5 with Hot-Deck imputation
    countdata5_imputedKNN_evaluation$MAE,  #MAE for countdata5 with KNN imputation
    countdata5_imputedRF_evaluation$MAE,  #MAE for countdata5 with Random Forest imputation
    countdata10_imputedMean_evaluation$MAE,  #MAE for countdata10 with Mean imputation
    countdata10_imputedHotDeck_evaluation$MAE,  #MAE for countdata10 with Hot-Deck imputation
    countdata10_imputedKNN_evaluation$MAE,  #MAE for countdata10 with KNN imputation
    countdata10_imputedRF_evaluation$MAE,  #MAE for countdata10 with Random Forest imputation
    countdata30_imputedMean_evaluation$MAE,  #MAE for countdata30 with Mean imputation
    countdata30_imputedHotDeck_evaluation$MAE,  #MAE for countdata30 with Hot-Deck imputation
    countdata30_imputedKNN_evaluation$MAE,  #MAE for countdata30 with KNN imputation
    countdata30_imputedRF_evaluation$MAE  #MAE for countdata30 with Random Forest imputation
  ),
  RMSE = c(
    countdata5_imputedMean_evaluation$RMSE,  #RMSE for countdata5 with Mean imputation
    countdata5_imputedHotDeck_evaluation$RMSE,  #RMSE for countdata5 with Hot-Deck imputation
    countdata5_imputedKNN_evaluation$RMSE,  #RMSE for countdata5 with KNN imputation
    countdata5_imputedRF_evaluation$RMSE,  #RMSE for countdata5 with Random Forest imputation
    countdata10_imputedMean_evaluation$RMSE,  #RMSE for countdata10 with Mean imputation
    countdata10_imputedHotDeck_evaluation$RMSE,  #RMSE for countdata10 with Hot-Deck imputation
    countdata10_imputedKNN_evaluation$RMSE,  #RMSE for countdata10 with KNN imputation
    countdata10_imputedRF_evaluation$RMSE,  #RMSE for countdata10 with Random Forest imputation
    countdata30_imputedMean_evaluation$RMSE,  #RMSE for countdata30 with Mean imputation
    countdata30_imputedHotDeck_evaluation$RMSE,  #RMSE for countdata30 with Hot-Deck imputation
    countdata30_imputedKNN_evaluation$RMSE,  #RMSE for countdata30 with KNN imputation
    countdata30_imputedRF_evaluation$RMSE  #RMSE for countdata30 with Random Forest imputation
  )
)

#Ensure the Dataset column is ordered
results$Dataset <- factor(results$Dataset, levels = c("countdata5", "countdata10", "countdata30"))  #Order the Dataset factor levels

#Ensure the Method column is ordered
results$Method <- factor(results$Method, levels = c("Mean", "Hot-Deck", "KNN", "Random Forest"))  #Order the Method factor levels

#Create a bar chart for MAE
mae_plot <- ggplot(results, aes(x = Dataset, y = MAE, fill = Method)) +  #Set up the plot with Dataset on x-axis, MAE on y-axis, and fill by Method
  geom_bar(stat = "identity", position = "dodge") +  #Create dodged bar chart for better comparison
  geom_text(aes(label = round(MAE, 2)), position = position_dodge(width = 0.9), vjust = -0.5) +  #Add MAE values on top of the bars
  labs(
    title = "Comparison of MAE for Different Imputation Methods",  #Title of the plot
    x = "Dataset",  #Label for x-axis
    y = "Mean Absolute Error (MAE)"  #Label for y-axis
  ) +
  theme_classic() +  #Use classic theme for a clean look
  theme(
    plot.title = element_text(size = 16, face = "bold"),  #Style the plot title
    axis.title.x = element_text(size = 14),  #Style x-axis title
    axis.title.y = element_text(size = 14),  #Style y-axis title
    axis.text.x = element_text(size = 12),  #Style x-axis text
    axis.text.y = element_text(size = 12)   #Style y-axis text
  )

#Display the MAE plot
print(mae_plot)

#Save the MAE plot as a high-quality image
output_path_mae <- paste0(mypath, "MAE_comparison_plot.png")
ggsave(filename = output_path_mae, plot = mae_plot, device = "png", dpi = 300, width = 10, height = 8)

#DESeq2 differential expression
#Function for differential expression analysis
DE_analysis <- function(countdata_file, plot_title = NULL, mypath) {
  
  #Read in the count data and metadata
  countdata <- read.table(countdata_file, header = TRUE, row.names = 1)  #Load count data
  coldata <- read.csv(paste0(mypath, "metadata.csv"))  #Load metadata
  
  #Define groups based on latitude
  coldata$group <- "North"  #Initialize all groups as "North"
  coldata$group[coldata$Latitude < 35] <- "South"  #Set "South" for latitudes below 35
  
  #Keep only the distant populations
  coldata <- rbind(coldata[coldata$Latitude < 35, ], coldata[coldata$Latitude > 40, ])  #Subset metadata for distant populations
  
  #Subset countdata to match coldata
  countdata <- countdata[, colnames(countdata) %in% coldata$RNAseq_Run_NCBI]  #Keep only columns in countdata that match the metadata
  
  #Ensure the columns of countdata correspond to the rows of coldata
  coldata <- coldata[match(colnames(countdata), coldata$RNAseq_Run_NCBI), ]  #Reorder rows of coldata to match columns of countdata
  
  #Check if the columns match after subsetting
  if (!all(colnames(countdata) == coldata$RNAseq_Run_NCBI)) {
    stop("Count data columns do not match sample information rows.")  #Stop if columns do not match
  }
  
  #Remove rows with NA values
  countdata <- na.omit(countdata)  #Remove rows with any NA values
  
  #Round the count data to ensure all values are integers
  countdata <- round(countdata)  #Round count data to nearest integers
  
  #Perform DESeq analysis with warnings suppressed
  suppressWarnings({
    
    #Create a DESeq object
    dds <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ Sex + group)  #Define the design formula
    #Perform differential expression analysis
    dds <- DESeq(dds)  #Run DESeq2 analysis
    
    # Retrieve results
    res <- results(dds)  #Extract results from DESeq2 object
    results_df <- as.data.frame(res)  #Convert results to a data frame
    results_df$Ensembl_id <- row.names(results_df)  #Add Ensembl IDs as a column
    results_df <- results_df[order(results_df$padj), ]  #Order results by adjusted p-value
    
    #Convert Ensembl names to gene names (if they are known)
    results_genes <- gconvert(row.names(res), organism = "mmusculus", target = "ENTREZGENE_ACC", filter_na = FALSE)  #Convert Ensembl to gene names
    
    #Add the gene names
    results_df <- merge(results_df,
                        results_genes[, c("input", "target", "name", "description")],
                        by.x = "Ensembl_id", by.y = "input", all.x = TRUE)  #Merge gene names with results
    
    #Ensure the 'Name' column is properly populated
    results_df$Name <- ifelse(is.na(results_df$name), results_df$Ensembl_id, results_df$name)  #Populate 'Name' column
    
    #Check for missing values in the 'Name' column
    sum(is.na(results_df$Name))  #Count missing values in 'Name' column
    
    #Subset the results to keep only significant genes
    significant_genes <- results_df[results_df$padj < 0.05 & !is.na(results_df$padj), ]  #Keep only significant genes (padj < 0.05)
    
    #Identify outliers among significant genes
    if (nrow(significant_genes) > 0) {
      
      #Filter for upregulated and downregulated genes based on log2 fold change and adjusted p-value
      upregulated_genes <- significant_genes[significant_genes$log2FoldChange > 1 & -log10(significant_genes$padj) > 5, ]
      downregulated_genes <- significant_genes[significant_genes$log2FoldChange < -1 & -log10(significant_genes$padj) > 5, ]
      
      #Print the top upregulated and downregulated genes
      cat("Top Upregulated Genes:\n")
      print(upregulated_genes)
      cat("\nTop Downregulated Genes:\n")
      print(downregulated_genes)
      
      #Plot significant genes
      message("Plotting significant genes")  #Message to indicate plotting
      
      #Volcano plot
      plot <- EnhancedVolcano(significant_genes,
                              lab = significant_genes$Name,
                              x = 'log2FoldChange',
                              y = 'padj',
                              drawConnectors = TRUE,
                              title = plot_title)  #Create volcano plot
      print(plot)  #Explicitly print the plot
      
      # Save the plot as a high-quality image
      output_path <- paste0(mypath, plot_title, "_volcano.png")  #Define output path
      ggsave(filename = output_path, plot = plot, device = "png", dpi = 300, width = 10, height = 8)  #Save the plot
      
      return(list(upregulated = upregulated_genes, downregulated = downregulated_genes))  #Return significant genes
    } else {
      message("No significant genes found.")  #Message if no significant genes are found
      return(list(upregulated = data.frame(), downregulated = data.frame()))  #Return empty data frames
    }
  })
}

#Countdata5 RF imputation
countdata5_imputedRF_DE <- DE_analysis(paste0(mypath, "countdata5_imputedRF.txt"), plot_title = "Volcano Plot of Countdata5 RF Imputation Results", mypath = mypath)

#Countdata10 RF imputation
countdata10_imputedRF_DE <- DE_analysis(paste0(mypath, "countdata10_imputedRF.txt"), plot_title = "Volcano Plot of Countdata10 RF Imputation Results", mypath = mypath)

#Countdata30 RF Imputation
countdata30_imputedRF_DE <- DE_analysis(paste0(mypath, "countdata30_imputedRF.txt"), plot_title = "Volcano Plot of Countdata30 RF Imputation Results", mypath = mypath)

#Extract gene names
upregulated_genes_5 <- countdata5_imputedRF_DE$upregulated$Name  #Extract upregulated genes from countdata5
downregulated_genes_5 <- countdata5_imputedRF_DE$downregulated$Name  #Extract downregulated genes from countdata5

upregulated_genes_10 <- countdata10_imputedRF_DE$upregulated$Name  #Extract upregulated genes from countdata10
downregulated_genes_10 <- countdata10_imputedRF_DE$downregulated$Name  #Extract downregulated genes from countdata10

#Create lists for Venn diagram
venn_list_up <- list(Countdata5 = upregulated_genes_5, Countdata10 = upregulated_genes_10)  #List for upregulated genes
venn_list_down <- list(Countdata5 = downregulated_genes_5, Countdata10 = downregulated_genes_10)  #List for downregulated genes

#Plot Venn diagram for upregulated genes
venn_plot_up <- ggvenn(venn_list_up, show_elements = TRUE, label_sep = "\n",
                       fill_color = c("red", "orange")) +  #Create Venn diagram for upregulated genes
  ggtitle("Upregulated Genes") +  #Title for the Venn diagram
  theme_minimal() +  #Use minimal theme
  theme(
    panel.grid = element_blank(),  #Remove panel grid
    panel.background = element_blank(),  #Remove panel background
    plot.background = element_blank(),  #Remove plot background
    axis.line = element_blank(),  #Remove axis line
    axis.text = element_blank(),  #Remove axis text
    axis.ticks = element_blank(),  #Remove axis ticks
    axis.title = element_blank(),  #Remove axis title
    legend.position = "none",  #Remove legend
    plot.title = element_text(size = 30)  #Increase title size
  ) +
  scale_size_area(max_size = 45)  #Increase circle size

#Plot Venn diagram for downregulated genes
venn_plot_down <- ggvenn(venn_list_down, show_elements = TRUE, label_sep = "\n",
                         fill_color = c("green", "blue")) +  #Create Venn diagram for downregulated genes
  ggtitle("Downregulated Genes") +  #Title for the Venn diagram
  theme_minimal() +  #Use minimal theme
  theme(
    panel.grid = element_blank(),  #Remove panel grid
    panel.background = element_blank(),  #Remove panel background
    plot.background = element_blank(),  #Remove plot background
    axis.line = element_blank(),  #Remove axis line
    axis.text = element_blank(),  #Remove axis text
    axis.ticks = element_blank(),  #Remove axis ticks
    axis.title = element_blank(),  #Remove axis title
    legend.position = "none",  #Remove legend
    plot.title = element_text(size = 30)  #Increase title size
  ) +
  scale_size_area(max_size = 45)  #Increase circle size

#Combine the two Venn diagrams into one plot
combined_plot <- grid.arrange(venn_plot_up, venn_plot_down, ncol = 2)

#Print and save the combined Venn diagram
output_path_combined <- paste0(mypath, "combined_venn_diagram.png")  #Define output path
ggsave(filename = output_path_combined, plot = combined_plot, device = "png", dpi = 300, width = 20, height = 10)  #Save the combined Venn diagram

#The End!