## Written by Murshid Saqlain
## Script to prepare ST dataset to be used for analysis in ST_All_Years_NoPrecip.R
## Change working directory to where the scripts are stored in Line 10.

# Load necessary packages
library(Matrix)
library(lgcp) #For circulant() function

#Where the functions are located
setwd("~/Paper 3/Version 6/")

source("Create_Z_functions_V6.R") #Source file where the functions to create Z are
source("Create_Q_functions_V5.R") #Source file where the functions to create Q are
source("Estimate_function_V6_binomial.R") #Source file where the functions to estimate the parameters are
source("Check_Q_and_Matern_functions_V4.R") #Source file where the functions to create regular lattice is
source("LogL_function_V7_binom.R") #Source file for the Log likelihood function
source("HL_Correction.R") #HL Correction file

### CREATE THE Q MATRIX ###
ST_data <- ST_data_cleaned[!(is.na(ST_data_cleaned$Hog)),] # Remove rows with NA valus in Hog

# Average the coordinates for each unique ID in each year. This is because 
# the coordinates change slightly year to year for each ID, due to measurement error
for (i in 1:(length(ST_data$x))){
  ST_data$x2[i] <- mean(ST_data$x[which(ST_data$ID == ST_data$ID[i])])
  ST_data$y2[i] <- mean(ST_data$y[which(ST_data$ID == ST_data$ID[i])])
}

## In this step, we find out if there are duplicates.
library(tidyverse)
library(dplyr)
ST_data3 <- as_tibble(ST_data)
ST_data3 <- ST_data3 %>% distinct(ID, .keep_all = TRUE)
ST_data3 <- ST_data3 %>% arrange(ID)
ST_data4 <- ST_data3 %>% distinct(x2, y2, .keep_all = TRUE)

IDs <- c(unique(ST_data$ID))
length(IDs) #357 unique IDs
ranked_IDs <- rank(IDs)
max(ranked_IDs)
# Replace the old ID with a sequence of IDs from 1 to length(unique_IDs).
# This is done to ease construction of the Z matrix for the augmented part.
for (i in 1:(length(ST_data$x))){
  for (j in 1:length(IDs)){
    if (ST_data$ID[i] == IDs[j]){
      ST_data$ranked_id[i] <- ranked_IDs[j]
    }
  }
}

IDs_check <- c(unique(ST_data$ranked_id))
length(IDs_check) #357 unique IDs
max(IDs_check)

plot(ST_data$x2,ST_data$y2) #Averaged coordinates for each ID

# Coordinates
x.coords <- ST_data$x2
y.coords <- ST_data$y2
uid <- ST_data$ranked_id
## Data is now ready to be used for parameter estimation.
