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
JO_data <- JO_data[!(is.na(JO_data$Hog)),] # Remove rows with NA valus in Hog

# Average the coordinates for each unique ID in each year. This is because 
# the coordinates change slightly year to year for each ID, due to measurement error
for (i in 1:(length(JO_data$x))){
  JO_data$x2[i] <- mean(JO_data$x[which(JO_data$ID == JO_data$ID[i])])
  JO_data$y2[i] <- mean(JO_data$y[which(JO_data$ID == JO_data$ID[i])])
}

## In this step, we find out if there are duplicates.
library(tidyverse)
library(dplyr)
JO_data3 <- as_tibble(JO_data)
JO_data3 <- JO_data3 %>% distinct(ID, .keep_all = TRUE)
JO_data3 <- JO_data3 %>% arrange(ID)
JO_data4 <- JO_data3 %>% distinct(x2, y2, .keep_all = TRUE)

JO_data3[which(duplicated(JO_data3$x2)),]$ID
JO_data3$ID
JO_data4$ID

### CREATE MODEL FRAME / SETUP MATRICES ###

# We remove the IDs that are in the same location as another
JO_data <- JO_data[-which(JO_data$ID == 179),] #179 is a duplicate location.

IDs <- c(unique(JO_data$ID))
length(IDs) #410 unique IDs
ranked_IDs <- rank(IDs)

# Replace the old ID with a sequence of IDs from 1 to length(unique_IDs).
# This is done to ease construction of the Z matrix for the augmented part.
for (i in 1:(length(JO_data$x))){
  for (j in 1:length(IDs)){
    if (JO_data$ID[i] == IDs[j]){
      JO_data$ranked_id[i] <- ranked_IDs[j]
    }
  }
}

IDs_check <- c(unique(JO_data$ranked_id))
length(IDs_check) #410 unique IDs
max(IDs_check)

plot(JO_data$x2,JO_data$y2) #Averaged coordinates for each ID

# Coordinates
x.coords <- JO_data$x2
y.coords <- JO_data$y2
uid <- JO_data$ranked_id

aggr_JO_data <- aggregate(x = JO_data, by = list(JO_data$ranked_id), FUN = "mean")
