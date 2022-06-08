## Written by Murshid Saqlain
## Script to load the Reindeer data into memory and perform some preliminary cleaning.
## Change Line 21 to working directory where the Pellet_2010_2015_2.RData 
## file is stored. 
## Rest of the script performs some preliminary cleaning and subsets the dataset to 
## ST_data which has data only from ST region, 
## JO_data which has data only from JO region, and
## REG_data which has only the regional data. 
## ST_data is further cleaned to form the ST_data_cleaned.

#I plan to subset the dataset (Pellet_2010_2015_2.RData) into
#the three different regions. The regions are specified by the 
#"omr" variable. ST for Storliden, JO for Jokkmokksliden, and 
#REG for the regional observations (square plots on the map). 
#For ST and JO, the observations are in transects that are 300m apart, 
#and the plots in each transects are 100m apart. 

rm(list=ls())

#Where the functions are located
setwd("~/Paper 3/Reindeer_Paper/")

#Load the dataset
load("~/Pellet_2010_2015_2.RData", envir = parent.frame(), verbose = FALSE)
dataset <- pellet.vari

# Check how many observations for each region type:
table(dataset$omr)
#ST   JO  REG 
#2209 2459 3221 


## Add PRECIPITATION TO DATASET
#Year	AvgPrec
#2010	71.68333333
#2011	70.2
#2012	93.21666667
#2013	87.31666667
#2014	65.41666667
#2015	46.66666667

dataset$precip[dataset$year==2010] <- 71.68333333
dataset$precip[dataset$year==2011] <- 70.2
dataset$precip[dataset$year==2012] <- 93.21666667
dataset$precip[dataset$year==2013] <- 87.31666667
dataset$precip[dataset$year==2014] <- 65.41666667
dataset$precip[dataset$year==2015] <- 46.66666667




xtabs(~dataset$precip)

# Subset the data for only ST region
ST_data <- subset(dataset, omr %in% c("ST"))
table(ST_data$omr)


# Subset the data for only JO region
JO_data <- subset(dataset, omr %in% c("JO"))

# Subset the data for only REG 
REG_data <- subset(dataset, omr %in% c("REG"))

plot(ST_data$x,ST_data$y)
plot(JO_data$x,JO_data$y)
plot(REG_data$x,REG_data$y)


## Clean the ST Subset
ST_data_cleaned <- ST_data[which(ST_data$x < 685000 & ST_data$y < 7239700),]
plot(ST_data_cleaned$x,ST_data_cleaned$y)
