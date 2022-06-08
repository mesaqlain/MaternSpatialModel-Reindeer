# MaternSpatialModel-Reindeer
R Scripts to estimate parameters for Reindeer Dataset (Skarin and Alam 2017), results included in original research manuscript by Saqlain, Rönnegårda, Alam, and Skarin

Spatial Analysis of Reindeer data (Skarin & Alam 2017)

To perform analysis on the Reindeer data, need the following:

## R packages used:
- Matrix
- lgcp (for circulant() function)
- tidyverse (for data cleaning)
- dplyr (for data cleaning)

## Scripts with functions:

- **Create_Z_functions_V6.R**: Source file where the functions to create Z are
- **Create_Q_functions_V5.R**: Source file where the functions to create Q are
- **Estimate_function_V6_binomial.R**: Source file where the functions to estimate the parameters are
- **Check_Q_and_Matern_functions_V4.R**: Source file where the functions to create regular lattice is
- **LogL_function_V7_binom.R**: Source file for the Log likelihood function
- **HL_Correction.R**: HL Correction file

## Scripts to run analysis on Reindeer dataset:
- ReindeerDataset.R: Script to load and clean the dataset
- ST_Prepare_Data.R: Script to clean and setup the dataset for ST region only, details in the script
- ST_All_Years_NoPrecip.R: Script to estimate parameters for ST region only, details in the script
- JO_Prepare_Data.R: Script to clean and setup the dataset for JO region only, details in the script
- JO_All_Years_NoPrecip.R: Script to estimate parameters for JO region only, details in the script

## Dataset file:
- Pellet_2010_2015_2.RData

## Directions:
**Step 1:**
- Load the ReindeerDataset.R file. 
- Setwd() in Line 21 to the location where the Pellet_2010_2015_2.RData file is stored. 
- The script subsets the dataset to ST_data which has data only from ST region, JO_data which has data only from JO region, and REG_data which has only the regional data. ST_data is further cleaned to form the ST_data_cleaned.
- No need to change anything, run the entire script.

**Step 2:**
- Load the ST_Prepare_Data.R file. 
- Setwd() in Line 10 to the location where the necessary scripts (shown in Lines 12-17) are stored.
- No need to change anything, run the entire script.

**Step 3:**
- Load the ST_All_Years_NoPrecip.R file.
- Line 6: lattice_spacing can be changed to define how fine the underlying mesh should be.
- Line 7: extension_points can be extend the mesh location beyond the observation area to avoid boundary effects.
- Line 16: Change range to a list of ranges to perform analysis for, currently set to c(350, 400). 
- Line 17. Alpha = 2 or 3, analysis was performed for 2.
- Line 50: Initial values b0, taken from Skarin and Alam 2017 results
- Outputs: ST_results2.csv and ST_results_best2.csv

**Step 4:**
- Load file ST_results_best2.csv
- This has the results of our analysis based on which range value yielded the maximum likelihood from the list of range values we had provided. Open ST_results2.csv to check results for ALL range values.
- Column 1: b0_best are the estimated parameters, these are in order:
  
  #[1] "(Intercept)"                                                                    
  #[2] "sqrt(distvindkraft/100)"                                          
  #[3] "factor(year_descrip)construction"                         
  #[4] "factor(year_descrip)operation"                          
  #[5] "factor(smd50upveg3)Forest"                             
  #[6] "factor(smd50upveg3)Clear"                              
  #[7] "factor(smd50upveg3)Young"                                                   
  #[8] "sqrt(distvindkraft/100):factor(year_descrip)construction" 
  #[9] "sqrt(distvindkraft/100):factor(year_descrip)operation"     
- Column 2: Estimate of s0 / tau parameter
- Column 3: log likelihood value, this is the maximum likelihood from the ranges we had provided in ST_All_Years_NoPrecip.R Line 16.
- Column 4: Range value where maximum likelihood was found

## JO DATASET Analysis:
- For JO Dataset, repeat Steps 2-4 but using **JO_Prepare_Data.R** and **JO_All_Years_NoPrecip.R** instead.

## Contact: 
- Murshid Saqlain (msq@du.se) for any questions related to the code.

## References:
- Lindgren, F., Rue, H. and Lindström, J., 2011. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(4), pp.423-498.
- Skarin, A. and Alam, M., 2017. Reindeer habitat use in relation to two small wind farms, during preconstruction, construction, and operation. Ecology and Evolution, 7(11), pp.3870-3882.
