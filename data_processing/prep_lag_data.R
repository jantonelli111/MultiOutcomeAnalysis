# Load packages
library(tidyverse)
library(arrow)

## Read in fake data
#dataOriginal = read.csv("~/Documents/Research/Suyeon/FakeMedicareData/outcomes_merged_all_years.csv")
## Read data
dataOriginal = read_parquet("data/processed/outcomes_merged_all_years.parquet")

## Convert all variables to numeric (some are integers)
data = data.frame(sapply(dataOriginal, as.numeric))

## Need to create future exposures and prior outcomes variables for 
## Use as negative controls
data$FuturePM25 = rep(NA, nrow(data))
data$FutureOzone = rep(NA, nrow(data))
data$FutureNH4 = rep(NA, nrow(data))
data$FutureSO4 = rep(NA, nrow(data))
data$FutureOC = rep(NA, nrow(data))
data$FutureEC = rep(NA, nrow(data))
data$FutureNO3 = rep(NA, nrow(data))

data$PriorMortality = rep(NA, nrow(data))
data$PriorCOPD = rep(NA, nrow(data))

uniqueZips = unique(data$zip)
uniqueYears = unique(data$year)

for (tempZip in uniqueZips) {
  for (tempYear in uniqueYears) {
    
    ## Finding indices of data frame that correspond to this 
    ## particular zip code and this particular year as well as
    ## previous and future year
    wCurrent = which(data$zip == tempZip &
                       data$year == tempYear)
    wPast = which(data$zip == tempZip &
                    data$year == tempYear-3)
    wFuture = which(data$zip == tempZip &
                      data$year == tempYear+3)
    
    ## only proceed if each of these is in the data set
    if (length(wCurrent) == 1 &
        length(wPast) == 1 &
        length(wFuture) == 1) {
      
      data$PriorMortality[wCurrent] = data$mort_rate[wPast]
      data$PriorCOPD[wCurrent] = data$copd_rate[wPast]
      
      data$FutureEC[wCurrent] = data$ec[wFuture]
      data$FutureOzone[wCurrent] = data$ozone[wFuture]
      data$FutureOC[wCurrent] = data$oc[wFuture]
      data$FutureNH4[wCurrent] = data$nh4[wFuture]
      data$FutureSO4[wCurrent] = data$so4[wFuture]
      data$FutureNO3[wCurrent] = data$no3[wFuture]
      data$FuturePM25[wCurrent] = data$pm25[wFuture]
    }
  }
}

## Remove rows that have NAs for variables of interest
variablesOfInterest = c("zip",
                        "year",
                        "mean_bmi",
                        "smoke_rate",
                        "medhouseholdincome",
                        "medianhousevalue",
                        "poverty",
                        "education",
                        "pct_owner_occ",
                        "summer_tmmx",
                        "winter_tmmx",
                        "summer_rmax",
                        "winter_rmax",
                        "pm25",
                        "ozone",
                        "ec",
                        "nh4",
                        "no3",
                        "oc",
                        "so4",
                        "mort_rate",
                        "female_pct",
                        "dual_pct",
                        "mean_age",
                        "race_white_pct",
                        "race_black_pct",
                        "race_hispanic_pct",
                        "anemia_rate",
                        "copd_rate",
                        "stroke_rate",
                        "lungCancer_rate",
                        "asthma_rate",
                        "hypert_rate",
                        "FuturePM25",
                        "FutureOzone",
                        "FutureNH4",
                        "FutureSO4",
                        "FutureOC",
                        "FutureEC",
                        "FutureNO3",
                        "PriorMortality",
                        "PriorCOPD")


## Remove NAs
wNotNA = which(complete.cases(data[,variablesOfInterest]) == TRUE)

## Store data frame without NAs
finalData = data[wNotNA,]
finalData = t(data.frame(do.call(rbind, finalData))) # convert list to dataframe

write_rds(finalData, "data/processed/data.rds")
