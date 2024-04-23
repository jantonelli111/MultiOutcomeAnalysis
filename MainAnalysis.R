
rm(list=ls())

###########################################################
## Code for analysis of multiple treatments              ##
## and multiple outcomes simultaneously                  ##
##                                                       ##
##                                                       ##
##                                                       ##
## Inputs to the function:                               ##
##                                                       ##
## 1: Y: an N by Q matrix of outcomes                    ##
## 2: Tr: an N by P matrix of exposures                  ##
## 3: X: an N by C matrix of confounders                 ##
## 4: b: J by Q matrix for NC contrasts                  ##
## 5: t1: exposure level 1 for estimand                  ##
## 6: t2: exposure level 2 for estimand                  ##
## 7: t1NC: list of exposure 1 values for NC             ##
## 8: t2NC: list of exposure 2 values for NC             ##
###########################################################

# Load packages
require(earth) # for nonlinear multivariate regression model (MARS; Friedman, 1991)
require(MASS)  # for generalized inverse of matrix: ginv()
require(psych) # for the estimation of unobserved confounders: vss(), fa.parallel()

## Read data -----

## when using the fake data 
dataOriginal = read.csv("data/aux/fake_data.csv")

## when using the real data
# library(arrow)
# dataOriginal = arrow::read_parquet("data/data.parquet")

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
                    data$year == tempYear-1)
    wFuture = which(data$zip == tempZip &
                      data$year == tempYear+1)
    
    ## only proceed if each of these is in the data set
    if (length(wCurrent) == 1 &
        length(wPast) == 1 &
        length(wFuture) == 1) {
      
      data$PriorMortality[wCurrent] = data$death_rate[wPast]
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
                        "death_rate",
                        "female_pct",
                        "dual_pct",
                        "mean_age",
                        "race_black_pct",
                        "race_white_pct",
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



##----------------------------------------------------------------------------##
## ANALYSIS 0: Original setting
##----------------------------------------------------------------------------##
## Now create variables needed for our analysis
Y = finalData[,c("death_rate",
                 "anemia_rate",
                 "copd_rate",
                 "stroke_rate",
                 "lungCancer_rate",
                 "asthma_rate",
                 "hypert_rate",
                 "PriorMortality",
                 "PriorCOPD")]

Tr = finalData[,c("pm25",
                  "ozone",
                  "ec",
                  "nh4",
                  "no3",
                  "oc",
                  "so4",
                  "FuturePM25",
                  "FutureOzone")]

X = finalData[,c("year",
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
                 "female_pct",
                 "dual_pct",
                 "mean_age",
                 "race_black_pct",
                 "race_white_pct",
                 "race_hispanic_pct")]

## estimand of interest
t1 = list()
t2 = list()
t1[[1]] = c(apply(Tr[,1:7], 2, quantile, 0.75), apply(Tr[,8:9], 2, median))
t2[[1]] = c(apply(Tr[,1:7], 2, quantile, 0.25), apply(Tr[,8:9], 2, median))

for (tt in 2 : 8) {
  t1[[tt]] = apply(Tr, 2, median)
  t1[[tt]][tt-1] = quantile(Tr[,tt-1], 0.75)
  t2[[tt]] = apply(Tr, 2, median)
  t2[[tt]][tt-1] = quantile(Tr[,tt-1], 0.25)
}

## Negative control outcomes
b = matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 1),
           byrow = FALSE, nrow = ncol(Y))

## Exposure contrasts for negative controls
t1NC = list()
t2NC = list()

## First use future exposures as NC for current mortality
t1NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 2), ncol=2))
t2NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 2), ncol=2))

for (tt in 8 : 9) {
  t1NC[[1]][tt-7,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[1]][tt-7,tt] = quantile(Tr[,tt], 0.25)
}

## Now use current exposures with prior outcomes
t1NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[2]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[2]][tt,tt] = quantile(Tr[,tt], 0.25)
}

t1NC[[3]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[3]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[3]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[3]][tt,tt] = quantile(Tr[,tt], 0.25)
}

maxMfunc = function(qp, m){
  (qp - m)^2 - qp - m
}

for(pq in 3:min(ncol(Y), ncol(Tr))){
  m = 0
  while(maxMfunc(pq, m) >= 0){m = m + 1}
  m = m - 1
}

A0 = list(Y = Y, Tr = Tr, X = X,
          t1 = t1, t2 = t2,
          t1NC = t1NC, t2NC = t2NC, b = b,
          maxM = m)


##----------------------------------------------------------------------------##
## ANALYSIS 1: No NCEs, only keeping NCOs
##----------------------------------------------------------------------------##
## Now create variables needed for our analysis
Y = finalData[,c("death_rate",
                 "anemia_rate",
                 "copd_rate",
                 "stroke_rate",
                 "lungCancer_rate",
                 "asthma_rate",
                 "hypert_rate",
                 "PriorMortality",
                 "PriorCOPD")]

Tr = finalData[,c("pm25",
                  "ozone",
                  "ec",
                  "nh4",
                  "no3",
                  "oc",
                  "so4")]

X = finalData[,c("year",
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
                 "female_pct",
                 "dual_pct",
                 "mean_age",
                 "race_black_pct",
                 "race_white_pct",
                 "race_hispanic_pct")]

## estimand of interest
t1 = list()
t2 = list()
t1[[1]] = apply(Tr[,1:7], 2, quantile, 0.75)
t2[[1]] = apply(Tr[,1:7], 2, quantile, 0.25)

for (tt in 2 : 8) {
  t1[[tt]] = apply(Tr, 2, median)
  t1[[tt]][tt-1] = quantile(Tr[,tt-1], 0.75)
  t2[[tt]] = apply(Tr, 2, median)
  t2[[tt]][tt-1] = quantile(Tr[,tt-1], 0.25)
}

## Negative control outcomes
b = matrix(c(0, 0, 0, 0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 1),
           byrow = FALSE, nrow = ncol(Y))

## Exposure contrasts for negative controls
t1NC = list()
t2NC = list()

## Now use current exposures with prior outcomes
t1NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[1]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[1]][tt,tt] = quantile(Tr[,tt], 0.25)
}

t1NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[2]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[2]][tt,tt] = quantile(Tr[,tt], 0.25)
}

maxMfunc = function(qp, m){
  (qp - m)^2 - qp - m
}

for(pq in 3:min(ncol(Y), ncol(Tr))){
  m = 0
  while(maxMfunc(pq, m) >= 0){m = m + 1}
  m = m - 1
}

A1 = list(Y = Y, Tr = Tr, X = X,
          t1 = t1, t2 = t2,
          t1NC = t1NC, t2NC = t2NC, b = b,
          maxM = m)

##----------------------------------------------------------------------------##
## ANALYSIS 2: ANALYSIS 0 but 3 years ahead and back in time
##----------------------------------------------------------------------------##
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
      
      data$PriorMortality[wCurrent] = data$death_rate[wPast]
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

## Remove NAs
wNotNA = which(complete.cases(data[,variablesOfInterest]) == TRUE)

## Store data frame without NAs
finalData = data[wNotNA,]
finalData = t(data.frame(do.call(rbind, finalData))) # convert list to dataframe

## Now create variables needed for our analysis
Y = finalData[,c("death_rate",
                 "anemia_rate",
                 "copd_rate",
                 "stroke_rate",
                 "lungCancer_rate",
                 "asthma_rate",
                 "hypert_rate",
                 "PriorMortality",
                 "PriorCOPD")]

Tr = finalData[,c("pm25",
                  "ozone",
                  "ec",
                  "nh4",
                  "no3",
                  "oc",
                  "so4",
                  "FuturePM25",
                  "FutureOzone")]

X = finalData[,c("year",
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
                 "female_pct",
                 "dual_pct",
                 "mean_age",
                 "race_black_pct",
                 "race_white_pct",
                 "race_hispanic_pct")]

## estimand of interest
t1 = list()
t2 = list()
t1[[1]] = c(apply(Tr[,1:7], 2, quantile, 0.75), apply(Tr[,8:9], 2, median))
t2[[1]] = c(apply(Tr[,1:7], 2, quantile, 0.25), apply(Tr[,8:9], 2, median))

for (tt in 2 : 8) {
  t1[[tt]] = apply(Tr, 2, median)
  t1[[tt]][tt-1] = quantile(Tr[,tt-1], 0.75)
  t2[[tt]] = apply(Tr, 2, median)
  t2[[tt]][tt-1] = quantile(Tr[,tt-1], 0.25)
}

## Negative control outcomes
b = matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 1),
           byrow = FALSE, nrow = ncol(Y))

## Exposure contrasts for negative controls
t1NC = list()
t2NC = list()

## First use future exposures as NC for current mortality
t1NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 2), ncol=2))
t2NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 2), ncol=2))

for (tt in 8 : 9) {
  t1NC[[1]][tt-7,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[1]][tt-7,tt] = quantile(Tr[,tt], 0.25)
}

## Now use current exposures with prior outcomes
t1NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[2]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[2]][tt,tt] = quantile(Tr[,tt], 0.25)
}

t1NC[[3]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[3]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[3]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[3]][tt,tt] = quantile(Tr[,tt], 0.25)
}

maxMfunc = function(qp, m){
  (qp - m)^2 - qp - m
}

for(pq in 3:min(ncol(Y), ncol(Tr))){
  m = 0
  while(maxMfunc(pq, m) >= 0){m = m + 1}
  m = m - 1
}

A2 = list(Y = Y, Tr = Tr, X = X,
          t1 = t1, t2 = t2,
          t1NC = t1NC, t2NC = t2NC, b = b,
          maxM = m)

##----------------------------------------------------------------------------##
## ANALYSIS 3: Spatial NCs
##----------------------------------------------------------------------------##
data = data.frame(sapply(dataOriginal, as.numeric))

## Need to create future exposures and prior outcomes variables for 
## Use as negative controls
data$DistantPM25 = rep(NA, nrow(data))
data$DistantOzone = rep(NA, nrow(data))
data$DistantNH4 = rep(NA, nrow(data))
data$DistantSO4 = rep(NA, nrow(data))
data$DistantOC = rep(NA, nrow(data))
data$DistantEC = rep(NA, nrow(data))
data$DistantNO3 = rep(NA, nrow(data))

data$DistantMortality = rep(NA, nrow(data))
data$DistantCOPD = rep(NA, nrow(data))

uniqueZips = unique(data$zip)
uniqueYears = unique(data$year)

zipBins = substr(data$zip, 1, nchar(data$zip) - 4)
zipBins = as.numeric(replace(zipBins, zipBins == "", "0"))
#DistantzipBins = as.numeric(zipBins) + 5
#DistantzipBins = substr(DistantzipBins,
#                        nchar(DistantzipBins),
#                        nchar(DistantzipBins))

set.seed(325)
for (tempZip in uniqueZips) {
  
  ZipBinCurrent = substr(tempZip, 1, nchar(tempZip) - 4)
  ZipBinCurrent = as.numeric(replace(ZipBinCurrent, ZipBinCurrent == "", "0"))
  ZipBinDistant = ifelse(ZipBinCurrent >= 5, ZipBinCurrent - 5, ZipBinCurrent + 5)
  wDistantCandidate = which(ZipBinDistant == zipBins)
  wDistantRd = sample(wDistantCandidate, 1)
  
  for (tempYear in uniqueYears) {
    
    #dataZipBin = substr(data$zip, 1, nchar(data$zip) - 4)
    #dataZipBin = as.numeric(replace(dataZipBin, dataZipBin == "", "0"))
    wCurrent = which(data$zip == tempZip &
                       data$year == tempYear)
    wDistant = which(data$zip == data$zip[wDistantRd] &
                       data$year == tempYear)
    #print(c(tempZip, tempYear, data$zip[wDistantRd]))
    #print(c(tempZip, tempYear, data$copd_rate[wCurrent], data$copd_rate[wDistant]))
    
    ## only proceed if each of these is in the data set
    if (length(wCurrent) == 1 &
        length(wDistant) == 1) {
      
      data$DistantMortality[wCurrent] = data$death_rate[wDistant]
      data$DistantCOPD[wCurrent] = data$copd_rate[wDistant]
      
      data$DistantEC[wCurrent] = data$ec[wDistant]
      data$DistantOzone[wCurrent] = data$ozone[wDistant]
      data$DistantOC[wCurrent] = data$oc[wDistant]
      data$DistantNH4[wCurrent] = data$nh4[wDistant]
      data$DistantSO4[wCurrent] = data$so4[wDistant]
      data$DistantNO3[wCurrent] = data$no3[wDistant]
      data$DistantPM25[wCurrent] = data$pm25[wDistant]
    }
  }
}
#View(finalData[,c("zip","year","copd_rate","DistantCOPD")])

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
                        "death_rate",
                        "female_pct",
                        "dual_pct",
                        "mean_age",
                        "race_black_pct",
                        "race_white_pct",
                        "race_hispanic_pct",
                        "anemia_rate",
                        "copd_rate",
                        "stroke_rate",
                        "lungCancer_rate",
                        "asthma_rate",
                        "hypert_rate",
                        "DistantPM25",
                        "DistantOzone",
                        "DistantNH4",
                        "DistantSO4",
                        "DistantOC",
                        "DistantEC",
                        "DistantNO3",
                        "DistantMortality",
                        "DistantCOPD")

## Remove NAs
wNotNA = which(complete.cases(data[, variablesOfInterest]) == TRUE)

## Store data frame without NAs
finalData = data[wNotNA,]
finalData = t(data.frame(do.call(rbind, finalData))) # convert list to dataframe

## Now create variables needed for our analysis
Y = finalData[,c("death_rate",
                 "anemia_rate",
                 "copd_rate",
                 "stroke_rate",
                 "lungCancer_rate",
                 "asthma_rate",
                 "hypert_rate",
                 "DistantMortality",
                 "DistantCOPD")]

Tr = finalData[,c("pm25",
                  "ozone",
                  "ec",
                  "nh4",
                  "no3",
                  "oc",
                  "so4",
                  "DistantPM25",
                  "DistantOzone")]

X = finalData[,c("year",
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
                 "female_pct",
                 "dual_pct",
                 "mean_age",
                 "race_black_pct",
                 "race_white_pct",
                 "race_hispanic_pct")]

## estimand of interest
t1 = list()
t2 = list()
t1[[1]] = c(apply(Tr[,1:7], 2, quantile, 0.75), apply(Tr[,8:9], 2, median))
t2[[1]] = c(apply(Tr[,1:7], 2, quantile, 0.25), apply(Tr[,8:9], 2, median))

for (tt in 2 : 8) {
  t1[[tt]] = apply(Tr, 2, median)
  t1[[tt]][tt-1] = quantile(Tr[,tt-1], 0.75)
  t2[[tt]] = apply(Tr, 2, median)
  t2[[tt]][tt-1] = quantile(Tr[,tt-1], 0.25)
}

## Negative control outcomes
b = matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 1),
           byrow = FALSE, nrow = ncol(Y))

## Exposure contrasts for negative controls
t1NC = list()
t2NC = list()

## First use future exposures as NC for current mortality
t1NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 2), ncol=2))
t2NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 2), ncol=2))

for (tt in 8 : 9) {
  t1NC[[1]][tt-7,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[1]][tt-7,tt] = quantile(Tr[,tt], 0.25)
}

## Now use current exposures with prior outcomes
t1NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[2]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[2]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[2]][tt,tt] = quantile(Tr[,tt], 0.25)
}

t1NC[[3]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[3]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[3]][tt,tt] = quantile(Tr[,tt], 0.75)
  t2NC[[3]][tt,tt] = quantile(Tr[,tt], 0.25)
}

maxMfunc = function(qp, m){
  (qp - m)^2 - qp - m
}

for(pq in 3:min(ncol(Y), ncol(Tr))){
  m = 0
  while(maxMfunc(pq, m) >= 0){m = m + 1}
  m = m - 1
}

A3 = list(Y = Y, Tr = Tr, X = X,
          t1 = t1, t2 = t2,
          t1NC = t1NC, t2NC = t2NC, b = b,
          maxM = m)

source("MainFunctions.R")
system.time(
  testALL <- wrapFunc(SettingList = list(ANALYSIS0 = A0,
                                         ANALYSIS1 = A1,
                                         ANALYSIS2 = A2,
                                         ANALYSIS3 = A3,
                                         scaleData = TRUE,
                                         nB = 50))
)

## When using the fake data
save(testALL, file="data/output/fake_OutputSaved.dat")

## When using the real data
#save(testALL, file="data/output/OutputSaved.dat")
