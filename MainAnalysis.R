
rm(list=ls())

# define parser arguments ----
# args <- list()
# args$fake <- TRUE
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--fake", action = "store_true", help = "use fake data")               
args = parser$parse_args()

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
require(data.table) # for data.table functions to improve efficiency

## Read data -----
if(args$fake){
  dataOriginal = fread("data/aux/fake_data.csv")
} else {
  library(arrow)
  dataOriginal = arrow::read_parquet("data/processed/data.parquet")
}
data = dataOriginal

## Convert data to a data.table
setDT(data)
## Convert all columns to numeric (handles integers fine)
data[, (names(data)) := lapply(.SD, as.numeric)]

## Need to create future exposures and prior outcomes variables for 
## Use as negative controls
newCols <- c("FuturePM25", "FutureOzone",
             "PriorMortality", "PriorCOPD")
data[, (newCols) := .(rep(NA_real_, .N), rep(NA_real_, .N),
                      rep(NA_real_, .N), rep(NA_real_, .N))]

## Create shifts in year columns for past and future matching
data[, c("YearPlus1", "YearMinus1") := .(year + 1, year - 1)]

## Create a key for faster subsetting and joins
setkey(data, zip, year)

## Join data with itself to fetch future and past values
dataFuture <- data[, .(zip, year = YearPlus1,
                       FuturePM25 = pm25, FutureOzone = ozone)]
dataPast <- data[, .(zip, year = YearMinus1,
                     PriorMortality = death_rate, PriorCOPD = copd_rate)]

## Merge future and past data into the original dataset by reference
data[dataFuture, on = .(zip, year),
     `:=` (FuturePM25 = i.FuturePM25, FutureOzone = i.FutureOzone)]
data[dataPast, on = .(zip, year),
     `:=` (PriorMortality = i.PriorMortality, PriorCOPD = i.PriorCOPD)]

## Clean up helper columns
data[, `:=` (YearPlus1 = NULL, YearMinus1 = NULL)]

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
                        "PriorMortality",
                        "PriorCOPD")


## Store data frame without NAs
finalData = na.omit(data, cols = variablesOfInterest)
setDT(finalData)



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
data = dataOriginal

## Convert data to a data.table
setDT(data)
## Convert all columns to numeric (handles integers fine)
data[, (names(data)) := lapply(.SD, as.numeric)]

## Need to create future exposures and prior outcomes variables for 
## Use as negative controls
newCols <- c("FuturePM25", "FutureOzone",
             "PriorMortality", "PriorCOPD")
data[, (newCols) := .(rep(NA_real_, .N), rep(NA_real_, .N),
                      rep(NA_real_, .N), rep(NA_real_, .N))]

## Create shifts in year columns for past and future matching
data[, c("YearPlus", "YearMinus") := .(year + 3, year - 3)]

## Create a key for faster subsetting and joins
setkey(data, zip, year)

## Join data with itself to fetch future and past values
dataFuture <- data[, .(zip, year = YearPlus,
                       FuturePM25 = pm25, FutureOzone = ozone)]
dataPast <- data[, .(zip, year = YearMinus,
                     PriorMortality = death_rate, PriorCOPD = copd_rate)]

## Merge future and past data into the original dataset by reference
data[dataFuture, on = .(zip, year),
     `:=` (FuturePM25 = i.FuturePM25, FutureOzone = i.FutureOzone)]
data[dataPast, on = .(zip, year),
     `:=` (PriorMortality = i.PriorMortality, PriorCOPD = i.PriorCOPD)]

## Clean up helper columns
data[, `:=` (YearPlus = NULL, YearMinus = NULL)]

## Store data frame without NAs
finalData = na.omit(data, cols = variablesOfInterest)
setDT(finalData)

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
data = dataOriginal

## Convert data to a data.table
setDT(data)
## Convert all columns to numeric (handles integers fine)
data[, (names(data)) := lapply(.SD, as.numeric)]

## Convert ZIP code information to a numeric bin for matching
data[, ZipBin := as.numeric(substr(zip, 1, nchar(as.character(zip)) - 4))]
data[is.na(ZipBin), ZipBin := 0] # Handle any NAs which might have been introduced

## New columns for distant exposure and outcome variables
data[, `:=` (
  DistantPM25 = NA_real_,
  DistantOzone = NA_real_,
  DistantMortality = NA_real_,
  DistantCOPD = NA_real_
)]

## Calculate Distant ZipBins
data[, DistantZipBin := fifelse(ZipBin >= 5, ZipBin - 5, ZipBin + 5)]

## Perform a self-join to find potential distant matches
distant_matches <- data[data, on = .(ZipBin = DistantZipBin, year),
                        .(zip, year, zip_distant = i.zip,
                          pm25_distant = i.pm25, ozone_distant = i.ozone,
                          death_rate_distant = i.death_rate,
                          copd_rate_distant = i.copd_rate),
                        allow.cartesian = TRUE, nomatch = 0L]

## Randomly select one distant zip for each (zip, year)
set.seed(325) # Set random seed for reproducibility 
distant_selected <- distant_matches[, .SD[sample(.N, 1)], by = .(zip, year)]

## Update data with the randomly selected distant values
data[distant_selected, on = .(zip, year), `:=` (
  DistantPM25 = pm25_distant,
  DistantOzone = ozone_distant,
  DistantMortality = death_rate_distant,
  DistantCOPD = copd_rate_distant
)]

## Clean up temporary columns if they are no longer needed
data[, `:=` (ZipBin = NULL, DistantZipBin = NULL)]

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
                        "DistantPM25",
                        "DistantOzone",
                        "DistantMortality",
                        "DistantCOPD")


## Store data frame without NAs
finalData = na.omit(data, cols = variablesOfInterest)
setDT(finalData)


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

##==========================================================================##
## Four ANALYSES
##==========================================================================##
source("MainFunctions.R")

## ANALYSES 0: Original setting
cat("ANALYSIS0","\n")
A0 = SettingList[["ANALYSIS0"]]
test0 = multiFunc(Y=A0$Y, 
                  Tr=A0$Tr, 
                  X=A0$X, 
                  b=A0$b, 
                  t1=A0$t1, 
                  t2=A0$t2, 
                  t1NC=A0$t1NC, 
                  t2NC=A0$t2NC,
                  maxM = A0$m,
                  scaleData = TRUE,
                  nB = 50)

## ANALYSIS 1: ANALYSIS 0 but No NCEs, only keeping NCOs
cat("ANALYSIS1","\n")
A1 = SettingList[["ANALYSIS1"]]
test1 = multiFunc(Y=A1$Y, 
                  Tr=A1$Tr, 
                  X=A1$X, 
                  b=A1$b, 
                  t1=A1$t1, 
                  t2=A1$t2, 
                  t1NC=A1$t1NC, 
                  t2NC=A1$t2NC,
                  maxM = A1$m,
                  scaleData = TRUE,
                  nB = 50)

## ANALYSIS 2: ANALYSIS 0 but 3 years ahead and back in time
cat("ANALYSIS2","\n")
A2 = SettingList[["ANALYSIS2"]]
test2 = multiFunc(Y=A2$Y, 
                  Tr=A2$Tr, 
                  X=A2$X, 
                  b=A2$b, 
                  t1=A2$t1, 
                  t2=A2$t2, 
                  t1NC=A2$t1NC, 
                  t2NC=A2$t2NC,
                  maxM = A2$m,
                  scaleData = TRUE,
                  nB = 50)

## ANALYSIS 3: Spatial NCs
cat("ANALYSIS3","\n")
A3 = SettingList[["ANALYSIS3"]]
test3 = multiFunc(Y=A3$Y, 
                  Tr=A3$Tr, 
                  X=A3$X, 
                  b=A3$b, 
                  t1=A3$t1, 
                  t2=A3$t2, 
                  t1NC=A3$t1NC, 
                  t2NC=A3$t2NC,
                  maxM = A3$m,
                  scaleData = TRUE,
                  nB = 50)

if(args$fake){
  save(test1, test2, test3, test4, file="data/output/fake_OutputSaved.dat")
} else {
  save(test1, test2, test3, test4, file="data/output/OutputSaved.dat")
}
