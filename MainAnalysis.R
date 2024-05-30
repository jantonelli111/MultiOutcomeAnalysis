
rm(list=ls())

# define parser arguments ----
# args <- list()
# args$fake <- TRUE
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--fake", action = "store_true", help = "use fake data")  
parser$add_argument("--analysis", default = "test1", help = "analysis to run")             
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
newCols <- c("PriorCOPD")
data[, (newCols) := .(rep(NA_real_, .N))]

## Create shifts in year columns for past  matching
data[, c("YearPlus1") := .(year + 1)]

## Create a key for faster subsetting and joins
setkey(data, zip, year)

## Join data with itself to fetch future and past values
dataCurrent <- data[, .(zip, year = year,
                        copd_rate = copd_rate)]
dataPast <- data[, .(zip, year = YearPlus1,
                     PriorCOPD = copd_rate)]

## Merge \past data into the original dataset by reference\
data[dataPast, on = .(zip, year),
     `:=` (PriorCOPD = i.PriorCOPD)]

## Clean up helper columns
data[, `:=` (YearPlus1 = NULL)]

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
                        "PriorCOPD")


## Store data frame without NAs
finalData = na.omit(data, cols = variablesOfInterest)
setDT(finalData)



##----------------------------------------------------------------------------##
## ANALYSIS 0: No NCEs, No NCOs
##----------------------------------------------------------------------------##
## Now create variables needed for our analysis
Y = finalData[,c("death_rate",
                 "anemia_rate",
                 "copd_rate",
                 "stroke_rate",
                 "lungCancer_rate",
                 "asthma_rate",
                 "hypert_rate",
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

## Scale data
Tr[] = lapply(Tr, scale, center = TRUE, scale = TRUE)
Y[] = lapply(Y, scale, center = TRUE, scale = TRUE)

## estimand of interest
t1 = list()
t2 = list()
t1[[1]] = apply(Tr[,1:7], 2, quantile, 0.75)
t2[[1]] = apply(Tr[,1:7], 2, quantile, 0.25)

for (tt in 2 : 8) {
  t1[[tt]] = apply(Tr, 2, median)
  t1[[tt]][tt-1] = apply(Tr[, .SD, .SDcols = c(tt-1)], 2, quantile, 0.75)
  t2[[tt]] = apply(Tr, 2, median)
  t2[[tt]][tt-1] = apply(Tr[, .SD, .SDcols = c(tt-1)], 2, quantile, 0.25)
}

## No negative controls
b = t1NC = t2NC = NULL

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
## ANALYSIS 1: No NCEs, One NCO
##----------------------------------------------------------------------------##
## Now create variables needed for our analysis
Y = finalData[,c("death_rate",
                 "anemia_rate",
                 "copd_rate",
                 "stroke_rate",
                 "lungCancer_rate",
                 "asthma_rate",
                 "hypert_rate",
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

## Scale data
Tr[] = lapply(Tr, scale, center = TRUE, scale = TRUE)
Y[] = lapply(Y, scale, center = TRUE, scale = TRUE)

## estimand of interest
t1 = list()
t2 = list()
t1[[1]] = apply(Tr[,1:7], 2, quantile, 0.75)
t2[[1]] = apply(Tr[,1:7], 2, quantile, 0.25)

for (tt in 2 : 8) {
  t1[[tt]] = apply(Tr, 2, median)
  t1[[tt]][tt-1] = apply(Tr[, .SD, .SDcols = c(tt-1)], 2, quantile, 0.75)
  t2[[tt]] = apply(Tr, 2, median)
  t2[[tt]][tt-1] = apply(Tr[, .SD, .SDcols = c(tt-1)], 2, quantile, 0.25)
}

## Negative control outcomes
b = matrix(c(0, 0, 0, 0, 0, 0, 0, 1),
           byrow = FALSE, nrow = ncol(Y))

## Exposure contrasts for negative controls
t1NC = list()
t2NC = list()

## Now use current exposures with prior outcomes
t1NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))
t2NC[[1]] = t(matrix(rep(apply(Tr, 2, median), 7), ncol=7))

for (tt in 1 : 7) {
  t1NC[[1]][tt,tt] = apply(Tr[, .SD, .SDcols = c(tt)], 2, quantile, 0.75)
  t2NC[[1]][tt,tt] = apply(Tr[, .SD, .SDcols = c(tt)], 2, quantile, 0.25)
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


##==========================================================================##
## TWO ANALYSES
##==========================================================================##
source("MainFunctions.R")



if(args$analysis == "test0"){
  ## ANALYSES 0: Original setting
  cat("ANALYSIS0","\n")
  #A0 = SettingList[["ANALYSIS0"]]
  test = multiFunc(Y=A0$Y, 
                   Tr=A0$Tr, 
                   X=A0$X, 
                   b=A0$b, 
                   t1=A0$t1, 
                   t2=A0$t2, 
                   t1NC=A0$t1NC, 
                   t2NC=A0$t2NC,
                   maxM = A0$m,
                   nB = 50)
  
} else if(args$analysis == "test1"){
  ## ANALYSIS 1: ANALYSIS 0 but No NCEs, only keeping NCOs
  cat("ANALYSIS1","\n")
  #A1 = SettingList[["ANALYSIS1"]]
  test = multiFunc(Y=A1$Y, 
                   Tr=A1$Tr, 
                   X=A1$X, 
                   b=A1$b, 
                   t1=A1$t1, 
                   t2=A1$t2, 
                   t1NC=A1$t1NC, 
                   t2NC=A1$t2NC,
                   maxM = A1$m,
                   nB = 50)
  
}

if(args$fake){
  output_file = paste0("data/output/fake_output_", args$analysis, ".dat")
} else {
  output_file = paste0("data/output/output_", args$analysis, ".dat")
}

save(test, file=output_file)
print(paste0("Output saved to: ", output_file))
