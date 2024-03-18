library(arrow)

rm(list=ls())

## Read in main function needed for analysis
source("MainFunctions.R")

## Load in the necessary libraries
require(earth)
require(MASS)
require(psych)

data <- read_parquet("data/processed/data.parquet")

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
                        "PriorMortality",
                        "PriorCOPD")


## Remove NAs
wNotNA = which(complete.cases(data[,variablesOfInterest]) == TRUE)

## Store data frame without NAs
finalData = data[wNotNA,]

## Now create variables needed for our analysis
Y = finalData[,c("mort_rate",
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
                 "race_white_pct",
                 "race_black_pct",
                 "race_hispanic_pct")]

t1 = c(apply(Tr[,1:7], 2, quantile, 0.75), apply(Tr[,8:9], 2, median))
t2 = c(apply(Tr[,1:7], 2, quantile, 0.25), apply(Tr[,8:9], 2, median))

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

## Find maximum number of factors allowed
maxMfunc = function(qp, m){
  (qp - m)^2 - qp - m
}

for(pq in 3:min(ncol(Y), ncol(Tr))){
  m = 0
  while(maxMfunc(pq, m) >= 0){m = m + 1}
  m = m - 1
}

## Run main function on cleaned data
test = multiFunc(Y=Y, 
                 Tr=Tr, 
                 X=X, 
                 b=b, 
                 t1=t1, 
                 t2=t2, 
                 t1NC=t1NC, 
                 t2NC=t2NC,
                 maxM = m)

## Save results
save(test, file="data/output/OutputSaved.dat")



