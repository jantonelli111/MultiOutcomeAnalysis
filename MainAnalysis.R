rm(list=ls())

## Read in main function needed for analysis
source("MainFunctions.R")

## Load in the necessary libraries
require(earth)
require(MASS)
require(psych)

## Read in fake data
dataOriginal = read.csv("Data/outcomes_merged_all_years.csv")

## Convert all variables to numeric (some are integers)
data = data.frame(sapply(dataOriginal, as.numeric))

## Need to create future exposures and prior outcomes variables for 
## Use as negative controls
data$FuturePM25 = rep(NA, nrow(data))
data$FutureOzone = rep(NA, nrow(data))

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
      
      data$FutureOzone[wCurrent] = data$ozone[wFuture]
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
                        "female_percentage",
                        "dual_percentage",
                        "mean_age",
                        "percentage_race_labelBlack",
                        "percentage_race_labelWhite",
                        "percentage_race_labelHispanic",
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
                 "female_percentage",
                 "dual_percentage",
                 "mean_age",
                 "percentage_race_labelBlack",
                 "percentage_race_labelWhite",
                 "percentage_race_labelHispanic")]

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
save(test, file="OutputSaved.dat")



