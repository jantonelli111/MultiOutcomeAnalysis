library(arrow)

## Read data
dataOriginal = read_parquet("data/processed/outcomes_merged_all_years.parquet")

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
      
      data$PriorMortality[wCurrent] = data$mort_rate[wPast]
      data$PriorCOPD[wCurrent] = data$copd_rate[wPast]
      
      data$FutureOzone[wCurrent] = data$ozone[wFuture]
      data$FuturePM25[wCurrent] = data$pm25[wFuture]
    }
  }
}

write_parquet(data, "data/processed/data.parquet")