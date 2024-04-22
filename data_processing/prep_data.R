library(tidyverse)
library(arrow)

# read covars
covars_df_list <- list()
for(year in 2000:2016){
  covars_df_list[[year]] <- read_rds(paste0("data/input/covars/aggregate_data_", year, ".rds"))
}
covars_df <- bind_rows(covars_df_list)
covars_df$zip <- as.integer(covars_df$zip)

# read multioutcome
multioutcome_df <- read_parquet("data/input/MultiOutcomeAnalysis/multioutcome.parquet")

# identify set of zipcodes
covars_zips <- unique(covars_df$zip)
print(length(covars_zips))

outcome_cols <- c("alzh_rate", "alzhdmta_rate", "ami_rate", "anemia_rate", "asthma_rate",
          "atrialfb_rate", "breastCancer_rate", "cataract_rate", "chf_rate",
          "chrnkidn_rate", "colorectalCancer_rate", "copd_rate", "depressn_rate",
          "diabetes_rate", "endometrialCancer_rate", "glaucoma_rate", "hipfrac_rate",
          "hyperl_rate", "hyperp_rate", "hypert_rate", "hypoth_rate", "ischmcht_rate",
          "lungCancer_rate", "osteoprs_rate", "prostateCancer_rate", "ra_oa_rate",
          "stroke_rate")

na_count_vector <- as.integer(rowSums(is.na(multioutcome_df[outcome_cols])))
print(table(na_count_vector))
index <- which(na_count_vector == 0)
multioutcome_zips <- unique(multioutcome_df[index,]$zip)
print(length(multioutcome_zips))

zips <- intersect(covars_zips, multioutcome_zips)
print(length(zips))

# merge data
df <- covars_df %>% 
  inner_join(multioutcome_df[index,])

dim(df)
colnames(df)
# [1] "zip"                    "year"                   "pm25_ensemble"          "mean_bmi"              
# [5] "smoke_rate"             "hispanic"               "pct_blk"                "medhouseholdincome"    
# [9] "medianhousevalue"       "poverty"                "education"              "popdensity"            
# [13] "pct_owner_occ"          "summer_tmmx"            "winter_tmmx"            "summer_rmax"           
# [17] "winter_rmax"            "pm25"                   "ozone"                  "br"                    
# [21] "ca"                     "cu"                     "ec"                     "fe"                    
# [25] "k"                      "nh4"                    "ni"                     "no3"                   
# [29] "oc"                     "pb"                     "si"                     "so4"                   
# [33] "v"                      "z"                      "female_pct"             "dual_pct"              
# [37] "mean_age"               "race_white_pct"         "race_black_pct"         "race_hispanic_pct"     
# [41] "race_asian_pct"         "death_rate"             "alzh_rate"              "alzhdmta_rate"         
# [45] "ami_rate"               "anemia_rate"            "asthma_rate"            "atrialfb_rate"         
# [49] "breastCancer_rate"      "cataract_rate"          "chf_rate"               "chrnkidn_rate"         
# [53] "colorectalCancer_rate"  "copd_rate"              "depressn_rate"          "diabetes_rate"         
# [57] "endometrialCancer_rate" "glaucoma_rate"          "hipfrac_rate"           "hyperl_rate"           
# [61] "hyperp_rate"            "hypert_rate"            "hypoth_rate"            "ischmcht_rate"         
# [65] "lungCancer_rate"        "osteoprs_rate"          "prostateCancer_rate"    "ra_oa_rate"            
# [69] "stroke_rate"  

# write processed data
write_parquet(df, "data/processed/data.parquet")
