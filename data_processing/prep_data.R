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

# merge data
df <- covars_df %>% 
  left_join(multioutcome_df)

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
# [41] "race_asian_pct"         "mort_rate"              "alzh_rate"              "alzhdmta_rate"         
# [45] "ami_rate"               "anemia_rate"            "asthma_rate"            "atrialfb_rate"         
# [49] "breastCancer_rate"      "cataract_rate"          "chf_rate"               "chrnkidn_rate"         
# [53] "colorectalCancer_rate"  "copd_rate"              "depressn_rate"          "diabetes_rate"         
# [57] "endometrialCancer_rate" "glaucoma_rate"          "hipfrac_rate"           "hyperl_rate"           
# [61] "hyperp_rate"            "hypert_rate"            "hypoth_rate"            "ischmcht_rate"         
# [65] "lungCancer_rate"        "osteoprs_rate"          "prostateCancer_rate"    "ra_oa_rate"            
# [69] "stroke_rate"

# write processed data
write_parquet(df, "data/processed/outcomes_merged_all_years.parquet")

# compare columns with fake data
fake_df <- read_csv("data/fake_outcomes_merged_all_years.csv")
dim(fake_df)
colnames(fake_df)
# [1] "zip"                                        "year"                                       "pm25_ensemble"                              "mean_bmi"                                  
# [5] "smoke_rate"                                 "hispanic"                                   "pct_blk"                                    "medhouseholdincome"                        
# [9] "medianhousevalue"                           "poverty"                                    "education"                                  "pct_owner_occ"                             
# [13] "summer_tmmx"                                "winter_tmmx"                                "summer_rmax"                                "winter_rmax"                               
# [17] "pm25"                                       "ozone"                                      "br"                                         "ca"                                        
# [21] "cu"                                         "ec"                                         "fe"                                         "k"                                         
# [25] "nh4"                                        "ni"                                         "no3"                                        "oc"                                        
# [29] "pb"                                         "si"                                         "so4"                                        "v"                                         
# [33] "z"                                          "death_rate"                                 "female_percentage"                          "dual_percentage"                           
# [37] "mean_age"                                   "percentage_race_labelAsian"                 "percentage_race_labelBlack"                 "percentage_race_labelOther"                
# [41] "percentage_race_labelWhite"                 "percentage_race_labelHispanic"              "percentage_race_labelNorth_American_Native" "anemia_rate"                               
# [45] "ami_rate"                                   "alzh_rate"                                  "alzhdmta_rate"                              "atrialfb_rate"                             
# [49] "cataract_rate"                              "chrnkidn_rate"                              "copd_rate"                                  "chf_rate"                                  
# [53] "diabetes_rate"                              "stroke_rate"                                "breastCancer_rate"                          "colorectalCancer_rate"                     
# [57] "prostateCancer_rate"                        "lungCancer_rate"                            "endometrialCancer_rate"                     "hyperp_rate"                               
# [61] "glaucoma_rate"                              "hipfrac_rate"                               "ischmcht_rate"                              "depressn_rate"                             
# [65] "osteoprs_rate"                              "ra_oa_rate"                                 "asthma_rate"                                "hyperl_rate"                               
# [69] "hypert_rate"                                "hypoth_rate"  