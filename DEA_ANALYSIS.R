library(tidyverse)
library(bnlearn)
library(mice)
library(survival)
library(mediation)
library(randomForest)
library(dplyr)

setwd("/Users/alikmofrad/UCLA PhD/Research/Kidney_Exchange")

# Min-Max Scaling
scale_min_max <- function(x) {
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  return(x_scaled)
}

# Load and preprocess data
kidney_data <-  read.csv("KIDPAN_DECEASED.csv")

kidney_data <- kidney_data %>% filter(ORGAN == "KI")

kidney_data$DAYSWAIT_CHRON_KI <- as.numeric(kidney_data$DAYSWAIT_CHRON_KI)

kidney_data <- kidney_data %>% 
  mutate(ETHCAT = case_when(
    ETHCAT == 1 ~ "White",
    ETHCAT == 2 ~ "Black",
    ETHCAT == 4 ~ "Hispanic",
    ETHCAT == 5 ~ "Asian",
    ETHCAT == 6 ~ "Other",
    ETHCAT == 7 ~ "Other",
    ETHCAT == 9 ~ "Other",
    ETHCAT == 998 ~ "Other",
    TRUE ~ as.character(ETHCAT) 
  ),
  ETHCAT_DON = case_when(
    ETHCAT_DON == 1 ~ "White",
    ETHCAT_DON == 2 ~ "Black",
    ETHCAT_DON == 4 ~ "Hispanic",
    ETHCAT_DON == 5 ~ "Asian",
    ETHCAT_DON == 6 ~ "Other",
    ETHCAT_DON == 7 ~ "Other",
    ETHCAT_DON == 9 ~ "Other",
    ETHCAT_DON == 998 ~ "Other",
    TRUE ~ as.character(ETHCAT_DON) 
  )
  )

kidney_data$ETHCAT <- as.factor(kidney_data$ETHCAT)
kidney_data$ETHCAT_DON <- as.factor(kidney_data$ETHCAT_DON)

kidney_data$KDPI <-(sapply(kidney_data$KDPI, function(x) {
  if (is.na(x) || x == ".") {
    NA
  } else {
    as.numeric(sub("%", "", x)) / 100
  }
}))

kidney_data$GTIME_KI <- as.numeric(kidney_data$GTIME_KI)
kidney_data$ON_DIALYSIS <- as.factor(kidney_data$ON_DIALYSIS)
kidney_data$ABO <- as.factor(kidney_data$ABO)
kidney_data$INIT_CPRA <- as.numeric(kidney_data$INIT_CPRA)
kidney_data$AGE <- as.numeric(kidney_data$AGE)
kidney_data$GENDER <- as.factor(kidney_data$GENDER)
kidney_data$PRA <- as.numeric(kidney_data$END_CPRA)
kidney_data$WGT_KG_CALC <- as.numeric(kidney_data$WGT_KG_CALC)
kidney_data$BMI <- as.numeric(kidney_data$BMI_CALC)
# 
# kidney_data_temp <- kidney_data
# kidney_data <- kidney_data_temp

kidney_data <- kidney_data[kidney_data$ON_DIALYSIS != "1",]




kidney_data <- kidney_data %>% filter(ETHCAT != "Other", ETHCAT_DON != "Other")
kidney_data$ETHCAT <- droplevels(kidney_data$ETHCAT)
kidney_data$ETHCAT_DON <- droplevels((kidney_data$ETHCAT_DON))
# kidney_data <- data.frame(Ethnicity = kidney_data$ETHCAT, Waitlist_Duration = kidney_data$DAYSWAIT_CHRON_KI, Graft_Lifespan = kidney_data$GTIME_KI, Dialysis_Status = kidney_data$ON_DIALYSIS, ABO_Blood_Type = kidney_data$ABO, PRA_Score = kidney_data$INIT_CPRA, KDPI = kidney_data$KDPI)

# kidney_data_temp <- kidney_data
#kidney_data <- kidney_data_temp


#selected_vars <- c("DAYSWAIT_CHRON_KI", "GTIME_KI", "ON_DIALYSIS", "ABO", "INIT_CPRA", "AGE", "GENDER", "PRA", "KDPI", "ETHCAT", "ETHCAT_DON", "WGT_KG_CALC", "BMI")


#kidney_data <- kidney_data[, selected_vars]
#kidney_data <- kidney_data[complete.cases(kidney_data), ]

#kidney_data <- kidney_data[sample(1:nrow(kidney_data), size = 50000, replace = F),]
# library(MatchIt)
kidney_data <- kidney_data %>% mutate(PRI_PAYMENT_TRR_KI = case_when(
  PRI_PAYMENT_TRR_KI %in% 1:14 ~ as.factor(PRI_PAYMENT_TRR_KI),
  !(PRI_PAYMENT_TRR_KI %in% 1:14) ~ NA),
  
  ACUTE_REJ_EPI_KI = case_when(
    ACUTE_REJ_EPI_KI %in% 1:3 ~ as.factor(ACUTE_REJ_EPI_KI),
  !(ACUTE_REJ_EPI_KI %in% 1:3) ~ NA),
  GRF_FAIL_CAUSE_TY_KI = case_when(
    GRF_FAIL_CAUSE_TY_KI %in% c(1:12, 999) ~ as.factor(GRF_FAIL_CAUSE_TY_KI),
    !(GRF_FAIL_CAUSE_TY_KI %in% c(1:12, 999)) ~ NA),
  EDUCATION = case_when(
    EDUCATION %in% c(1:6, 996, 998) ~ as.factor(EDUCATION),
    !(EDUCATION %in% c(1:6, 996, 998)) ~ NA),
  URINE_INF_DON = case_when(
    URINE_INF_DON %in% c(0,1) ~ as.factor(URINE_INF_DON),
    !(URINE_INF_DON %in% c(0,1)) ~ NA)
  )


kidney_data$PRI_PAYMENT_TRR_KI <- droplevels(kidney_data$PRI_PAYMENT_TRR_KI)
kidney_data$ACUTE_REJ_EPI_KI <- droplevels(kidney_data$ACUTE_REJ_EPI_KI)
kidney_data$EDUCATION <- droplevels(kidney_data$EDUCATION)
kidney_data$URINE_INF_DON <- droplevels(kidney_data$URINE_INF_DON)




kidney_data$REJTRT_KI = as.factor(kidney_data$REJTRT_KI)
kidney_data$CREAT_DON = as.numeric(kidney_data$CREAT_DON)
kidney_data$AGE_DON = as.numeric(kidney_data$AGE_DON)
kidney_data$HGT_CM_DON_CALC = as.numeric(kidney_data$HGT_CM_DON_CALC)
kidney_data$WGT_KG_DON_CALC = as.numeric(kidney_data$WGT_KG_DON_CALC)
kidney_data$HIST_HYPERTENS_DON = as.factor(kidney_data$HIST_HYPERTENS_DON)
kidney_data$HYPERTENS_DUR_DON = as.factor(kidney_data$HYPERTENS_DUR_DON)
kidney_data$GFR = as.numeric(kidney_data$GFR)
kidney_data$DISTANCE = as.numeric(kidney_data$DISTANCE)
kidney_data$PREV_KI_TX = as.factor(kidney_data$PREV_KI_TX)
kidney_data$HIST_DIABETES_DON = as.factor(kidney_data$HIST_DIABETES_DON)
kidney_data$HIST_CANCER_DON = as.factor(kidney_data$HIST_CANCER_DON)
kidney_data$HIST_CIG_DON = as.factor(kidney_data$HIST_CIG_DON)
kidney_data$HIST_COCAINE_DON = as.factor(kidney_data$HIST_COCAINE_DON)
kidney_data$COD_CAD_DON = as.factor(kidney_data$COD_CAD_DON)
kidney_data$DIAG_KI = as.factor(kidney_data$DIAG_KI)
kidney_data$CREAT_TRR = as.numeric(kidney_data$CREAT_TRR)
kidney_data$ACADEMIC_PRG_TCR = as.factor(kidney_data$ACADEMIC_PRG_TCR)
kidney_data$ACADEMIC_LEVEL_TCR = as.factor(kidney_data$ACADEMIC_LEVEL_TCR)
kidney_data$HLAMIS = as.numeric(kidney_data$HLAMIS)
kidney_data$SERUM_CREAT = as.numeric(kidney_data$SERUM_CREAT)
kidney_data$DIAB = as.factor(kidney_data$DIAB)
kidney_data$EDUCATION = as.factor(kidney_data$EDUCATION)
kidney_data$ABO_MAT = as.factor(kidney_data$ABO_MAT)
kidney_data$PREV_TX_ANY = as.factor(kidney_data$PREV_TX_ANY)
kidney_data$URINE_INF_DON = as.factor(kidney_data$URINE_INF_DON)
kidney_data$CDC_RISK_HIV_DON = as.factor(kidney_data$CDC_RISK_HIV_DON)
kidney_data$DRUGTRT_COPD = as.factor(kidney_data$DRUGTRT_COPD)
kidney_data$REGION = as.factor(kidney_data$REGION)
kidney_data$BUN_DON = as.numeric(kidney_data$BUN_DON)
kidney_data$TOT_SERUM_ALBUM = as.numeric(kidney_data$TOT_SERUM_ALBUM)
kidney_data$WORK_INCOME_TCR = as.factor(kidney_data$WORK_INCOME_TCR)
kidney_data$MED_COND_TRR = as.factor(kidney_data$MED_COND_TRR)
kidney_data$LOS = as.numeric(kidney_data$LOS)
kidney_data$TRTREJ1Y_KI = as.factor(kidney_data$TRTREJ1Y_KI)
kidney_data$TRTREJ6M_KI = as.factor(kidney_data$TRTREJ6M_KI)
kidney_data$TX_PROCEDUR_TY_KI = as.factor(kidney_data$TX_PROCEDUR_TY_KI)
kidney_data$INIT_STAT = as.factor(kidney_data$INIT_STAT)
kidney_data$END_STAT = as.factor(kidney_data$END_STAT)
kidney_data$COLD_ISCH_KI = as.numeric(kidney_data$COLD_ISCH_KI)
kidney_data$HGT_CM_CALC = as.numeric(kidney_data$HGT_CM_CALC)

# ############## Imputations (doesnt seem to help) ##############
# library(mice)
# 
# # Define the features with missing values you want to impute
# features_to_impute <- c("TRTREJ6M_KI", "TRTREJ1Y_KI", "LOS", "WORK_INCOME_TCR", "TOT_SERUM_ALBUM", "BUN_DON", "DRUGTRT_COPD", "CDC_RISK_HIV_DON", "EDUCATION",
#                         "SERUM_CREAT", "HLAMIS", "ACUTE_REJ_EPI_KI", "CREAT_TRR", "ACUTE_REJ_EPI_KI", "REJTRT_KI", "CREAT_DON", "HIST_CIG_DON", "HIST_CANCER_DON", "HIST_HYPERTENS_DON", 
#                         "WGT_KG_DON_CALC", "HGT_CM_DON_CALC")
# # Perform imputation on the dataset
# set.seed(14)  # for reproducibility
# mice_data <- mice(kidney_data[, features_to_impute], m = 5, method = 'pmm', maxit = 5)
# 
# # Obtain the completed data
# completed_data <- complete(mice_data, 1)  
# completed_data_pooled <- pool(mice_data) 
# #################################################



dmu_data <- data.frame(
  Group = kidney_data$ETHCAT,
  
  WaitlistDuration = kidney_data$DAYSWAIT_CHRON_KI, # Input 1
  QualityScore = kidney_data$KDPI, # Input 2
  OutcomeScore = kidney_data$GTIME_KI,# Output 1
  
  # Confounders:
  BMI = kidney_data$BMI, # 8083 NA
  WGT_KG_CALC = kidney_data$WGT_KG_CALC, # 4113 NA
  ABO = kidney_data$ABO, # NO NA
  ON_DIALYSIS = kidney_data$ON_DIALYSIS,# NO NA
  AGE = kidney_data$AGE, # NO NA
  AGE_DON = kidney_data$AGE_DON, # NO NA
  HGT_CM_DON_CALC =  kidney_data$HGT_CM_DON_CALC, #29118 NA
  WGT_KG_DON_CALC = kidney_data$WGT_KG_DON_CALC, # 15931 NA 
  ETHCAT_DON = kidney_data$ETHCAT_DON, # NO NA
  HIST_HYPERTENS_DON = kidney_data$HIST_HYPERTENS_DON, # 44311 NA
  HYPERTENS_DUR_DON = kidney_data$HYPERTENS_DUR_DON, # NO NA
  PRA = kidney_data$PRA, # 138370 NA
  GENDER = kidney_data$GENDER, # NO NA
  GFR = kidney_data$GFR, # 238971 NA
  DISTANCE = kidney_data$DISTANCE, # 3295 NA
  PRI_PAYMENT_TRR_KI = kidney_data$PRI_PAYMENT_TRR_KI, # 45354 NA
  PREV_KI_TX = kidney_data$PREV_KI_TX, # NO NA
  
  HIST_DIABETES_DON = kidney_data$HIST_DIABETES_DON, # NO NA
  HIST_CANCER_DON = kidney_data$HIST_CANCER_DON, # 44452 NA
  HIST_CIG_DON = kidney_data$HIST_CIG_DON, # 44311 NA
  HIST_COCAINE_DON = kidney_data$HIST_COCAINE_DON, # 600 NA

  CREAT_DON = kidney_data$CREAT_DON, # 45188 NA
  COD_CAD_DON = kidney_data$COD_CAD_DON, # NO NA
  
  DIAG_KI = kidney_data$DIAG_KI, # NO NA
  
  #.     REJTRT_KI = kidney_data$REJTRT_KI, # USEFUL BUT CAUSES LOT OF NA (223446 NA)
  
  #GRF_FAIL_CAUSE_TY_KI = as.factor(kidney_data$GRF_FAIL_CAUSE_TY_KI),
  
  ACUTE_REJ_EPI_KI = kidney_data$ACUTE_REJ_EPI_KI, # 103385 NA
  CREAT_TRR = kidney_data$CREAT_TRR, # 66847 NA
  
  # Waitlist relevant Confounders
  ACADEMIC_PRG_TCR = kidney_data$ACADEMIC_PRG_TCR, # NO NA
  ACADEMIC_LEVEL_TCR = kidney_data$ACADEMIC_LEVEL_TCR, # NO NA
  HLAMIS = kidney_data$HLAMIS, # 2231 NA
  
  SERUM_CREAT = kidney_data$SERUM_CREAT, # 8507 NA
  
  DIAB = kidney_data$DIAB, # NO NA
  
  EDUCATION = kidney_data$EDUCATION, # 50314 NA
  ABO_MAT = kidney_data$ABO_MAT, # 3 NA
  
  COLD_ISCH_KI = kidney_data$COLD_ISCH_KI, # 14684 NA
  
  PREV_TX_ANY = kidney_data$PREV_TX_ANY, # NO NA
  
  URINE_INF_DON = kidney_data$URINE_INF_DON, # 2 NA
  
  CDC_RISK_HIV_DON = kidney_data$CDC_RISK_HIV_DON, # 102812 NA
  
  DRUGTRT_COPD = kidney_data$DRUGTRT_COPD, # 129654 NA
  
  REGION = kidney_data$REGION, # NO NA
  
  BUN_DON = kidney_data$BUN_DON, # 45209 NA
  # LIPASE = as.numeric(kidney_data$LIPASE) # too much NA
  # AMYLASE = as.numeric(kidney_data$AMYLASE) # too much NA
  TOT_SERUM_ALBUM = kidney_data$TOT_SERUM_ALBUM, # 105676 NA
  
  WORK_INCOME_TCR = kidney_data$WORK_INCOME_TCR, # 130557 NA
  
  MED_COND_TRR = kidney_data$MED_COND_TRR, # 0 NA
  
  #FUNC_STAT_TRR = as.factor(kidney_data$FUNC_STAT_TRR) # 0 NA
  
  LOS = kidney_data$LOS,
  
  TRTREJ1Y_KI = kidney_data$TRTREJ1Y_KI,
  
  TRTREJ6M_KI = kidney_data$TRTREJ6M_KI,
  
  TX_PROCEDUR_TY_KI = kidney_data$TX_PROCEDUR_TY_KI,
  
  INIT_STAT = kidney_data$INIT_STAT, #  marginal improvement
  END_STAT = kidney_data$END_STAT, #  marginal improvement
  
  NPKID = kidney_data$NPKID,
  
  HGT_CM_CALC = kidney_data$HGT_CM_CALC
)





dim(dmu_data[complete.cases(dmu_data),])
dmu_data <- dmu_data[complete.cases(dmu_data),]


# library(randomForest)
# 
# # Prepare the data
# confounders <- c("BMI", "WGT_KG_CALC", "ABO", "ON_DIALYSIS", "AGE", "AGE_DON", 
#                  "HGT_CM_DON_CALC", "WGT_KG_DON_CALC", "ETHCAT_DON", 
#                  "HIST_HYPERTENS_DON", "HYPERTENS_DUR_DON", "PRA", "GENDER", 
#                  "GFR", "DISTANCE", "PRI_PAYMENT_TRR_KI", "PREV_KI_TX", 
#                  "HIST_DIABETES_DON", "HIST_CANCER_DON", "HIST_CIG_DON", 
#                  "HIST_COCAINE_DON", "CREAT_DON", "COD_CAD_DON", "DIAG_KI", 
#                  "ACUTE_REJ_EPI_KI", "CREAT_TRR", "ACADEMIC_PRG_TCR", 
#                  "ACADEMIC_LEVEL_TCR", "HLAMIS", "SERUM_CREAT", "DIAB",
#                  "EDUCATION", "ABO_MAT", "PREV_TX_ANY", "URINE_INF_DON", 
#                  "CDC_RISK_HIV_DON", "DRUGTRT_COPD", "REGION", "BUN_DON", 
#                  "TOT_SERUM_ALBUM","WORK_INCOME_TCR", "MED_COND_TRR", "LOS",
#                  "TRTREJ1Y_KI", "TRTREJ6M_KI", "TX_PROCEDUR_TY_KI",
#                  "INIT_STAT", "END_STAT", "NPKID", "COLD_ISCH_KI")


 # KDPI_confounders <- c("REJTRT_KI")
# dmu_data <- dmu_data[ , c("Group", "WaitlistDuration", "QualityScore", "OutcomeScore", KDPI_confounders)]



# dmu_data <- dmu_data[complete.cases(dmu_data), c("Group", "WaitlistDuration", "QualityScore", "OutcomeScore", confounders)]



## Setting the data subsample ##
# set.seed(14)
# # sample_index <- sample(1:nrow(dmu_data), replace = F, size = 10000)
# # dmu_data <- dmu_data[sample_index,]
# dmu_data <- data.frame(
#   Group = kidney_data$ETHCAT,
#   
#   WaitlistDuration = kidney_data$DAYSWAIT_CHRON_KI, # Input 1
#   QualityScore = kidney_data$KDPI, # Input 2
#   OutcomeScore = kidney_data$GTIME_KI,# Output 1
#   
#   # Confounders:
#   BMI = kidney_data$BMI, # 8083 NA
#   WGT_KG_CALC = kidney_data$WGT_KG_CALC, # 4113 NA
#   ABO = kidney_data$ABO, # NO NA
#   ON_DIALYSIS = kidney_data$ON_DIALYSIS,# NO NA
#   AGE = kidney_data$AGE, # NO NA
#   AGE_DON = kidney_data$AGE_DON, # NO NA
#   HGT_CM_DON_CALC =  kidney_data$HGT_CM_DON_CALC, #29118 NA
#   WGT_KG_DON_CALC = kidney_data$WGT_KG_DON_CALC, # 15931 NA 
#   ETHCAT_DON = kidney_data$ETHCAT_DON, # NO NA
#   HIST_HYPERTENS_DON = kidney_data$HIST_HYPERTENS_DON, # 44311 NA
#   HYPERTENS_DUR_DON = kidney_data$HYPERTENS_DUR_DON, # NO NA
#   PRA = kidney_data$PRA, # 138370 NA
#   GENDER = kidney_data$GENDER, # NO NA
#   GFR = kidney_data$GFR, # 238971 NA
#   DISTANCE = kidney_data$DISTANCE, # 3295 NA
#   PRI_PAYMENT_TRR_KI = kidney_data$PRI_PAYMENT_TRR_KI, # 45354 NA
#   PREV_KI_TX = kidney_data$PREV_KI_TX, # NO NA
#   
#   HIST_DIABETES_DON = kidney_data$HIST_DIABETES_DON, # NO NA
#   HIST_CANCER_DON = kidney_data$HIST_CANCER_DON, # 44452 NA
#   HIST_CIG_DON = kidney_data$HIST_CIG_DON, # 44311 NA
#   HIST_COCAINE_DON = kidney_data$HIST_COCAINE_DON, # 600 NA
#   
#   CREAT_DON = kidney_data$CREAT_DON, # 45188 NA
#   COD_CAD_DON = kidney_data$COD_CAD_DON, # NO NA
#   
#   DIAG_KI = kidney_data$DIAG_KI, # NO NA
#   
#   #.     REJTRT_KI = kidney_data$REJTRT_KI, # USEFUL BUT CAUSES LOT OF NA (223446 NA)
#   
#   #GRF_FAIL_CAUSE_TY_KI = as.factor(kidney_data$GRF_FAIL_CAUSE_TY_KI),
#   
#   ACUTE_REJ_EPI_KI = kidney_data$ACUTE_REJ_EPI_KI, # 103385 NA
#   CREAT_TRR = kidney_data$CREAT_TRR, # 66847 NA
#   
#   # Waitlist relevant Confounders
#   ACADEMIC_PRG_TCR = kidney_data$ACADEMIC_PRG_TCR, # NO NA
#   ACADEMIC_LEVEL_TCR = kidney_data$ACADEMIC_LEVEL_TCR, # NO NA
#   HLAMIS = kidney_data$HLAMIS, # 2231 NA
#   
#   SERUM_CREAT = kidney_data$SERUM_CREAT, # 8507 NA
#   
#   DIAB = kidney_data$DIAB, # NO NA
#   
#   EDUCATION = kidney_data$EDUCATION, # 50314 NA
#   ABO_MAT = kidney_data$ABO_MAT, # 3 NA
#   
#   COLD_ISCH_KI = kidney_data$COLD_ISCH_KI, # 14684 NA
#   
#   PREV_TX_ANY = kidney_data$PREV_TX_ANY, # NO NA
#   
#   URINE_INF_DON = kidney_data$URINE_INF_DON, # 2 NA
#   
#   CDC_RISK_HIV_DON = kidney_data$CDC_RISK_HIV_DON, # 102812 NA
#   
#   DRUGTRT_COPD = kidney_data$DRUGTRT_COPD, # 129654 NA
#   
#   REGION = kidney_data$REGION, # NO NA
#   
#   BUN_DON = kidney_data$BUN_DON, # 45209 NA
#   # LIPASE = as.numeric(kidney_data$LIPASE) # too much NA
#   # AMYLASE = as.numeric(kidney_data$AMYLASE) # too much NA
#   TOT_SERUM_ALBUM = kidney_data$TOT_SERUM_ALBUM, # 105676 NA
#   
#   WORK_INCOME_TCR = kidney_data$WORK_INCOME_TCR, # 130557 NA
#   
#   MED_COND_TRR = kidney_data$MED_COND_TRR, # 0 NA
#   
#   #FUNC_STAT_TRR = as.factor(kidney_data$FUNC_STAT_TRR) # 0 NA
#   
#   LOS = kidney_data$LOS,
#   
#   TRTREJ1Y_KI = kidney_data$TRTREJ1Y_KI,
#   
#   TRTREJ6M_KI = kidney_data$TRTREJ6M_KI,
#   
#   TX_PROCEDUR_TY_KI = kidney_data$TX_PROCEDUR_TY_KI,
#   
#   INIT_STAT = kidney_data$INIT_STAT, #  marginal improvement
#   END_STAT = kidney_data$END_STAT, #  marginal improvement
# 
#   HGT_CM_CALC = kidney_data$HGT_CM_CALC
# )


####
total_samples <- nrow(dmu_data)  # or set a specific number if you prefer
desired_props <- c(Asian = 0.064, Black = 0.139, Hispanic = 0.195, White = 0.602)
desired_counts <- round(desired_props * total_samples)

balanced_data <- data.frame()

# Loop over each group and sample accordingly
set.seed(14)
for (group in names(desired_props)) {
  group_data <- dmu_data[dmu_data$Group == group, ]
  
  if (nrow(group_data) < desired_counts[group]) {
    # Oversample if there are fewer rows than needed
    sampled_data <- group_data[sample(nrow(group_data), desired_counts[group], replace = TRUE), ]
  } else {
    # Undersample if there are more rows than needed
    sampled_data <- group_data[sample(nrow(group_data), desired_counts[group], replace = FALSE), ]
  }
  
  # Bind the sampled data
  balanced_data <- rbind(balanced_data, sampled_data)
}

set.seed(12)  # for reproducibility
balanced_data <- balanced_data[sample(nrow(balanced_data)), ]
library(mosaic)
tally(balanced_data$Group, "proportion")
####




rm(desired_counts)
rm(desired_props)
#rm(selected_vars)



# # Train random forest models for treatment and inputs/outputs
# treatment_model <- randomForest(Group ~ ., data = balanced_data[, c("Group", confounders)])
# input1_model <- randomForest(WaitlistDuration ~ ., data = balanced_data[, c("WaitlistDuration", confounders)])
# input2_model <- randomForest(QualityScore ~ ., data = balanced_data[, c("QualityScore", confounders)])
# output_model <- randomForest(OutcomeScore ~ ., data = balanced_data[, c("OutcomeScore", confounders)])
# 
# ### taking subset of test_data cause its so big
# #set.seed(14)
# #test_data <- test_data[sample(1:nrow(test_data), replace = F, size = 10000),]
# #test_data <- dmu_data[-train_indices, ]
# 
# 
# # Estimate the treatment effect
# treatment_pred <- predict(treatment_model, newdata = test_data[, confounders])
# input1_pred <- predict(input1_model, newdata = test_data[, confounders])
# input2_pred <- predict(input2_model, newdata = test_data[, confounders])
# output_pred <- predict(output_model, newdata = test_data[, confounders])
# 
# treatment_effect <- data.frame(
#   Group = test_data$Group,
#   WaitlistDuration = (test_data$WaitlistDuration - input1_pred),
#   QualityScore = (test_data$QualityScore - input2_pred),
#   OutcomeScore = (test_data$OutcomeScore - output_pred)
# )
###############################################################################################
########################################## Mediation ##########################################
###############################################################################################
# c("BMI", "WGT_KG_CALC", "ABO", "ON_DIALYSIS", "AGE", "AGE_DON", 
#   "HGT_CM_DON_CALC", "WGT_KG_DON_CALC", "ETHCAT_DON", 
#   "HIST_HYPERTENS_DON", "HYPERTENS_DUR_DON", "PRA", "GENDER", 
#   "GFR", "DISTANCE", "PRI_PAYMENT_TRR_KI", "PREV_KI_TX", 
#   "HIST_DIABETES_DON", "HIST_CANCER_DON", "HIST_CIG_DON", 
#   "HIST_COCAINE_DON", "CREAT_DON", "COD_CAD_DON", "DIAG_KI", 
#   "ACUTE_REJ_EPI_KI", "CREAT_TRR", "ACADEMIC_PRG_TCR", 
#   "ACADEMIC_LEVEL_TCR", "HLAMIS", "SERUM_CREAT", "DIAB",
#   "EDUCATION", "ABO_MAT", "PREV_TX_ANY", "URINE_INF_DON", 
#   "CDC_RISK_HIV_DON", "DRUGTRT_COPD", "REGION", "BUN_DON", 
#   "TOT_SERUM_ALBUM","WORK_INCOME_TCR", "MED_COND_TRR", "LOS",
#   "TRTREJ1Y_KI", "TRTREJ6M_KI", "TX_PROCEDUR_TY_KI",
#   "INIT_STAT", "END_STAT", "NPKID")
# 
# library(mediation)
# library(lavaan)
# 
# multimed("WaitlistDuration", c("BMI", "ABO", "ABO_MAT"), treat = "Group", data = balanced_data[1:100,])
# 
# # Specify the mediation model
# waitlist_mediation <- "# a path
#          BMI ~ a * Group
# 
#          # b path
#          WaitlistDuration ~ b * BMI
# 
#          # c prime path 
#          WaitlistDuration ~ cp * Group
# 
#          # indirect and total effects
#          ab := a * b
#          total := cp + ab"
# 
# # Fit the mediation model
# mediation_fit <- sem(waitlist_mediation, data = balanced_data)
# # 
# # # Examine the model results
# # summary(mediation_fit)
# 
# 
# library(mma)
# 
# #which(colnames(balanced_data)=="WaitlistDuration")
# 
# balanced_data$ABO = droplevels(balanced_data$ABO)
# balanced_data$ON_DIALYSIS = droplevels(balanced_data$ON_DIALYSIS)
# 
# pred_ETH <- balanced_data$Group
# 
# x <- data.frame(balanced_data[, c(5:10)])
# y <- data.frame(balanced_data[,2])
# 
# data.org(x, y, pred_ETH, alpha = .05)
# 
# mma_waitlist <- mma(x, y, pred_ETH, mediator = 2:6, para = F, n2=2)
# 
# 
# 
# 
# library(boot)
# library(marginaleffects)
# 
# 
# 
# balanced_data <- balanced_data %>% mutate(Asian = case_when(
#   Group == "Asian" ~ 1,
#   Group != "Asian" ~ 0
# ), Black = case_when(
#   Group == "Black" ~ 1,
#   Group != "Black" ~ 0
# ),
#   Hispanic = case_when(
#     Group == "Hispanic" ~ 1,
#     Group != "Hispanic" ~ 0
#   ),
#   White = case_when(
#     Group == "White" ~ 1,
#     Group != "White" ~ 0
#   )
# )
# 
# #cbind(balanced_data$Asian, balanced_data$Black, balanced_data$Hispanic, balanced_data$White)
# 
# 
# 
# 
# WL_Asian <- WaitlistDuration ~ Asian + BMI + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI
# WL_Black <- WaitlistDuration ~ Black + BMI + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI
# WL_Hispanic <- WaitlistDuration ~ Hispanic + BMI + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI
# WL_White <- WaitlistDuration ~ White + BMI + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI
# 
# Asian_fit <- glm(WL_Asian, data = balanced_data)
# Black_fit <- glm(WL_Black, data = balanced_data)
# Hispanic_fit <- glm(WL_Hispanic, data = balanced_data)
# White_fit <- glm(WL_White, data = balanced_data)
# 
# comp_Asian <- comparisons(Asian_fit, variables = list(Asian = 0:1))
# comp_Black <- comparisons(Black_fit, variables = list(Black = 0:1))
# comp_Hispanic <- comparisons(Hispanic_fit, variables = list(Hispanic = 0:1))
# comp_White <- comparisons(White_fit, variables = list(White = 0:1))
# 
# 
# 
# 
# Group <- balanced_data$Group
# WL_duration <- c()
# 
# for (i in 1:length(Group)){
#   if (Group[i] == "Asian"){
#     WL_duration[i] = comp_Asian$predicted[i]
#   } else if (Group[i] == "Black"){
#     WL_duration[i] = comp_Black$predicted[i]
#   } else if (Group[i] == "Hispanic"){
#     WL_duration[i] = comp_Hispanic$predicted[i]
#   } else if (Group[i] == "White"){
#     WL_duration[i] = comp_White$predicted[i]
#   }
# }
# 
# summary(WL_duration)
# 
# balanced_data$Adjusted_WL <- WL_duration


###############################################################################################
############################# Neyman Orthogonality Debiasing (not this) #######################
###############################################################################################
# 
# library(randomForest)
# library(caret)
# 
# ######################################## Waitlist #############################################
# set.seed(123) # For reproducibility
# # Create 5 equally sized folds
# folds <- createFolds(balanced_data$Group, k = 5, list = TRUE, returnTrain = TRUE)
# # Prepare a container for the nuisance estimates
# predicted_waitlist <- rep(NA, nrow(balanced_data))
# # Cross-fitting: Estimate nuisance parameters
# for(k in seq_along(folds)) {
#   # Split the data
#   train_indices <- folds[[k]]
#   test_indices <- setdiff(seq_len(nrow(balanced_data)), train_indices)
#   
#   train_data <- balanced_data[train_indices, ]
#   test_data <- balanced_data[test_indices, ]
#   
#   # Train the model on training data excluding the k-th fold
#   model_WL <- randomForest(WaitlistDuration ~ Group + BMI + PRA + ABO + ABO_MAT + ON_DIALYSIS +
#                              PRI_PAYMENT_TRR_KI + GENDER + GFR + DISTANCE + WGT_KG_CALC + HGT_CM_CALC,
#                            data = train_data)
#   
#   # Predict on the left-out fold
#   predicted_waitlist[test_indices] <- predict(model_WL, newdata = test_data)
# }
# # Compute residuals for the outcome
# residuals_waitlist <- balanced_data$WaitlistDuration - predicted_waitlist
# 
# mean_residual <- mean(residuals_waitlist, na.rm = TRUE)
# adjusted_waitlist_times <- predicted_waitlist + mean_residual
# 
# # Model for the treatment effect using the residuals
# # effect_model <- lm(residuals_waitlist ~ balanced_data$Group)
# # # Output the summary of the model to see the coefficients
# # summary(effect_model)
# 
# ######################################## KDPI #############################################
# set.seed(123) # For reproducibility
# 
# folds <- createFolds(balanced_data$Group, k = 5, list = TRUE, returnTrain = TRUE)
# 
# # Prepare a container for the nuisance estimates
# predicted_quality <- rep(NA, nrow(balanced_data))
# 
# # Cross-fitting: Estimate nuisance parameters
# for(k in seq_along(folds)) {
#   # Split the data
#   train_indices <- folds[[k]]
#   test_indices <- setdiff(seq_len(nrow(balanced_data)), train_indices)
#   
#   train_data <- balanced_data[train_indices, ]
#   test_data <- balanced_data[test_indices, ]
#   
#   # Train the model on training data excluding the k-th fold
#   model_kdpi <- randomForest(QualityScore ~ Group + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI +
#                                AGE_DON + AGE + CREAT_DON + HIST_DIABETES_DON + HIST_HYPERTENS_DON + GENDER,
#                              data = train_data)
#   
#   # Predict on the left-out fold
#   predicted_quality[test_indices] <- predict(model_kdpi, newdata = test_data)
# }
# 
# residuals_quality <- balanced_data$QualityScore - predicted_quality
# 
# mean_residual <- mean(residuals_quality, na.rm = TRUE)
# adjusted_quality_times <- predicted_quality + mean_residual
# 
# # # Model for the treatment effect using the residuals
# # effect_model_quality <- lm(residuals_quality ~ balanced_data$Group)
# # 
# # # Output the summary of the model to see the coefficients
# # summary(effect_model_quality)
# 
# ######################################## Graft Lifespan #############################################
# set.seed(123)
# 
# folds <- createFolds(balanced_data$Group, k = 5, list = TRUE, returnTrain = TRUE)
# 
# predicted_outcome <- rep(NA, nrow(balanced_data))
# # Cross-fitting: Estimate nuisance parameters
# for (k in seq_along(folds)) {
#   # Split the data
#   train_indices <- folds[[k]]
#   test_indices <- setdiff(seq_len(nrow(balanced_data)), train_indices)
#   train_data <- balanced_data[train_indices, ]
#   test_data <- balanced_data[test_indices, ]
#   # Train the random forest model on training data excluding the k-th fold
#   model_lifespan <- randomForest(OutcomeScore ~ Group + BMI + PRA + ABO + ON_DIALYSIS +
#                                    PRI_PAYMENT_TRR_KI + GENDER + COLD_ISCH_KI + AGE + GFR + 
#                                    SERUM_CREAT + PREV_KI_TX + LOS,
#                                  data = train_data)
#   # Predict on the left-out fold
#   predicted_outcome[test_indices] <- predict(model_lifespan, newdata = test_data)
# }
# 
# residuals_outcome <- balanced_data$OutcomeScore - predicted_outcome
# mean_residual <- mean(residuals_outcome, na.rm = TRUE)
# adjusted_outcome_scores <- predicted_outcome + mean_residual
# 
# ## Putting it together ##
# 
# treatment_effect <- data.frame(
#   Group = balanced_data$Group,
#   WaitlistDuration = adjusted_waitlist_times,
#   QualityScore = adjusted_quality_times,
#   OutcomeScore = adjusted_outcome_scores
# )
# 
# treatment_effect$WaitlistDuration <- scale_min_max(adjust_value(treatment_effect$WaitlistDuration))
# treatment_effect$QualityScore <- scale_min_max(adjust_value(treatment_effect$QualityScore))
# treatment_effect$OutcomeScore <- scale_min_max(adjust_value(treatment_effect$OutcomeScore))
# 
# 
# treatment_effect %>% group_by(Group) %>% 
#   summarise(median(WaitlistDuration), median(QualityScore), median(OutcomeScore))

###############################################################################################
###############################################################################################
###############################################################################################







###############################################################################################
##################################### Parametric G Formula ####################################
###############################################################################################
# library(randomForest)
# rownames(balanced_data) <- NULL
# 
# ######################################## Waitlist #############################################
# model_WL <- randomForest(WaitlistDuration ~ (Group) + BMI + PRA + ABO + ABO_MAT + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER + GFR + DISTANCE + WGT_KG_CALC + HGT_CM_CALC, 
#                          data = balanced_data)
# 
# 
# ethnic_groups <- unique(balanced_data$Group)
# # Create an empty vector to store adjusted waitlist durations
# balanced_data$Adjusted_WL <- NA
# # Predict waitlist duration for all observations
# predictions_WL <- predict(model_WL, newdata = balanced_data)
# 
# # Calculate the R-squared value
# (r_squared <- 1 - sum((balanced_data$WaitlistDuration - predictions_WL)^2) / sum((balanced_data$WaitlistDuration - mean(balanced_data$WaitlistDuration))^2))
# 
# for (eth in ethnic_groups) {
#   # Get the indices of observations belonging to the current ethnicity
#   eth_indices <- which(balanced_data$Group == eth)
#   # Assign the predicted waitlist durations to the Adjusted_WL column for the current ethnicity
#   balanced_data$Adjusted_WL[eth_indices] <- predictions_WL[eth_indices]
# }
# 
# balanced_data %>% 
#   group_by(Group) %>% 
#   summarise(median(Adjusted_WL), mean(Adjusted_WL), median(WaitlistDuration), mean(WaitlistDuration))
# 
# 
# 
# ######################################## KDPI #############################################
# model_kdpi <- randomForest(QualityScore ~ (Group) + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
#                              HIST_DIABETES_DON + HIST_HYPERTENS_DON + GENDER,
#                           data = balanced_data)
# ethnic_groups <- unique(balanced_data$Group)
# # Create an empty vector to store adjusted waitlist durations
# balanced_data$Adjusted_KDPI <- NA
# # Predict waitlist duration for all observations
# predictions_kdpi <- predict(model_kdpi, newdata = balanced_data)
# # Calculate the R-squared value
# (r_squared <- 1 - sum((balanced_data$QualityScore - predictions_kdpi)^2) / sum((balanced_data$QualityScore - mean(balanced_data$QualityScore))^2))
# 
# for (eth in ethnic_groups) {
#   # Get the indices of observations belonging to the current ethnicity
#   eth_indices <- which(balanced_data$Group == eth)
#   # Assign the predicted waitlist durations to the Adjusted_KDPI column for the current ethnicity
#   balanced_data$Adjusted_KDPI[eth_indices] <- predictions_kdpi[eth_indices]
# }
# 
# balanced_data %>% 
#   group_by(Group) %>% 
#   summarise(median(Adjusted_KDPI), mean(Adjusted_KDPI), median(QualityScore), mean(QualityScore))
# 
# 
# ######################################## Graft Lifespan #############################################
# model_lifespan <- randomForest(OutcomeScore ~ (Group) + BMI + PRA + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER + COLD_ISCH_KI + AGE + GENDER + GFR + SERUM_CREAT + PREV_KI_TX + LOS,
#                            data = balanced_data)
# ethnic_groups <- unique(balanced_data$Group)
# # Create an empty vector to store adjusted waitlist durations
# balanced_data$Adjusted_Lifespan <- NA
# # Predict waitlist duration for all observations
# predictions_Lifespan <- predict(model_lifespan, newdata = balanced_data)
# # Calculate the R-squared value
# (r_squared <- 1 - sum((balanced_data$OutcomeScore - predictions_Lifespan)^2) / sum((balanced_data$OutcomeScore - mean(balanced_data$OutcomeScore))^2))
# 
# for (eth in ethnic_groups) {
#   # Get the indices of observations belonging to the current ethnicity
#   eth_indices <- which(balanced_data$Group == eth)
#   # Assign the predicted waitlist durations to the Adjusted_Lifespan column for the current ethnicity
#   balanced_data$Adjusted_Lifespan[eth_indices] <- predictions_Lifespan[eth_indices]
# }
# 
# balanced_data %>% 
#   group_by(Group) %>% 
#   summarise(median(Adjusted_Lifespan), mean(Adjusted_Lifespan), median(OutcomeScore), mean(OutcomeScore))
# 
# 
# balanced_data %>% 
#   group_by(Group) %>% 
#   summarise(median(Adjusted_WL), median(Adjusted_KDPI), median(Adjusted_Lifespan))

###############################################################################################
############################### Neyman Orthogonalization (w/o cross-fits) #####################
###############################################################################################

library(randomForest)
library(caret)
library(dplyr)
library(nnet)

Groups <- balanced_data$Group

balanced_data$Group <- Groups
set.seed(123)

# Define the treatment variable and confounders
treatment_var <- "Group"

confounders_WL <- c("BMI", "PRA", "ABO", "ABO_MAT", "ON_DIALYSIS", "PRI_PAYMENT_TRR_KI",
                 "GENDER", "GFR", "DISTANCE", "WGT_KG_CALC", "HGT_CM_CALC")

confounders_KDPI <- c("BMI", "ABO", "ABO_MAT", "PRI_PAYMENT_TRR_KI", "AGE_DON", "AGE", "CREAT_DON", "HIST_DIABETES_DON",
                      "HIST_HYPERTENS_DON", "GENDER")

confounders_Lifespan <- c("BMI", "PRA", "ABO", "ON_DIALYSIS", "PRI_PAYMENT_TRR_KI", "GENDER", "COLD_ISCH_KI", "AGE", "GFR",
                          "SERUM_CREAT", "PREV_KI_TX", "LOS")

# Convert the treatment variable to numeric representation
balanced_data[[treatment_var]] <- as.numeric(balanced_data[[treatment_var]])

# Estimate the propensity scores using multinomial logistic regression
prop_model_WL <- multinom(paste0(treatment_var, " ~ ", paste(confounders_WL, collapse = " + ")),
                       data = balanced_data)

prop_model_KDPI <- multinom(paste0(treatment_var, " ~ ", paste(confounders_KDPI, collapse = " + ")),
                          data = balanced_data)

prop_model_Lifespan <- multinom(paste0(treatment_var, " ~ ", paste(confounders_Lifespan, collapse = " + ")),
                          data = balanced_data)

propensity_scores_WL <- prop_model_WL$fitted.values
propensity_scores_KDPI <- prop_model_KDPI$fitted.values
propensity_scores_Lifespan <- prop_model_Lifespan$fitted.values




# Estimate the outcome regression function using Random Forest
model_WL <- randomForest(WaitlistDuration ~ Group + BMI + PRA + ABO + ABO_MAT + ON_DIALYSIS +
                                PRI_PAYMENT_TRR_KI + GENDER + GFR + DISTANCE + WGT_KG_CALC + HGT_CM_CALC, data = balanced_data)

model_kdpi <- randomForest(QualityScore ~ (Group) + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
                             HIST_DIABETES_DON + HIST_HYPERTENS_DON + GENDER, data = balanced_data)

model_lifespan <- randomForest(OutcomeScore ~ Group + BMI + PRA + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER + COLD_ISCH_KI + AGE + GFR + 
                                  SERUM_CREAT + PREV_KI_TX + LOS, data = balanced_data)



# Calculate the residuals from the outcome model
residuals_WL <- balanced_data$WaitlistDuration - predict(model_WL)
residuals_KDPI <- balanced_data$QualityScore - predict(model_kdpi)
residuals_Lifespan <- balanced_data$OutcomeScore - predict(model_lifespan)


### THIS WORKED BUT HAD ISSUE WITH balanced_data[[treatment_var]] being 1,2,3,4 so subtracting prop scores was problematic
# # Construct the orthogonal score function
# score_function_WL <- (balanced_data[[treatment_var]] - propensity_scores_WL) * residuals_WL
# score_function_KDPI <- (balanced_data[[treatment_var]] - propensity_scores_KDPI) * residuals_KDPI
# score_function_Lifespan <- (balanced_data[[treatment_var]] - propensity_scores_Lifespan) * residuals_Lifespan
# 
# 
# Waitlist <- c()
# KDPI <- c()
# Lifespan <- c()
# 
# for (i in 1:nrow(score_function_WL)){
#   Waitlist[i] <- score_function_WL[i, balanced_data[[treatment_var]][i]]
#   KDPI[i] <- score_function_KDPI[i, balanced_data[[treatment_var]][i]]
#   Lifespan[i] <- score_function_Lifespan[i, balanced_data[[treatment_var]][i]]
# }
# 
# balanced_data$Adjusted_WL_Neyman <- scale_min_max(adjust_value(Waitlist))
# balanced_data$Adjusted_KDPI_Neyman <- scale_min_max(adjust_value(KDPI))
# balanced_data$Adjusted_Lifespan_Neyman <- scale_min_max(adjust_value(Lifespan))


# Convert the treatment variable to a series of binary indicators
group_indicators <- model.matrix(~ factor(balanced_data[[treatment_var]]) - 1, data = balanced_data)
colnames(group_indicators) <- paste("Group", c(1:4), sep="")

# Construct the orthogonal score function for each group, adjusting for binary indicators
score_function_WL <- group_indicators * (group_indicators - propensity_scores_WL) * residuals_WL
score_function_KDPI <- groupp_indicators * (group_indicators - propensity_scores_KDPI) * residuals_KDPI
score_function_Lifespan <- group_indicators * (group_indicators - propensity_scores_Lifespan) * residuals_Lifespan

# Aggregate scores (sum over columns since each row now contains a score for each group)
final_scores_WL <- rowSums(score_function_WL)
final_scores_KDPI <- rowSums(score_function_KDPI)
final_scores_Lifespan <- rowSums(score_function_Lifespan)

# Store the final adjusted scores back to the balanced_data dataframe
balanced_data$Adjusted_WL_Neyman <- scale_min_max(final_scores_WL)
balanced_data$Adjusted_KDPI_Neyman <- scale_min_max(final_scores_KDPI)
balanced_data$Adjusted_Lifespan_Neyman <- scale_min_max(final_scores_Lifespan)

####

# balanced_data$Adjusted_WL_Neyman <- Waitlist
# balanced_data$Adjusted_KDPI_Neyman <- KDPI
# balanced_data$Adjusted_Lifespan_Neyman <- Lifespan

balanced_data$Group <- Groups

balanced_data %>% group_by(Group) %>% 
  summarise(median(Adjusted_WL_Neyman), median(Adjusted_KDPI_Neyman), median(Adjusted_Lifespan_Neyman))

balanced_data %>% group_by(Group) %>%
  summarise(mean(Adjusted_WL_Neyman), mean(Adjusted_KDPI_Neyman), mean(Adjusted_Lifespan_Neyman))









# #Estimate the debiased treatment effect
# debiased_effect <- mean(score_function)
# 
# # Calculate the standard error of the debiased effect
# n <- nrow(balanced_data)
# se_debiased <- sqrt(var(score_function) / n)
# 
# # Print the debiased treatment effect and its standard error
# cat("Debiased Treatment Effect:", debiased_effect, "\n")
# cat("Standard Error:", se_debiased, "\n")

###############################################################################################
############################## Neyman Orthogonalization w/ 5-folds ############################
###############################################################################################

library(randomForest)
library(caret)
library(dplyr)
library(nnet)

Groups <- balanced_data$Group
balanced_data$Group <- Groups
set.seed(123)

# Define the treatment variable and confounders
treatment_var <- "Group"

confounders_WL <- c("BMI", "PRA", "ABO", "ABO_MAT", "ON_DIALYSIS", "PRI_PAYMENT_TRR_KI",
                    "GENDER", "GFR", "DISTANCE", "WGT_KG_CALC", "HGT_CM_CALC")

confounders_KDPI <- c("BMI", "ABO", "ABO_MAT", "PRI_PAYMENT_TRR_KI", "AGE_DON", "AGE", "CREAT_DON", "HIST_DIABETES_DON",
                      "HIST_HYPERTENS_DON", "GENDER")

confounders_Lifespan <- c("BMI", "PRA", "ABO", "ON_DIALYSIS", "PRI_PAYMENT_TRR_KI", "GENDER", "COLD_ISCH_KI", "AGE", "GFR",
                          "SERUM_CREAT", "PREV_KI_TX", "LOS")

# Convert the treatment variable to numeric representation
balanced_data[[treatment_var]] <- as.numeric(balanced_data[[treatment_var]])

library(caret)
library(dplyr)
library(randomForest)
library(nnet)

set.seed(123)
k <- 5  # Number of folds
folds <- createFolds(balanced_data$Group, k = k)

################################ WAITLIST ################################
# Placeholder lists to store models and scores
prop_models <- list()
outcome_models <- list()
scores <- list()

for(i in seq_along(folds)) {
  # Split the data
  training <- balanced_data[-folds[[i]],]
  testing <- balanced_data[folds[[i]],]
  # Fit propensity score model on training data
  prop_model <- multinom(paste0(treatment_var, " ~ ", paste(confounders_WL, collapse = " + ")), data = training)
  prop_models[[i]] <- prop_model
  # Predict on testing data
  propensity_scores <- predict(prop_model, newdata = testing, type = "probs")
  # Fit outcome model on training data
  model_WL <- randomForest(WaitlistDuration ~ ., data = training[, c("WaitlistDuration", treatment_var, confounders_WL)])
  outcome_models[[i]] <- model_WL
  # Predict outcomes and calculate residuals on testing data
  residuals_WL <- testing$WaitlistDuration - predict(model_WL, newdata = testing)
  # Compute scores
  group_indicators <- model.matrix(~ factor(testing[[treatment_var]]) - 1, data = testing)
  score_function_WL <- group_indicators * (group_indicators - propensity_scores) * residuals_WL
  scores[[i]] <- rowSums(score_function_WL)
}
final_scores_WL <- numeric(length = nrow(balanced_data))
# Aggregate scores from all folds
for(i in seq_along(folds)) {
  final_scores_WL[folds[[i]]] <- scores[[i]]
}
balanced_data$WL_Adjusted <- scale_min_max(final_scores_WL)

################################ KDPI ################################
prop_models_KDPI <- list()
kdpi_models <- list()
scores_kdpi <- list()
for(i in seq_along(folds)) {
  # Split the data
  training <- balanced_data[-folds[[i]],]
  testing <- balanced_data[folds[[i]],]
  # Fit propensity score model on training data
  prop_model_KDPI <- multinom(paste0(treatment_var, " ~ ", paste(confounders_KDPI, collapse = " + ")), data = training)
  prop_models_KDPI[[i]] <- prop_model_KDPI
  # Predict on testing data
  propensity_scores_KDPI <- predict(prop_model_KDPI, newdata = testing, type = "probs")
  
  # Fit outcome model on training data
  model_kdpi <- randomForest(QualityScore ~ ., data = training[, c("QualityScore", treatment_var, confounders_KDPI)])
  kdpi_models[[i]] <- model_kdpi
  
  # Predict outcomes and calculate residuals on testing data
  residuals_KDPI <- testing$QualityScore - predict(model_kdpi, newdata = testing)
  # Compute scores
  group_indicators <- model.matrix(~ factor(testing[[treatment_var]]) - 1, data = testing)
  score_function_KDPI <- group_indicators * (group_indicators - propensity_scores_KDPI) * residuals_KDPI
  scores_kdpi[[i]] <- rowSums(score_function_KDPI)
}
final_scores_KDPI <- numeric(length = nrow(balanced_data))
# Aggregate scores from all folds
for(i in seq_along(folds)) {
  final_scores_KDPI[folds[[i]]] <- scores_kdpi[[i]]
}
balanced_data$KDPI_Adjusted <- scale_min_max(final_scores_KDPI)

################################ Graft Lifespan ################################
prop_models_Lifespan <- list()
Lifespan_models <- list()
scores_Lifespan <- list()

for(i in seq_along(folds)) {
  # Split the data
  training <- balanced_data[-folds[[i]],]
  testing <- balanced_data[folds[[i]],]
  
  # Fit propensity score model on training data
  prop_model_Lifespan <- multinom(paste0(treatment_var, " ~ ", paste(confounders_Lifespan, collapse = " + ")), data = training)
  prop_models_Lifespan[[i]] <- prop_model_Lifespan
  
  # Predict on testing data
  propensity_scores_Lifespan <- predict(prop_model_Lifespan, newdata = testing, type = "probs")
  
  # Fit outcome model on training data
  model_Lifespan <- randomForest(OutcomeScore ~ ., data = training[, c("OutcomeScore", treatment_var, confounders_Lifespan)])
  Lifespan_models[[i]] <- model_Lifespan
  
  # Predict outcomes and calculate residuals on testing data
  residuals_Lifespan <- testing$OutcomeScore - predict(model_Lifespan, newdata = testing)
  
  # Compute scores
  group_indicators <- model.matrix(~ factor(testing[[treatment_var]]) - 1, data = testing)
  score_function_Lifespan <- group_indicators * (group_indicators - propensity_scores_Lifespan) * residuals_Lifespan
  scores_Lifespan[[i]] <- rowSums(score_function_Lifespan)
}

final_scores_Lifespan <- numeric(length = nrow(balanced_data))

# Aggregate scores from all folds
for(i in seq_along(folds)) {
  final_scores_Lifespan[folds[[i]]] <- scores_Lifespan[[i]]
}

balanced_data$Lifespan_Adjusted <- scale_min_max(final_scores_Lifespan)

# Summarize
balanced_data %>% group_by(Group) %>% summarise(median(WL_Adjusted), median(KDPI_Adjusted), median(Lifespan_Adjusted))
balanced_data %>% group_by(Group) %>% summarise(mean(WL_Adjusted), mean(KDPI_Adjusted), mean(Lifespan_Adjusted))



# # Fit the propensity score model (multinomial logistic regression)
# #prop_model <- multinom(Group ~ BMI + PRA + ABO + ABO_MAT + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER + GFR + DISTANCE + WGT_KG_CALC + HGT_CM_CALC, data = balanced_data)
# prop_model_KDPI <- multinom(paste0(treatment_var, " ~ ", paste(confounders_KDPI, collapse = " + ")),
#                             data = balanced_data)
# # Fit the outcome model (Random Forest)
# #outcome_model <- randomForest(Y ~ Group + BMI + PRA + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER + COLD_ISCH_KI + AGE + GFR + SERUM_CREAT + PREV_KI_TX + LOS, data = balanced_data, importance = TRUE, keep.forest = TRUE)
# model_kdpi <- randomForest(QualityScore ~ (Group) + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
#                              HIST_DIABETES_DON + HIST_HYPERTENS_DON + GENDER, data = balanced_data , importance = TRUE, keep.forest = TRUE)
# 
# # Calculate predicted propensity scores and predicted outcomes
# balanced_data$predicted_propensity <- predict(prop_model_KDPI, newdata = balanced_data, type = "probs")
# balanced_data$predicted_outcome <- predict(model_kdpi, newdata = balanced_data)
# 
# # Calculate the orthogonalized treatments and residuals
# balanced_data$V <- with(balanced_data, as.numeric(Group) - predicted_propensity[, as.numeric(Group)])
# balanced_data$R <- with(balanced_data, QualityScore - predicted_outcome)
# 
# # Adjust the outcomes
# balanced_data$adjusted_Y <- with(balanced_data, QualityScore + V * R)
# 
# balanced_data %>% group_by(Group) %>% summarise(median(adjusted_Y), median(Adjusted_KDPI_Neyman))
# Nuisance function for Waitlist duration
# model_WL <- randomForest(WaitlistDuration ~ Group + BMI + PRA + ABO + ABO_MAT + ON_DIALYSIS +
#                            PRI_PAYMENT_TRR_KI + GENDER + GFR + DISTANCE + WGT_KG_CALC + HGT_CM_CALC, data = balanced_data)
# balanced_data$predicted_WL <- predict(model_WL, newdata = balanced_data)
# 
# # Nuisance function for KDPI score
# model_kdpi <- randomForest(QualityScore ~ (Group) + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
#                              HIST_DIABETES_DON + HIST_HYPERTENS_DON + GENDER, data = balanced_data)
# balanced_data$predicted_KDPI <- predict(model_kdpi, newdata = balanced_data)
# 
# # Nuisance function for Graft Lifespan
# model_lifespan <- randomForest(OutcomeScore ~ Group + BMI + PRA + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER + COLD_ISCH_KI + AGE + GFR + 
#                                  SERUM_CREAT + PREV_KI_TX + LOS, data = balanced_data)balanced_data$predicted_KDPI <- predict(model_g_KDPI, newdata = balanced_data)
# balanced_data$predicted_Lifespan <- predict(model_lifespan, newdata = balanced_data)
# 
# 
# # Estimate m_0(X) for orthogonalized treatment
# model_m_WL <- randomForest(Group ~ BMI + PRA + ABO + ABO_MAT + ON_DIALYSIS +
#                              PRI_PAYMENT_TRR_KI + GENDER + GFR + DISTANCE + WGT_KG_CALC + HGT_CM_CALC, data = balanced_data)
# balanced_data$V_WL <- balanced_data$Group - predict(model_m_WL, newdata = balanced_data)
# 
# model_m_KDPI <- randomForest(Group ~ BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
#                                HIST_DIABETES_DON + HIST_HYPERTENS_DON + GENDER, data = balanced_data)
# balanced_data$V_KDPI <- balanced_data$Group - predict(model_m_KDPI, newdata = balanced_data)
# 
# model_m_Lifespan <- randomForest(Group ~ BMI + PRA + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER + COLD_ISCH_KI + AGE + GFR + 
#                                    SERUM_CREAT + PREV_KI_TX + LOS, data = balanced_data)
# balanced_data$V_Lifespan <- balanced_data$Group - predict(model_m_Lifespan, newdata = balanced_data)
# 
# # Calculate residuals
# balanced_data$resid_WL <- balanced_data$WaitlistDuration - balanced_data$predicted_WL
# balanced_data$resid_KDPI <- balanced_data$QualityScore - balanced_data$predicted_KDPI
# balanced_data$resid_Lifespan <- balanced_data$OutcomeScore - balanced_data$predicted_Lifespan
# 
# # Calculate adjusted scores
# balanced_data$adjusted_WL <- balanced_data$V_WL * balanced_data$resid_WL
# balanced_data$adjusted_KDPI <- balanced_data$V_KDPI * balanced_data$resid_KDPI
# balanced_data$adjusted_Lifespan <- balanced_data$V_Lifespan * balanced_data$resid_Lifespan
# 
# balanced_data$adjusted_WL <- scale_min_max(balanced_data$adjusted_WL)
# balanced_data$adjusted_KDPI <- scale_min_max(balanced_data$adjusted_KDPI)
# balanced_data$adjusted_Lifespan <- scale_min_max(balanced_data$adjusted_Lifespan)
# 
# balanced_data %>% group_by(Group) %>% summarise(median(adjusted_WL), median(adjusted_KDPI), median(adjusted_Lifespan))
# ################
# library(caret)
# library(nnet)  # For multinomial logistic regression
# set.seed(123)
# 
# # Split data into two parts
# index <- createDataPartition(balanced_data$Group, p = 0.5, list = FALSE)
# train_data <- balanced_data[index,]
# test_data <- balanced_data[-index,]
# 
# # Fit separate outcome models for each group
# models <- lapply(unique(train_data$Group), function(g) {
#   subset_data <- train_data[train_data$Group == g, ]
#   # randomForest(Y ~ ., data = subset_data, na.action = na.omit)
#   randomForest(QualityScore ~ (Group) + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
#                  HIST_DIABETES_DON + HIST_HYPERTENS_DON + GENDER, data = subset_data)
# })
# 
# # Predict the outcome for each individual under each treatment condition
# predictions <- lapply(models, function(model) {
#   predict(model, newdata = test_data, type = "response")
# })
# 
# # Combine predictions into a matrix
# predictions_matrix <- do.call(cbind, predictions)
# colnames(predictions_matrix) <- paste("Group", unique(train_data$Group), "Predicted_Y")
# 
# # Calculate individual treatment effects
# # Example: Calculate ITE for being in Group 1 vs each of the other groups
# ite_group_1_vs_others <- sweep(predictions_matrix, 2, predictions_matrix[, "Group 1"], "-")
# 
# # Store results
# test_data <- cbind(test_data, predictions_matrix, ite_group_1_vs_others)

# ######################################## Waitlist #############################################
# # model_all <- glm(WaitlistDuration ~ factor(Group) + PRA + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI, 
# #                  data = balanced_data, family = gaussian())
# 
# 
# model_all <- lm(WaitlistDuration ~ factor(Group) + BMI + PRA + ABO + ON_DIALYSIS + PRI_PAYMENT_TRR_KI + GENDER, 
#                  data = balanced_data)
# summary(model_all)
# 
# 
# # Create a new data frame for predictions
# new_data <- balanced_data
# 
# ethnic_groups <- unique(balanced_data$Group)
# 
# # Create an empty list to store predictions
# predictions <- list()
# 
# for (eth in ethnic_groups) {
#   # Set everyone to the current ethnicity in the loop
#   new_data$Group <- eth
#   
#   # Predict waitlist duration for this hypothetical scenario
#   predictions[[eth]] <- predict(model_all, newdata = new_data, type = "response")
# }
# Group <- balanced_data$Group
# WL_duration <- c()
# for (i in 1:length(Group)){
#   if (Group[i] == "Asian"){
#     WL_duration[i] = predictions$Asian[i]
#   } else if (Group[i] == "Black"){
#     WL_duration[i] = predictions$Black[i]
#   } else if (Group[i] == "Hispanic"){
#     WL_duration[i] = predictions$Hispanic[i]
#   } else if (Group[i] == "White"){
#     WL_duration[i] = predictions$White[i]
#   }
# }
# balanced_data$Adjusted_WL <- WL_duration
# balanced_data %>% group_by(Group) %>% summarise(median(Adjusted_WL), mean(Adjusted_WL), median(WaitlistDuration), mean(WaitlistDuration))
# 
# 
# ######################################## KDPI #############################################
# # model_all <- glm(QualityScore ~ factor(Group) + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
# #                    HIST_DIABETES_DON + HIST_HYPERTENS_DON, 
# #                  data = balanced_data, family = gaussian())
# 
# model_all <- lm(QualityScore ~ factor(Group) + BMI + ABO + ABO_MAT + PRI_PAYMENT_TRR_KI + AGE_DON + AGE + CREAT_DON +
#                    HIST_DIABETES_DON + HIST_HYPERTENS_DON, 
#                  data = balanced_data)
# # Create a new data frame for predictions
# new_data <- balanced_data
# ethnic_groups <- unique(balanced_data$Group)
# 
# # Create an empty list to store predictions
# predictions <- list()
# 
# for (eth in ethnic_groups) {
#   # Set everyone to the current ethnicity in the loop
#   new_data$Group <- eth
#   
#   # Predict waitlist duration for this hypothetical scenario
#   predictions[[eth]] <- predict(model_all, newdata = new_data, type = "response")
# }
# Group <- balanced_data$Group
# KDPI_adjusted <- c()
# for (i in 1:length(Group)){
#   if (Group[i] == "Asian"){
#     KDPI_adjusted[i] = predictions$Asian[i]
#   } else if (Group[i] == "Black"){
#     KDPI_adjusted[i] = predictions$Black[i]
#   } else if (Group[i] == "Hispanic"){
#     KDPI_adjusted[i] = predictions$Hispanic[i]
#   } else if (Group[i] == "White"){
#     KDPI_adjusted[i] = predictions$White[i]
#   }
# }
# balanced_data$Adjusted_KDPI <- KDPI_adjusted
# balanced_data %>% group_by(Group) %>% summarise(median(Adjusted_KDPI), mean(Adjusted_KDPI), median(QualityScore), mean(QualityScore))
# 
# 
# 
# ######################################## Graft Lifespan #############################################
# # model_all <- glm(OutcomeScore ~ factor(Group) + MED_COND_TRR + ABO + ABO_MAT, 
# #                  data = balanced_data, family = gaussian())
# 
# model_all <- lm(OutcomeScore ~ factor(Group) + MED_COND_TRR + ABO + ABO_MAT, 
#                  data = balanced_data)
# # Create a new data frame for predictions
# new_data <- balanced_data
# ethnic_groups <- unique(balanced_data$Group)
# 
# # Create an empty list to store predictions
# predictions <- list()
# 
# for (eth in ethnic_groups) {
#   # Set everyone to the current ethnicity in the loop
#   new_data$Group <- eth
#   
#   # Predict waitlist duration for this hypothetical scenario
#   predictions[[eth]] <- predict(model_all, newdata = new_data, type = "response")
# }
# 
# Group <- balanced_data$Group
# Lifespan_adjusted <- c()
# for (i in 1:length(Group)){
#   if (Group[i] == "Asian"){
#     Lifespan_adjusted[i] = predictions$Asian[i]
#   } else if (Group[i] == "Black"){
#     Lifespan_adjusted[i] = predictions$Black[i]
#   } else if (Group[i] == "Hispanic"){
#     Lifespan_adjusted[i] = predictions$Hispanic[i]
#   } else if (Group[i] == "White"){
#     Lifespan_adjusted[i] = predictions$White[i]
#   }
# }
# balanced_data$Adjusted_Lifespan <- Lifespan_adjusted
# balanced_data %>% group_by(Group) %>% summarise(median(Adjusted_Lifespan), mean(Adjusted_Lifespan), median(OutcomeScore), mean(OutcomeScore))
# 



#####
balanced_data$Waitlist_Score <- balanced_data$WL_Adjusted
balanced_data$Quality_Score <- balanced_data$KDPI_Adjusted
balanced_data$Outcome_Score <- balanced_data$Lifespan_Adjusted
  
#   
# balanced_data$Waitlist_Score <- balanced_data$Adjusted_WL
# balanced_data$Quality_Score <- balanced_data$Adjusted_KDPI
# balanced_data$Outcome_Score <- balanced_data$Adjusted_Lifespan

# ###############################################################################################
# ########################################## Double ML ##########################################
# ###############################################################################################
# 
# library(DoubleML)
# library(mlr3)
# 
# data_encoded <- model.matrix(~ Group - 1, data = balanced_data)
# data_combined <- cbind(balanced_data[, c("OutcomeScore", "WaitlistDuration", "QualityScore", confounders)], data_encoded)
# 
# 
# # balanced_data <- balanced_data[complete.cases(balanced_data), ]
# 
# 
# # Create a DoubleMLData object
# dml_data_Quality <- DoubleMLData$new(data_combined, y_col = "QualityScore", d_cols = c("GroupWhite","GroupAsian", "GroupBlack", "GroupHispanic"))
# dml_data_Waitlist <- DoubleMLData$new(data_combined, y_col = "WaitlistDuration", d_cols = c("GroupWhite","GroupAsian", "GroupBlack", "GroupHispanic"))
# dml_data_Outcome <- DoubleMLData$new(data_combined, y_col = "OutcomeScore", d_cols = c("GroupWhite","GroupAsian", "GroupBlack", "GroupHispanic"))
# 
# # Specify the learner for the outcome and input models
# learner <- lrn("regr.ranger", num.trees = 500, min.node.size = 2, max.depth = 5)
# 
# 
# # Estimate the causal effect using Double ML
# dml_plr_Quality <- DoubleMLPLR$new(dml_data_Quality, ml_l = learner, ml_m = learner)
# dml_plr_fit_Quality <- dml_plr_Quality$fit(store_predictions = T)
# #dml_plr_Quality$bootstrap()
# #dml_plr_Quality$confint(joint = T)
# 
# 
# dml_plr_Waitlist <- DoubleMLPLR$new(dml_data_Waitlist, ml_l = learner, ml_m = learner)
# dml_plr_fit_Waitlist <- dml_plr_Waitlist$fit(store_predictions = T)
# #dml_plr_Waitlist$bootstrap()
# #dml_plr_Waitlist$confint(joint = T)
# 
# dml_plr_Outcome <- DoubleMLPLR$new(dml_data_Outcome, ml_l = learner, ml_m = learner)
# dml_plr_fit_Outcome <- dml_plr_Outcome$fit(store_predictions = T)
# #dml_plr_Outcome$bootstrap()
# #dml_plr_Outcome$confint(joint = T)
# # scale_min_max(adjust_value(dml_plr_Outcome$confint(joint = T)))
# 
# # Predictions
# predicted_Waitlist <- dml_plr_fit_Waitlist$predictions$ml_l
# residuals_Waitlist <- data_combined$WaitlistDuration - predicted_Waitlist
# 
# predicted_Quality <- dml_plr_fit_Quality$predictions$ml_l
# residuals_Quality <- data_combined$QualityScore - predicted_Quality
# 
# predicted_outcome <- dml_plr_fit_Outcome$predictions$ml_l
# residuals_outcome <- data_combined$OutcomeScore - predicted_outcome
# 
# balanced_data <- balanced_data %>% mutate(Group_num = case_when(
#   Group == "White" ~ 1,
#   Group == "Asian" ~ 2,
#   Group == "Black" ~ 3,
#   Group == "Hispanic" ~ 4
# ))
# 
# Waitlist_Score <- c()
# Quality_Score <- c()
# Outcome_Score <- c()
# 
# for (i in 1:length(balanced_data$Group_num)){
#   Waitlist_Score[i] <- residuals_Waitlist[i, , balanced_data$Group_num[i]]
#   Quality_Score[i] <- residuals_Quality[i, , balanced_data$Group_num[i]]
#   Outcome_Score[i] <-  residuals_outcome[i, , balanced_data$Group_num[i]]
# }
# 
# 
# # Calculating R-squared
# (r2_Waitlist <- 1 - sum(Waitlist_Score^2) / sum((data_combined$WaitlistDuration - mean(data_combined$WaitlistDuration))^2))
# (r2_Quality <- 1 - sum(Quality_Score^2) / sum((data_combined$QualityScore - mean(data_combined$QualityScore))^2))
# (r2_Outcome <- 1 - sum(Outcome_Score^2) / sum((data_combined$OutcomeScore - mean(data_combined$OutcomeScore))^2))
# 
# 


###############################################################################################
###############################################################################################
###############################################################################################

#### Translation and Rescaling ####
# X_translated =X+∣Min(X)∣+ϵ
adjust_value <- function(x) {
  min_val <- min(x)
  x_translated <- x + abs(min_val) + 1 # adding 1 to ensure all values are positive
  return(x_translated)
}



treatment_effect <- data.frame(
  Group = balanced_data$Group,
  WaitlistDuration = balanced_data$Waitlist_Score,
  QualityScore = balanced_data$Quality_Score,
  OutcomeScore = balanced_data$Outcome_Score
)

# treatment_effect$WaitlistDuration <- scale_min_max(adjust_value(treatment_effect$WaitlistDuration))
# treatment_effect$QualityScore <- scale_min_max(adjust_value(treatment_effect$QualityScore))
# treatment_effect$OutcomeScore <- scale_min_max(adjust_value(treatment_effect$OutcomeScore))


### OR, alternatively ###
# Logarithmic Scaling (with Offset)
# log_transform_with_offset <- function(x) {
#   min_val <- min(x)
#   x_offset <- x + abs(min_val) + 1 # adding 1 to ensure all values are positive
#   x_log_transformed <- log(x_offset) # applying log transformation
#   return(x_log_transformed)
# }

# treatment_effect$WaitlistDuration <- log_transform_with_offset(treatment_effect$WaitlistDuration)
# treatment_effect$QualityScore <- log_transform_with_offset(treatment_effect$QualityScore)
# treatment_effect$OutcomeScore <- log_transform_with_offset(treatment_effect$OutcomeScore)

##########################################################################################
################################### DEA w/ Uncertainty ###################################
##########################################################################################
# set.seed(123)
# # Assuming uncertainty ranges for demonstration
# # Inputs
# waitlist_duration_min <- treatment_effect$WaitlistDuration * 0.9
# waitlist_duration_max <- treatment_effect$WaitlistDuration * 1.1
# quality_score_min <- treatment_effect$QualityScore * 0.9
# quality_score_max <- treatment_effect$QualityScore * 1.1
# 
# # Outputs
# outcome_score_min <- treatment_effect$OutcomeScore * 0.9
# outcome_score_max <- treatment_effect$OutcomeScore * 1.1
# 
# # Simulating uncertain variables using uniform distribution
# inputs_dml_uncertain <- cbind(
#   runif(n = nrow(treatment_effect), min = waitlist_duration_min, max = waitlist_duration_max),
#   runif(n = nrow(treatment_effect), min = quality_score_min, max = quality_score_max)
# )
# outputs_dml_uncertain <- matrix(
#   runif(n = nrow(treatment_effect), min = outcome_score_min, max = outcome_score_max),
#   ncol = 1
# )
# 
# library(Benchmarking)
# 
# dea_results_dml_uncertain <- Benchmarking::dea(X = inputs_dml_uncertain, Y = outputs_dml_uncertain, RTS = "VRS", ORIENTATION = "out")
# 
# balanced_data$Efficiency_uncertain <- dea_results_dml_uncertain$eff
# 
# 
# (DEA_summary <- balanced_data %>%
#     group_by(Group) %>%
#     filter(Efficiency_uncertain != -Inf) %>% 
#     summarise(
#       Median_Efficiency = median(Efficiency_uncertain),
#       Mean_Efficiency = mean(Efficiency_uncertain),
#       Variance_Efficiency = var(Efficiency_uncertain)
# ))





# Example: Simulating uncertainty based on normalized data

set.seed(123)  # Ensuring reproducibility

# Simulating uncertain variables using the normal distribution centered around the normalized scores
inputs_dml_uncertain <- cbind(
  WaitlistScore = rnorm(n = length(balanced_data$Waitlist_Score), mean = balanced_data$Waitlist_Score, sd = sqrt(var(balanced_data$Waitlist_Score))),
  QualityScore = rnorm(n = length(balanced_data$Quality_Score), mean = balanced_data$Quality_Score, sd = sqrt(var(balanced_data$Quality_Score)))
)
# Assuming outputs are also normalized and stored in balanced_data
outputs_dml_uncertain <- matrix(
  rnorm(n = length(balanced_data$Outcome_Score), mean = balanced_data$Outcome_Score, sd = sqrt(var(balanced_data$Outcome_Score))),
  ncol = 1
)

# Running DEA with uncertain data
library(Benchmarking)
dea_results_dml_uncertain <- Benchmarking::dea(X = inputs_dml_uncertain, Y = outputs_dml_uncertain, RTS = "VRS", ORIENTATION = "out")

balanced_data$Efficiency_uncertain <- dea_results_dml_uncertain$eff
# 
# 
(DEA_summary <- balanced_data %>%
    group_by(Group) %>%
    filter(Efficiency_uncertain != -Inf) %>%
    summarise(
      Median_Efficiency = median(Efficiency_uncertain),
      Mean_Efficiency = mean(Efficiency_uncertain),
      Variance_Efficiency = var(Efficiency_uncertain)
))


###

inputs_simulated <- cbind(
  WaitlistScore = rnorm(n = length(balanced_data$Waitlist_Score), mean = balanced_data$Waitlist_Score, sd = sqrt(var(balanced_data$Waitlist_Score))),
  QualityScore = rnorm(n = length(balanced_data$Quality_Score), mean = balanced_data$Quality_Score, sd = sqrt(var(balanced_data$Quality_Score)))
)
# Assuming outputs are also normalized and stored in balanced_data
outputs_simulated <- matrix(
  rnorm(n = length(balanced_data$Outcome_Score), mean = balanced_data$Outcome_Score, sd = sqrt(var(balanced_data$Outcome_Score))),
  ncol = 1
)


# Number of simulations
n_simulations <- 3
# Prepare to collect efficiency scores
efficiency_scores <- numeric(n_simulations)

# Run the Monte Carlo simulation
set.seed(123)
for (i in 1:n_simulations) {
  # Simulate inputs
  inputs_simulated <- cbind(
    WaitlistScore = rnorm(n = length(balanced_data$Waitlist_Score), mean = balanced_data$Waitlist_Score, sd = sqrt(var(balanced_data$Waitlist_Score))),
    QualityScore = rnorm(n = length(balanced_data$Quality_Score), mean = balanced_data$Quality_Score, sd = sqrt(var(balanced_data$Quality_Score)))
  )
  
  # Simulate outputs
  outputs_simulated <- matrix(
    rnorm(n = length(balanced_data$Outcome_Score), mean = balanced_data$Outcome_Score, sd = sqrt(var(balanced_data$Outcome_Score))),
    ncol = 1
  )
  
  # Run DEA
  dea_results <- dea(X = inputs_simulated, Y = outputs_simulated, RTS = "VRS", ORIENTATION = "out")
  
  effs <- dea_results$eff
  # Replace -Inf with NA
  effs[effs == -Inf] <- NA
  
  # Store the average efficiency score for this simulation
  efficiency_scores[i] <- mean(effs, na.rm = TRUE)
}

# Analyze the results
hist(efficiency_scores, main = "Distribution of Efficiency Scores from Monte Carlo Simulation", xlab = "Efficiency Scores", breaks = 30)


#####
library(Benchmarking)

# Number of bootstrap samples
n_bootstraps <- 2

# Prepare to collect efficiency scores
efficiency_scores <- list()

# Run the bootstrap
set.seed(123)  # For reproducibility
for (i in 1:n_bootstraps) {
  # Create a bootstrap sample of the indices
  boot_indices <- sample(1:nrow(balanced_data), replace = TRUE, size = nrow(balanced_data))
  
  # Generate the bootstrapped data sets
  boot_inputs <- as.matrix(balanced_data[boot_indices, c("Waitlist_Score", "Quality_Score")])
  boot_outputs <- matrix(balanced_data[boot_indices, "Outcome_Score"], ncol = 1)
  
  # Run DEA
  dea_results <- dea(X = boot_inputs, Y = boot_outputs, RTS = "VRS", ORIENTATION = "in")
  
  # Store the efficiency scores for this bootstrap
  efficiency_scores[[i]] <- dea_results$eff
}

# Analyze the results
all_efficiencies <- unlist(efficiency_scores)  
# Replace -Inf with NA
all_efficiencies[all_efficiencies == -Inf] <- NA

hist(all_efficiencies, main = "Distribution of Efficiency Scores from Bootstrapping", xlab = "Efficiency Scores", breaks = 30)
summary(all_efficiencies)


######
# Define uncertainty distributions for input and output variables
input1_dist <- runif(nrow(treatment_effect), min = 0.9, max = 1.1) # Example: Uniform distribution
input2_dist <- rnorm(nrow(treatment_effect), mean = 1, sd = 0.1) # Example: Normal distribution
output_dist <- runif(nrow(treatment_effect), min = 0.95, max = 1.05) # Example: Uniform distribution

# Modify input and output matrices with uncertainty distributions
inputs_dml_uncertain <- inputs_dml * cbind(input1_dist, input2_dist)
outputs_dml_uncertain <- outputs_dml * output_dist

# Implement the uncertain DEA model
dea_results_dml_uncertain <- Benchmarking::dea(X = inputs_dml_uncertain, Y = outputs_dml_uncertain, RTS = "VRS", ORIENTATION = "out")
##########################################################################################
##########################################################################################
##########################################################################################

inputs_dml <- as.matrix(treatment_effect[, c("WaitlistDuration", "QualityScore")])
outputs_dml <- as.matrix(treatment_effect[, c("OutcomeScore")])


dea_results_dml <- Benchmarking::dea(X = inputs_dml, Y = outputs_dml, RTS = "VRS", ORIENTATION = "out")
# install.packages("job")
job::job({dea_results_boot_out <- Benchmarking::dea.boot(X = inputs_dml, Y = outputs_dml, NREP = 100, RTS = "VRS", alpha = 0.05, ORIENTATION = "out")})


balanced_data$Efficiency <- dea_results_dml$eff
balanced_data$Efficiency <- dea_results_boot$eff.bc

balanced_data$boots_median <- apply(dea_results_boot$boot, 1, median)
balanced_data$bias <- dea_results_boot$bias

balanced_data$lower_ci <- dea_results_boot$conf.int[,1]
balanced_data$upper_ci <- dea_results_boot$conf.int[,2]

(DEA_summary <- balanced_data %>%
    group_by(Group) %>%
    filter(Efficiency != -Inf) %>% 
    summarise(
      Median_Efficiency = median(Efficiency),
      Mean_Efficiency = mean(Efficiency),
      Variance_Efficiency = var(Efficiency)
))


(DEA_summary <- balanced_data %>%
  group_by(Group) %>%
  summarise(
    Median_Efficiency = median(Efficiency),
    Mean_Efficiency = mean(Efficiency),
    Variance_Efficiency = var(Efficiency),
    # I THINK IT SHOULD BE THIS ***
    Lower_CI= median(lower_ci),
    Upper_CI = median(upper_ci)
    #Pop_prop = n()/nrow(balanced_data)
    # Lower_CI = sort(boots_median)[(length(boots_median) + 1)*0.025], # (N+1) * (0.025)
    # Upper_CI = sort(boots_median)[(length(boots_median) + 1)*0.975], # (N+1) * (0.975)
))

balanced_data$Group_factor <- Groups

balanced_data$Group_factor <- factor(balanced_data$Group_factor, levels = c("White","Black", "Hispanic", "Asian"))


# Efficiency Score Density Plots for each Group
ggplot(balanced_data, aes(x = Efficiency, fill = Group_factor, color = Group_factor)) +
  geom_density(alpha = 0.5) +  
  labs(title = "DEA Efficiency Scores by Ethnic Group",
       x = "Efficiency",
       y = "Density") +
  theme_minimal() +  # Minimal theme
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

# Efficiency Score Box Plots for each Group
ggplot(balanced_data, aes(x = Efficiency, fill = Group_factor, color = Group_factor)) +
  geom_boxplot(alpha = 0.5) +  
  labs(title = "Boxplot of Efficiency by Group",
       x = "Efficiency",
       y = "Density") +
  theme_minimal() +  # Minimal theme
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")  





#effs_ci <- cbind(test_data$lower_ci, test_data$Efficiency, test_data$upper_ci)
## CI?
Benchmarking::dea.plot(inputs_dml, outputs_dml, txt=F)

# Bias-corrected ( MAIN ONE )
Benchmarking::dea.plot.frontier(dea_results_boot$eff.bc*inputs_dml, outputs_dml, col = "black", lty = "dashed", xlab = "Inputs (Priority & Access Fairness)", ylab = "Output (Outcome Fairness)")
Benchmarking::dea.plot.frontier((dea_results_boot$conf.int[,1])*inputs_dml, outputs_dml, add=TRUE, lty="dotted", col = "red", lwd=3)
Benchmarking::dea.plot.frontier((dea_results_boot$conf.int[,2])*inputs_dml, outputs_dml, add=TRUE, lty="dotted", col = "red", lwd=3)
##
peer_bests$Efficiency
# points(dea_results_boot$eff.bc*inputs_dml[which(dea_results_boot$eff.bc >= 0.9),], outputs_dml[which(dea_results_boot$eff.bc >= 0.9), ], col = "red", pch = 21, cex = 2, lwd = 3)


# Not Bias-corrected
Benchmarking::dea.plot.frontier(dea_results_boot$eff*inputs_dml, outputs_dml, col = "black", lty = "dashed")
Benchmarking::dea.plot.frontier((dea_results_boot$conf.int[,1]+dea_results_boot$bias)*inputs_dml, outputs_dml, add=TRUE, lty="dotted", col = "red", lwd=3)
Benchmarking::dea.plot.frontier((dea_results_boot$conf.int[,2]+dea_results_boot$bias)*inputs_dml, outputs_dml, add=TRUE, lty="dotted", col = "red", lwd=3)
##


library(dunn.test)
kruskal.test(Efficiency ~ Group, data = balanced_data)
dunn.test(balanced_data$Efficiency, balanced_data$Group, method = "bonferroni")


### PEERS Analysis ### 
group_type <- treatment_effect[, "Group"]
inputs_dml_peer <- as.matrix(treatment_effect[, c("WaitlistDuration", "QualityScore")])
outputs_dml_peer <- as.matrix(treatment_effect[, c("OutcomeScore")])

dea_results_dml_peer <- Benchmarking::dea(X = inputs_dml_peer, Y = outputs_dml_peer, RTS = "VRS", SLACK = TRUE)


#**** Not needed?
dea_results_dml_peer_crs <- Benchmarking::dea(X = inputs_dml_peer, Y = outputs_dml_peer, RTS = "CRS", SLACK = TRUE)

### Scale efficiency ### 
scale_eff <- dea_results_dml_peer_crs$eff / dea_results_dml_peer$eff


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Slacks
slacks <- slack(inputs_dml_peer, outputs_dml_peer, dea_results_dml_peer)
slack_df <- data.frame(Group = group_type, Efficiency = dea_results_dml_peer$eff, Slack = slacks$slack, Waitlist_slack = slacks$sx[,1], Quality_slack = slacks$sx[,2], Outcome_slack = slacks$sy[,1])

slacks_real <- slack(inputs_dml, outputs_dml, dea_results_dml)
slack_df <- data.frame(Group = group_type, Efficiency = dea_results_dml$eff, Slack = slacks_real$slack, Waitlist_slack = slacks_real$sx[,1], Quality_slack = slacks_real$sx[,2], Outcome_slack = slacks_real$sy[,1])


( slack_by_Group <- slack_df %>%
  filter(Waitlist_slack != 0 | Quality_slack != 0 | Outcome_slack != 0) %>% 
  group_by(Group) %>%
  summarise(median(Waitlist_slack), mean(Quality_slack), mean(Outcome_slack)) ) 



ggplot(slack_by_Group, aes(x = Group, y = `mean(Outcome_slack)`, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Group", y = "Average Slack", title = "Average Slack by Group")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



# Probably not useful
peer_analysis <- data.frame(peers(dea_results_dml))
peer_analysis$Group <- group_type
peer_analysis$Efficiency <- dea_results_dml$eff
peer_analysis$scale_eff <- scale_eff
(peer_analysis)

# Find best peers that anchor frontier
peer_nums <- get.number.peers(dea_results_dml)


Peer <- c()
Peer_count <- c()
for (i in 1:nrow(treatment_effect)){
  if (i %in% peer_nums$peer){
    Peer[i] <- 1
    Peer_count[i] <- peer_nums$count[which(peer_nums$peer == i)]
  } else {
    Peer[i] <- 0
    Peer_count[i] <- 0
  }
}

# treatment_effect$Peer <- Peer
treatment_effect$Peer_count <- Peer_count

treatment_effect$scale_eff <- scale_eff

treatment_effect %>% group_by(Group) %>% summarise(mean(Peer_count))




peer_bests <- treatment_effect[peer_nums$peer,] # Should include Group, WL, Quality, Outcome, Efficiency, Scale Eff.

peer_bests %>% group_by(Group) %>% summarise(n(), mean(Efficiency), mean(scale_eff))

(peer_analysis %>% group_by(Group) %>% summarise(median(Efficiency), median(scale_eff), n()))

# lambda(dea_results_dml_peer)

# data.frame(dea_results_dml_truncated_crs$eff, dea_results_dml_truncated$eff, scale_eff, slacks$slack, slacks$sx, slacks$sy)
dea.plot(x = inputs_dml_peer, y = outputs_dml_peer, RTS = "vrs", ORIENTATION = "in")

########################################################################
########################################################################
######################## DEA Frontier (Isoquant) #######################
########################################################################
########################################################################
# Inputs (e.g., Waitlist Duration)
x1 <- matrix(inputs_dml_peer[, "WaitlistDuration"], ncol=1)
x2 <- matrix(inputs_dml_peer[, "QualityScore"], ncol=1)
# Outputs (e.g., Outcome Score)
y <- outputs_dml_peer

########################################################################
#### For only including efficient units in the plots ################### 
########################################################################
# Extracting efficient units
efficient_units <- which(eff(dea_results_dml_peer) >= 0.7)

#efficient_units <- peer_nums$peer
dea.plot.isoquant(x1[efficient_units], x2[efficient_units], RTS = "vrs")
########################################################################
########################################################################


######################################################################################################
######################################################################################################
#### GOOD ISOQUANT ####
unique_groups <- unique(treatment_effect$Group)
group_colors <- setNames(rainbow(5)[2:5], unique_groups)


dea.plot(x1[efficient_units], x2[efficient_units], ORIENTATION = "in", RTS = "vrs",
         col = "black", lwd = 4, cex = 1.2,
         main = "DEA Frontier", xlab = "Waitlist", ylab = "KDPI", xlim = c(0,0.5), ylim = c(0, 0.5))

# Add colored points for each group
points(x1[efficient_units], x2[efficient_units], pch = 19, col = group_colors[treatment_effect$Group], cex = 1.2)
#Add lines from the origin to each point (if needed)
for (i in efficient_units){
  lines(c(0, x1[i]), c(0, x2[i]), col = "black", lwd = 0.5, lty = "dotted")
}

points(peer_bests$WaitlistDuration, peer_bests$QualityScore, col = "red", pch = 21, cex = 2, lwd = 3)
legend("topright", legend = names(group_colors), col = group_colors, pch = 19, cex = 0.8, title = "Group")
######################################################################################################
######################################################################################################

library(RColorBrewer)

color_with_alpha <- mapply(function(color, alpha) {
  rgb(red = col2rgb(color)[1,]/255, green = col2rgb(color)[2,]/255, blue = col2rgb(color)[3,]/255, alpha = alpha)
}, color = group_colors[treatment_effect$Group], alpha = dea_results_boot$eff.bc, SIMPLIFY = TRUE)


# This plot will display the efficiency frontier showing the relationship between waitlist duration and 
# graft lifespan, highlighting the best-performing DMUs in balancing these factors. Efficient points will
# be those where no other DMU can achieve more graft lifespan with shorter or the same waitlist duration or
# less graft lifespan with a longer waitlist duration, under the assumption of VRS.
dea.plot(x1[efficient_units], y[efficient_units], ORIENTATION = "in-out", RTS = "vrs",
         col = "black", lwd = 1, cex = .5,
         main = "DEA Frontier", xlab = "Waitlist", ylab = "Graft Lifespan")
points(x1[efficient_units], y[efficient_units], pch = 20, cex = 1.2, lwd = 3, col = color_with_alpha)
legend("bottomright", legend = names(group_colors), col = group_colors, pch = 19, cex = 0.8, title = "Group")



# This plot shows how efficiently the KDPI can be managed to maximize graft 
# lifespan. The efficiency frontier in this plot will illustrate which DMUs are best at achieving higher 
# graft lifespan for given levels of KDPI, under VRS conditions.
dea.plot(x2[efficient_units], y[efficient_units], ORIENTATION = "in-out", RTS = "vrs",
         col = "black", lwd = 1, cex = .2,
         main = "DEA Frontier", xlab = "KDPI", ylab = "Graft Lifespan", xlim = c(0,1), ylim = c(0.5, 1.2))


points(x2[efficient_units], y[efficient_units], pch = 20, cex = 1.2, lwd = 3, col = color_with_alpha)
legend("bottomright", legend = names(group_colors), col = group_colors, pch = 19, cex = 0.8, title = "Group")




# This plot demonstrates how DMUs can manage both Waitlist and KDPI efficiently without considering the
# output explicitly. It visually represents the best practices for minimizing these inputs, showing which 
# combinations of Waitlist and KDPI are used most efficiently together by the best-performing DMUs under VRS.
dea.plot(x1[efficient_units], x2[efficient_units], ORIENTATION = "in", RTS = "vrs",
         col = "black", lwd = 4, cex = 1.2,
         main = "DEA Frontier", xlab = "Waitlist", ylab = "KDPI")






# Add colored points for each group
points(x1[efficient_units], x2[efficient_units], pch = 20, cex = 1.2, lwd = 3, col = color_with_alpha)

# Add lines from the origin to each point (if needed)
for (i in efficient_units){
  lines(c(0, x1[i]), c(0, x2[i]), col = "black", lwd = 0.5, lty = "dotted")
}

points(peer_bests$WaitlistDuration, peer_bests$QualityScore, col = "red", pch = 21, cex = 2, lwd = 3)
legend("topright", legend = names(group_colors), col = group_colors, pch = 19, cex = 0.8, title = "Group")

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
plot3d(x1, x2, y, col = group_colors[treatment_effect$Group], size = 5)

# ##
# library(rgl)
# 
# # Assuming you have x1, x2 as inputs, and y as output, with DEA efficiency scores calculated
# # Find the efficient frontier (assuming efficiency scores are normalized from 0 to 1)
# treatment_effect$Efficiency <- dea_results_dml_truncated$eff
# frontier_points <- treatment_effect[treatment_effect$Efficiency >= 0.95, ]  # Adjust threshold as needed
# 
# 
# 
# fit <- gam(OutcomeScore ~ s(WaitlistDuration, QualityScore, k = 5), data = frontier_points)
# frontier_points <- treatment_effect[treatment_effect$Efficiency >= 0.90, ] # Adjust threshold as needed
# fit <- lm(OutcomeScore ~ poly(WaitlistDuration, degree = 2) + poly(QualityScore, degree = 2) + 
#             I(WaitlistDuration^2) + I(QualityScore^2) + I(WaitlistDuration * QualityScore), 
#           data = frontier_points)
# 
# # Create a grid to plot the plane
# # Create a grid to plot the curved surface
# grid_size <- 30
# x1_range <- seq(min(treatment_effect$WaitlistDuration), max(treatment_effect$WaitlistDuration), length.out = grid_size)
# x2_range <- seq(min(treatment_effect$QualityScore), max(treatment_effect$QualityScore), length.out = grid_size)
# grid <- expand.grid(WaitlistDuration = x1_range, QualityScore = x2_range)
# 
# grid$y <- predict(fit, newdata = grid)
# 
# # Plot the curved frontier surface
# # Plot all points
# plot3d(treatment_effect$WaitlistDuration, treatment_effect$QualityScore, treatment_effect$OutcomeScore, col = "blue", size = 1, type = 's')
# 
# # Add frontier points in a different color
# points3d(frontier_points$WaitlistDuration, frontier_points$QualityScore, frontier_points$OutcomeScore, col = "red", size = 10)
# 
# surface3d(matrix(grid$WaitlistDuration, nrow = grid_size, ncol = grid_size),
#           matrix(grid$QualityScore, nrow = grid_size, ncol = grid_size),
#           matrix(grid$y, nrow = grid_size, ncol = grid_size),
#           color = "green", alpha = 0.5) # Semi-transparent surface



# or level/contour plot
library(akima)

# Interpolating the output based on inputs
interp <- with(treatment_effect, interp(x1, x2, y))

# Plotting
filled.contour(interp, color.palette = heat.colors, 
               xlab = "Input 1 (x1)", ylab = "Input 2 (x2)", 
               main = "Contour Plot of Output Based on Two Inputs")

#########################################################################################
################### Producing Efficiency Density Plots for each group ################### 
#########################################################################################

densityplot(~Efficiency, data=balanced_data,
            groups=Group,
            main="DEA Efficiency Densities by Ethnic Group",
            plot.points=FALSE,
            auto.key=TRUE)
#########################################################################################
#########################################################################################
#########################################################################################


# With Conf Int.
library(lattice)
library(latticeExtra)
# Generate density plots for each variable
plot_efficiency <- densityplot(~Efficiency, data = balanced_data,
                               groups=Group, plot.points=FALSE,
                               auto.key=TRUE, main="DEA Efficiency Densities by Ethnic Group",
                               xlab="Values", ylab="Density")

plot_lower_ci <- densityplot(~lower_ci, data=balanced_data,
                             groups=Group, plot.points=FALSE, main="", xlab="", ylab="", auto.key=TRUE, lty = "dashed", lwd = 3)

plot_upper_ci <- densityplot(~upper_ci, data=balanced_data,
                             groups=Group, plot.points=FALSE, main="", xlab="", ylab="", lty = "dashed", auto.key=TRUE, lwd = 3)

# Combine the plots into a single plot
combined_plot <- plot_efficiency +
  as.layer(plot_lower_ci) +
  as.layer(plot_upper_ci)

# Display the combined plot
print(combined_plot)

#########################################################################################



library(Benchmarking)
library(ggplot2)

# Assume 'dea_results' is a list containing DEA results for each group
dea_results <- list()

groups <- unique(test_data$Group)
colors <- rainbow(length(groups))  # Generate unique colors for each group

plot_list <- list()  # To store plots (if using ggplot)

for (i in seq_along(groups)) {
  # Filter data for the current group
  group_data <- balanced_data[balanced_data$Group == groups[i], ]
  
  
  inputs_dml <- as.matrix(group_data[, c("WaitlistDuration", "QualityScore")])
  outputs_dml <- as.matrix(group_data[, c("OutcomeScore")])
  
  # Perform DEA (assuming inputs and outputs are predefined columns)
  dea_results[[i]] <- Benchmarking::dea(X = inputs_dml, Y = outputs_dml, RTS = "VRS")
}

mean(dea_results[[1]]$eff) # White
mean(dea_results[[2]]$eff) # Black
mean(dea_results[[3]]$eff) # Hispanic
mean(dea_results[[4]]$eff) # Asian

### Using eff.dens plot ###
# plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 10), xlab = "Efficiency", ylab = "Density", 
#      main = "Density of DEA Efficiencies")

# Loop through each set of DEA results
for (i in 1:length(dea_results)) {
  # Calculate density
  dens <- eff.dens(dea_results[[i]])   # Adjust this depending on how your DEA results are stored
  
  # Use 'plot' for the first group, 'lines' for subsequent groups
  if (i == 1) {
    plot(dens, type = "l", lwd = 2, col = colors[i], xlab = "Efficiency", ylab = "Density", main = "DEA Efficiency Densities by Ethnic Group", ylim = c(0,4))
  } else {
    lines(dens, type = "l", lwd = 2, col = colors[i])
  }
}

# Add a legend to the plot to differentiate the groups
legend("topright", legend = groups, col = colors, lwd = 2, title = "Groups")

