library(tidyverse)
library(bnlearn)
library(mice)
library(survival)
library(mediation)
library(randomForest)
library(dplyr)

setwd("/Users/alikmofrad/UCLA PhD/Research/Kidney_Exchange")


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
  
  NPKID = kidney_data$NPKID
)



sum(is.na(kidney_data$LOS))

dim(dmu_data[complete.cases(dmu_data),])
dmu_data <- dmu_data[complete.cases(dmu_data),]


library(randomForest)

# Prepare the data
confounders <- c("BMI", "WGT_KG_CALC", "ABO", "ON_DIALYSIS", "AGE", "AGE_DON", 
                 "HGT_CM_DON_CALC", "WGT_KG_DON_CALC", "ETHCAT_DON", 
                 "HIST_HYPERTENS_DON", "HYPERTENS_DUR_DON", "PRA", "GENDER", 
                 "GFR", "DISTANCE", "PRI_PAYMENT_TRR_KI", "PREV_KI_TX", 
                 "HIST_DIABETES_DON", "HIST_CANCER_DON", "HIST_CIG_DON", 
                 "HIST_COCAINE_DON", "CREAT_DON", "COD_CAD_DON", "DIAG_KI", 
                 "ACUTE_REJ_EPI_KI", "CREAT_TRR", "ACADEMIC_PRG_TCR", 
                 "ACADEMIC_LEVEL_TCR", "HLAMIS", "SERUM_CREAT", "DIAB",
                 "EDUCATION", "ABO_MAT", "PREV_TX_ANY", "URINE_INF_DON", 
                 "CDC_RISK_HIV_DON", "DRUGTRT_COPD", "REGION", "BUN_DON", 
                 "TOT_SERUM_ALBUM","WORK_INCOME_TCR", "MED_COND_TRR", "LOS",
                 "TRTREJ1Y_KI", "TRTREJ6M_KI", "TX_PROCEDUR_TY_KI",
                 "INIT_STAT", "END_STAT", "NPKID")


# KDPI_confounders <- c("REJTRT_KI")
# dmu_data <- dmu_data[ , c("Group", "WaitlistDuration", "QualityScore", "OutcomeScore", KDPI_confounders)]



dmu_data <- dmu_data[complete.cases(dmu_data), c("Group", "WaitlistDuration", "QualityScore", "OutcomeScore", confounders)]



## Setting the data subsample ##
# set.seed(14)
# sample_index <- sample(1:nrow(dmu_data), replace = F, size = 10000)
# dmu_data <- dmu_data[sample_index,]


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
rm(selected_vars)



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
########################################## Double ML ##########################################
###############################################################################################

library(DoubleML)
library(mlr3)

data_encoded <- model.matrix(~ Group - 1, data = balanced_data)
data_combined <- cbind(balanced_data[, c("OutcomeScore", "WaitlistDuration", "QualityScore", confounders)], data_encoded)


# balanced_data <- balanced_data[complete.cases(balanced_data), ]


# Create a DoubleMLData object
dml_data_Quality <- DoubleMLData$new(data_combined, y_col = "QualityScore", d_cols = c("GroupWhite","GroupAsian", "GroupBlack", "GroupHispanic"))
dml_data_Waitlist <- DoubleMLData$new(data_combined, y_col = "WaitlistDuration", d_cols = c("GroupWhite","GroupAsian", "GroupBlack", "GroupHispanic"))
dml_data_Outcome <- DoubleMLData$new(data_combined, y_col = "OutcomeScore", d_cols = c("GroupWhite","GroupAsian", "GroupBlack", "GroupHispanic"))

# Specify the learner for the outcome and input models
learner <- lrn("regr.ranger", num.trees = 500, min.node.size = 2, max.depth = 5)


# Estimate the causal effect using Double ML
dml_plr_Quality <- DoubleMLPLR$new(dml_data_Quality, ml_l = learner, ml_m = learner)
dml_plr_fit_Quality <- dml_plr_Quality$fit(store_predictions = T)
#dml_plr_Quality$bootstrap()
#dml_plr_Quality$confint(joint = T)


dml_plr_Waitlist <- DoubleMLPLR$new(dml_data_Waitlist, ml_l = learner, ml_m = learner)
dml_plr_fit_Waitlist <- dml_plr_Waitlist$fit(store_predictions = T)
#dml_plr_Waitlist$bootstrap()
#dml_plr_Waitlist$confint(joint = T)

dml_plr_Outcome <- DoubleMLPLR$new(dml_data_Outcome, ml_l = learner, ml_m = learner)
dml_plr_fit_Outcome <- dml_plr_Outcome$fit(store_predictions = T)
#dml_plr_Outcome$bootstrap()
#dml_plr_Outcome$confint(joint = T)
# scale_min_max(adjust_value(dml_plr_Outcome$confint(joint = T)))

# Predictions
predicted_Waitlist <- dml_plr_fit_Waitlist$predictions$ml_l
residuals_Waitlist <- data_combined$WaitlistDuration - predicted_Waitlist

predicted_Quality <- dml_plr_fit_Quality$predictions$ml_l
residuals_Quality <- data_combined$QualityScore - predicted_Quality

predicted_outcome <- dml_plr_fit_Outcome$predictions$ml_l
residuals_outcome <- data_combined$OutcomeScore - predicted_outcome




balanced_data <- balanced_data %>% mutate(Group_num = case_when(
  Group == "White" ~ 1,
  Group == "Asian" ~ 2,
  Group == "Black" ~ 3,
  Group == "Hispanic" ~ 4
))

Waitlist_Score <- c()
Quality_Score <- c()
Outcome_Score <- c()

for (i in 1:length(balanced_data$Group_num)){
  Waitlist_Score[i] <- residuals_Waitlist[i, , balanced_data$Group_num[i]]
  Quality_Score[i] <- residuals_Quality[i, , balanced_data$Group_num[i]]
  Outcome_Score[i] <-  residuals_outcome[i, , balanced_data$Group_num[i]]
}


# Calculating R-squared
(r2_Waitlist <- 1 - sum(Waitlist_Score^2) / sum((data_combined$WaitlistDuration - mean(data_combined$WaitlistDuration))^2))
(r2_Quality <- 1 - sum(Quality_Score^2) / sum((data_combined$QualityScore - mean(data_combined$QualityScore))^2))
(r2_Outcome <- 1 - sum(Outcome_Score^2) / sum((data_combined$OutcomeScore - mean(data_combined$OutcomeScore))^2))







treatment_effect <- data.frame(
  Group = balanced_data$Group,
  WaitlistDuration = Waitlist_Score,
  QualityScore = Quality_Score,
  OutcomeScore = Outcome_Score
)

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

# Min-Max Scaling
scale_min_max <- function(x) {
  x_scaled <- (x - min(x)) / (max(x) - min(x))
  return(x_scaled)
}

treatment_effect$WaitlistDuration <- scale_min_max(adjust_value(treatment_effect$WaitlistDuration))
treatment_effect$QualityScore <- scale_min_max(adjust_value(treatment_effect$QualityScore))
treatment_effect$OutcomeScore <- scale_min_max(adjust_value(treatment_effect$OutcomeScore))


### OR, alternatively ###
# Logarithmic Scaling (with Offset)
# log_transform_with_offset <- function(x) {
#   min_val <- min(x)
#   x_offset <- x + abs(min_val) + 1 # adding 1 to ensure all values are positive
#   x_log_transformed <- log(x_offset) # applying log transformation
#   return(x_log_transformed)
# }
# 
# treatment_effect$WaitlistDuration <- log_transform_with_offset(treatment_effect$WaitlistDuration)
# treatment_effect$QualityScore <- log_transform_with_offset(treatment_effect$QualityScore)
# treatment_effect$OutcomeScore <- log_transform_with_offset(treatment_effect$OutcomeScore)


# set.seed(12)
# subsample <- sample(1:nrow(treatment_effect), replace = F, size = 10000)
# treatment_effect <- treatment_effect[subsample, ]
# tally(treatment_effect$Group, "proportion")
# Perform DEA analysis
#treatment_effect <- treatment_effect[1:1000,]

inputs_dml <- as.matrix(treatment_effect[, c("WaitlistDuration", "QualityScore")])
outputs_dml <- as.matrix(treatment_effect[, c("OutcomeScore")])


dea_results_dml <- Benchmarking::dea(X = inputs_dml, Y = outputs_dml, RTS = "VRS")

# install.packages("job")
job::job({dea_results_boot <- Benchmarking::dea.boot(X = inputs_dml, Y = outputs_dml, NREP = 1000, RTS = "VRS", alpha = 0.05)})




balanced_data$Efficiency <- dea_results_boot$eff.bc

balanced_data$boots_median <- apply(dea_results_boot$boot, 1, median)
balanced_data$bias <- dea_results_boot$bias

balanced_data$lower_ci <- dea_results_boot$conf.int[,1]
balanced_data$upper_ci <- dea_results_boot$conf.int[,2]


(DEA_summary <- balanced_data %>%
  group_by(Group) %>%
  summarise(
    #median(Bias_Corrected_Efficiency),
    Median_Efficiency = median(Efficiency),
    # Lower_CI = sort(boots_median)[(length(boots_median) + 1)*0.025], # (N+1) * (0.025)
    # Upper_CI = sort(boots_median)[(length(boots_median) + 1)*0.975], # (N+1) * (0.975)
    Variance_Efficiency = var(Efficiency),
    # I THINK IT SHOULD BE THIS ***
    Lower_CI= median(lower_ci),
    Upper_CI = median(upper_ci)
    #Pop_prop = n()/nrow(balanced_data)
))

balanced_data$Group <- factor(balanced_data$Group, levels = c("White","Black", "Hispanic", "Asian"))


# Efficiency Score Density Plots for each Group
ggplot(balanced_data, aes(x = Efficiency, fill = Group, color = Group)) +
  geom_density(alpha = 0.5) +  
  labs(title = "DEA Efficiency Scores by Ethnic Group",
       x = "Efficiency",
       y = "Density") +
  theme_minimal() +  # Minimal theme
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")


# Efficiency Score Box Plots for each Group
ggplot(balanced_data, aes(x = Efficiency, fill = Group, color = Group)) +
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

