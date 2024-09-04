library(tidyverse)
library(bnlearn)
library(mice)
library(dplyr)
library(lubridate)
library(Benchmarking)
library(ggplot2)
library(gridExtra)
library(grid)
library(purrr)
library(randomForest)
library(caret)
library(nnet)

################################################################################
################################# DATA PREP ####################################
################################################################################

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

kidney_data$ETH_White <- as.numeric(kidney_data$ETHCAT == "White")
kidney_data$ETH_Asian <- as.numeric(kidney_data$ETHCAT == "Asian")
kidney_data$ETH_Black <- as.numeric(kidney_data$ETHCAT == "Black")
kidney_data$ETH_Hispanic <- as.numeric(kidney_data$ETHCAT == "Hispanic")

kidney_data$ETHCAT <- as.factor(kidney_data$ETHCAT)
kidney_data$ETHCAT_DON <- as.factor(kidney_data$ETHCAT_DON)

kidney_data$KDPI <-(sapply(kidney_data$KDPI, function(x) {
  if (is.na(x) || x == ".") {
    NA
  } else {
    as.numeric(sub("%", "", x)) / 100
  }
}))

kidney_data <- kidney_data[kidney_data$ON_DIALYSIS != "1",]
kidney_data <- kidney_data %>%
  filter(!is.na(ON_DIALYSIS)) %>%
  mutate(ON_DIALYSIS = ifelse(ON_DIALYSIS == "Y", 1, 0))

kidney_data$GTIME_KI <- as.numeric(kidney_data$GTIME_KI)
kidney_data$ON_DIALYSIS <- as.factor(kidney_data$ON_DIALYSIS)
kidney_data$ABO <- as.factor(kidney_data$ABO)
kidney_data$INIT_CPRA <- as.numeric(kidney_data$INIT_CPRA)
kidney_data$AGE <- as.numeric(kidney_data$AGE)
kidney_data$GENDER <- as.factor(kidney_data$GENDER)
kidney_data$PRA <- as.numeric(kidney_data$END_CPRA)
kidney_data$WGT_KG_CALC <- as.numeric(kidney_data$WGT_KG_CALC)
kidney_data$BMI <- as.numeric(kidney_data$BMI_CALC)
kidney_data$DIABETES_DON <- as.factor(kidney_data$DIABETES_DON)




kidney_data <- kidney_data %>% filter(ETHCAT != "Other", ETHCAT_DON != "Other")

kidney_data$ETHCAT <- droplevels(kidney_data$ETHCAT)
kidney_data$ETHCAT_DON <- droplevels((kidney_data$ETHCAT_DON))

kidney_data <- kidney_data %>% mutate(PRI_PAYMENT_TRR_KI = case_when(
  PRI_PAYMENT_TRR_KI %in% 1:14 ~ as.factor(PRI_PAYMENT_TRR_KI),
  !(PRI_PAYMENT_TRR_KI %in% 1:14) ~ NA),
  
  ACUTE_REJ_EPI_KI = case_when(
    ACUTE_REJ_EPI_KI %in% 1:3 ~ as.factor(ACUTE_REJ_EPI_KI),
    !(ACUTE_REJ_EPI_KI %in% 1:3) ~ NA),
  # GRF_FAIL_CAUSE_TY_KI = case_when(
  #    GRF_FAIL_CAUSE_TY_KI %in% c(1:12, 999) ~ as.factor(GRF_FAIL_CAUSE_TY_KI),
  #   !(GRF_FAIL_CAUSE_TY_KI %in% c(1:12, 999)) ~ NA),
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



kidney_data$NON_HRT_DON <- as.factor(kidney_data$NON_HRT_DON)
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
kidney_data$TX_PROCEDUR_TY_KI = as.factor(kidney_data$TX_PROCEDUR_TY_KI)
kidney_data$INIT_STAT = as.factor(kidney_data$INIT_STAT)
kidney_data$END_STAT = as.factor(kidney_data$END_STAT)
kidney_data$COLD_ISCH_KI = as.numeric(kidney_data$COLD_ISCH_KI)
kidney_data$HGT_CM_CALC = as.numeric(kidney_data$HGT_CM_CALC)




kidney_data <- kidney_data %>% mutate(HEP_C_ANTI_DON = case_when(
  HEP_C_ANTI_DON == "P" ~ 1,
  HEP_C_ANTI_DON == "N" ~ 0
))
kidney_data <- kidney_data %>% filter(!is.na(HEP_C_ANTI_DON))
kidney_data$HEP_C_ANTI_DON <- as.factor(kidney_data$HEP_C_ANTI_DON)


kidney_data$TX_YEAR <- year(mdy(kidney_data$TX_DATE))

kidney_data <- within(kidney_data, {
  ABO_A <- as.numeric(ABO == "A")
  ABO_B <- as.numeric(ABO == "B")
  ABO_AB <- as.numeric(ABO == "AB")
  ABO_O <- as.numeric(ABO == "O") # O is the reference category
})


kidney_data <- kidney_data %>%
  mutate(GRF_FAIL_CAUSE_TY_KI = case_when(
    GRF_FAIL_CAUSE_TY_KI == 1 ~ "Hyperacute Rejection",
    GRF_FAIL_CAUSE_TY_KI == 2 ~ "Acute Rejection",
    GRF_FAIL_CAUSE_TY_KI == 3 ~ "Primary Failure",
    GRF_FAIL_CAUSE_TY_KI == 4 ~ "Graft Thrombosis",
    GRF_FAIL_CAUSE_TY_KI == 5 ~ "Infection",
    GRF_FAIL_CAUSE_TY_KI == 6 ~ "Surgical Complications",
    GRF_FAIL_CAUSE_TY_KI == 7 ~ "Urological Complications",
    GRF_FAIL_CAUSE_TY_KI == 8 ~ "Recurrent Disease",
    GRF_FAIL_CAUSE_TY_KI == 9 ~ "Primary Non-Function (Graft Never Functioned Post-Transplant)", # Note: This is repeated for code 12
    GRF_FAIL_CAUSE_TY_KI == 10 ~ "Chronic Rejection",
    GRF_FAIL_CAUSE_TY_KI == 11 ~ "BK (Polyoma) Virus",
    GRF_FAIL_CAUSE_TY_KI == 12 ~ "Primary Non-Function (Graft Never Functioned Post-Transplant)",
    GRF_FAIL_CAUSE_TY_KI == 999 ~ "Other",
    TRUE ~ "Unknown" 
  ))
# 

kidney_data$GRF_FAIL_CAUSE_TY_KI = as.factor(kidney_data$GRF_FAIL_CAUSE_TY_KI)

kidney_data <- kidney_data %>% 
  mutate(DWFG_KI = case_when(
    DWFG_KI == "Y" ~ 1,
    DWFG_KI == "N" ~ 0
  ))

kidney_data$DWFG_KI <- as.factor(kidney_data$DWFG_KI)
kidney_data$GSTATUS_KI <- as.numeric(kidney_data$GSTATUS_KI)

kidney_data <- kidney_data %>% filter(ABO %in% c("A", "B", "AB", "O"))
kidney_data$ABO <- droplevels(as.factor(kidney_data$ABO))

kidney_data$BMI_DON_CALC <- as.numeric(kidney_data$BMI_DON_CALC)



kidney_data$COD_KI <- as.factor(kidney_data$COD_KI)

kidney_data <- kidney_data %>% mutate(COD_KI = case_when(
  COD_KI %in% c("3200", "3201", "3202", "3203" , "3204", "3299") ~ "Graft Fail", 
  COD_KI %in% c("3300", "3301", "3302", "3303", "3304", "3305", "3306", "3307", "3308", "3399") ~ "Infection",
  COD_KI %in% c("3400", "3401", "3402", "3499") ~ "Cardiovascular",
  COD_KI %in% c("3500", "3599") ~ "Cerebrovascular",
  COD_KI %in% c("3600", "3601", "3699") ~ "Hemorrhage",
  COD_KI %in% c("3700","3701","3702","3799") ~ "Malignancy",
  COD_KI %in% c("3800","3899") ~ "Trauma",
  COD_KI %in% c("3900","3901","3902","3903", "3904", "3905", "3906", "3907","3908", "3909", "3910", "3911", "3912",
                "3913","3914") ~ "Misc.",
  COD_KI == "3915" ~ "Primary non-function",
  COD_KI %in% c("3916","3917") ~ "Viral Infection",
  COD_KI == "."  ~ "Unknown"
))

kidney_data$CTR_CODE <- as.factor(kidney_data$CTR_CODE)

dmu_data <- data.frame(
  Group = kidney_data$ETHCAT,
  ETHCAT = kidney_data$ETHCAT,
  WaitlistDuration = kidney_data$DAYSWAIT_CHRON_KI, # Input 1
  QualityScore = kidney_data$KDPI, # Input 2
  OutcomeScore = kidney_data$GTIME_KI,# Output 1
  
  WL_days = kidney_data$DAYSWAIT_CHRON_KI,
  KDPI = kidney_data$KDPI,
  GTIME_KI = kidney_data$GTIME_KI,
  Year = kidney_data$TX_YEAR,
  DIABETES_DON = kidney_data$DIABETES_DON,
  # Confounders:
  BMI = kidney_data$BMI,
  BMI_CALC = kidney_data$BMI,
  WGT_KG_CALC = kidney_data$WGT_KG_CALC, 
  ABO = kidney_data$ABO, 
  ON_DIALYSIS = kidney_data$ON_DIALYSIS,
  AGE = kidney_data$AGE, 
  AGE_DON = kidney_data$AGE_DON, 
  HGT_CM_DON_CALC =  kidney_data$HGT_CM_DON_CALC, 
  WGT_KG_DON_CALC = kidney_data$WGT_KG_DON_CALC, 
  ETHCAT_DON = kidney_data$ETHCAT_DON, 
  HIST_HYPERTENS_DON = kidney_data$HIST_HYPERTENS_DON,
  HYPERTENS_DUR_DON = kidney_data$HYPERTENS_DUR_DON, 
  PRA = kidney_data$PRA, 
  GENDER = kidney_data$GENDER, 
  GFR = kidney_data$GFR, 
  DISTANCE = kidney_data$DISTANCE, 
  PRI_PAYMENT_TRR_KI = kidney_data$PRI_PAYMENT_TRR_KI, 
  PREV_KI_TX = kidney_data$PREV_KI_TX, 
  
  HIST_DIABETES_DON = kidney_data$HIST_DIABETES_DON,
  CTR_CODE = kidney_data$CTR_CODE,
  
  CREAT_DON = kidney_data$CREAT_DON, 
  COD_CAD_DON = kidney_data$COD_CAD_DON, 
  
  DIAG_KI = kidney_data$DIAG_KI, 
  BMI_DON_CALC = kidney_data$BMI_DON_CALC,
  
  NON_HRT_DON = kidney_data$NON_HRT_DON,
  
  GRF_FAIL_CAUSE_TY_KI = kidney_data$GRF_FAIL_CAUSE_TY_KI,
  
  ACADEMIC_LEVEL_TCR = kidney_data$ACADEMIC_LEVEL_TCR,
  HLAMIS = kidney_data$HLAMIS,
  
  SERUM_CREAT = kidney_data$SERUM_CREAT, 
  
  DIAB = kidney_data$DIAB,
  
  EDUCATION = kidney_data$EDUCATION,
  ABO_MAT = kidney_data$ABO_MAT,
  
  COLD_ISCH_KI = kidney_data$COLD_ISCH_KI, 
  
  PREV_TX_ANY = kidney_data$PREV_TX_ANY, 
  
  URINE_INF_DON = kidney_data$URINE_INF_DON, 
  
  
  DRUGTRT_COPD = kidney_data$DRUGTRT_COPD, 
  
  REGION = kidney_data$REGION, 
  
  
  HEP_C_ANTI_DON = kidney_data$HEP_C_ANTI_DON,
  MED_COND_TRR = kidney_data$MED_COND_TRR, 
  
  
  LOS = kidney_data$LOS,
  
  TRTREJ1Y_KI = kidney_data$TRTREJ1Y_KI,
  
  TRTREJ6M_KI = kidney_data$TRTREJ6M_KI,
  
  TX_PROCEDUR_TY_KI = kidney_data$TX_PROCEDUR_TY_KI,
  
  INIT_STAT = kidney_data$INIT_STAT, 
  END_STAT = kidney_data$END_STAT, 
  
  NPKID = kidney_data$NPKID,
  
  HGT_CM_CALC = kidney_data$HGT_CM_CALC,
  
  ABO_A = kidney_data$ABO_A,
  ABO_B = kidney_data$ABO_B,
  ABO_AB = kidney_data$ABO_B,
  ABO_O = kidney_data$ABO_O, # O is the reference category
  
  ETH_Asian = kidney_data$ETH_Asian,
  ETH_Black = kidney_data$ETH_Black,
  ETH_Hispanic = kidney_data$ETH_Hispanic,
  ETH_White = kidney_data$ETH_White,
  DWFG_KI = kidney_data$DWFG_KI,
  GSTATUS_KI = kidney_data$GSTATUS_KI,
  COD_KI = kidney_data$COD_KI
)

dmu_data <- dmu_data[complete.cases(dmu_data),]
dmu_data <- dmu_data %>% filter(Year >= 2011)

####################################################################################################
##########################################  RESAMPLING  ############################################
####################################################################################################
# Desired proportions based on 2020 U.S. Census data (https://data.census.gov/table?g=010XX00US)
Asian_n <- 19886049
Black_n <- 41104200
Hispanic_n <- 62080044
White_n <- 204277273
Eth_total <- Asian_n + Black_n + Hispanic_n + White_n
Asian_prop <- Asian_n/Eth_total
Black_prop <- Black_n/Eth_total
Hispanic_prop <- Hispanic_n/Eth_total
White_prop <- White_n/Eth_total
desired_props <- c(Asian = Asian_prop , Black = Black_prop , Hispanic = Hispanic_prop , White = White_prop)

# Function to perform stratified sampling
stratified_sample <- function(data, group_var, desired_props, total_samples) {
  desired_counts <- floor(desired_props * total_samples)
  grouped_data <- group_split(data, !!sym(group_var))
  sampled_data <- map2_dfr(grouped_data, names(desired_counts), function(group_data, group_name) {
    slice_sample(group_data, n = desired_counts[group_name])
  })
  return(sampled_data)
}

# Determine total sample size
group_sizes <- table(dmu_data$Group)
max_possible_samples <- min(group_sizes / desired_props[names(group_sizes)])
total_samples <- floor(max_possible_samples)

# Perform stratified sampling
set.seed(14)
balanced_data <- stratified_sample(dmu_data, "Group", desired_props, total_samples)

###############################################################################################
############################## Neyman Orthogonalization w/ 5-folds ############################
###############################################################################################
Groups <- balanced_data$Group
balanced_data$Group <- Groups
set.seed(123)

# Define the treatment variable and confounders
treatment_var <- "Group"

confounders_WL <- c("ABO", "PRA", "AGE", "GENDER", "PRI_PAYMENT_TRR_KI", "REGION")

confounders_KDPI <- c("AGE_DON", "WGT_KG_CALC", "HGT_CM_CALC", "HIST_HYPERTENS_DON", "DIABETES_DON", "COD_CAD_DON", "SERUM_CREAT", "HEP_C_ANTI_DON") 

confounders_Lifespan <- c("PRA", "AGE", "GENDER", "COLD_ISCH_KI", "HLAMIS", "AGE_DON", "PREV_KI_TX") 

balanced_data[[treatment_var]] <- as.numeric(balanced_data[[treatment_var]])

library(caret)
library(dplyr)
library(randomForest)
library(nnet)

set.seed(123)
k <- 5  # Number of folds
folds <- createFolds(balanced_data$Group, k = k)

################################ WAITLIST ################################
prop_models <- list()
outcome_models <- list()
scores <- list()

for(i in seq_along(folds)) {
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

balanced_data$WL_Adjusted <- final_scores_WL

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
balanced_data$KDPI_Adjusted <- final_scores_KDPI

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

balanced_data$Lifespan_Adjusted <- final_scores_Lifespan


balanced_data$Waitlist_Score <- balanced_data$WL_Adjusted
balanced_data$Quality_Score <- balanced_data$KDPI_Adjusted
balanced_data$Outcome_Score <- balanced_data$Lifespan_Adjusted

write.csv(balanced_data, "balanced_data.csv", row.names = F)

################################################################################
############################## Exploratory Data Analysis #######################
################################################################################
# Data Summaries for Fairness Criterion
balanced_data %>% 
  group_by(Group) %>% 
  summarise(
    Mean = format(mean(WaitlistDuration, na.rm = TRUE), nsmall = 3),
    SD = format(sd(WaitlistDuration, na.rm = TRUE), nsmall = 3),
    Median = format(median(WaitlistDuration, na.rm = TRUE), nsmall = 3),
    Q1 = format(quantile(WaitlistDuration, 0.25, na.rm = TRUE), nsmall = 3),
    Q3 = format(quantile(WaitlistDuration, 0.75, na.rm = TRUE), nsmall = 3)
  )

balanced_data %>% 
  group_by(Group) %>% 
  summarise(
    Mean = format(mean(QualityScore, na.rm = TRUE), nsmall = 3),
    SD = format(sd(QualityScore, na.rm = TRUE), nsmall = 3),
    Median = format(median(QualityScore, na.rm = TRUE), nsmall = 3),
    Q1 = format(quantile(QualityScore, 0.25, na.rm = TRUE), nsmall = 3),
    Q3 = format(quantile(QualityScore, 0.75, na.rm = TRUE), nsmall = 3)
  )

balanced_data %>% 
  group_by(Group) %>% 
  summarise(
    Mean = format(mean(OutcomeScore, na.rm = TRUE), nsmall = 3),
    SD = format(sd(OutcomeScore, na.rm = TRUE), nsmall = 3),
    Median = format(median(OutcomeScore, na.rm = TRUE), nsmall = 3),
    Q1 = format(quantile(OutcomeScore, 0.25, na.rm = TRUE), nsmall = 3),
    Q3 = format(quantile(OutcomeScore, 0.75, na.rm = TRUE), nsmall = 3)
  )

################################################################################
################################################################################
################################################################################



################################################################################
######################### DATA IMBALANCE PLOT (Fig. 2) #########################
################################################################################
balanced_data <- balanced_data %>% mutate(Group = case_when(
  Group == 1 ~ "Asian",
  Group == 2 ~ "Black",
  Group == 3 ~ "Hispanic",
  Group == 4 ~ "White"
))

balanced_data$Group <- factor(
  balanced_data$Group, levels = c("Asian", "Black", "Hispanic", "White")
)

# Prepare the data
original_data <- dmu_data %>%
  count(Group) %>%
  mutate(Percentage = n / sum(n),
         Dataset = "Imbalanced Data")

resampled_data <- balanced_data %>%
  count(Group) %>%
  mutate(Percentage = n / sum(n),
         Dataset = "U.S. Population")

combined_data <- rbind(original_data, resampled_data)


color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4",
                   "Hispanic" = "#6baed6", "White" = "#b3cde3")

ggplot(combined_data, aes(x = Dataset, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage * 100)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 6, fontface = "bold") +  
  scale_fill_manual(values = color_palette) +
  labs(x = "", y = "Relative Frequency") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 0.65)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold"),  
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20, face = "bold", margin = ggplot2::margin(r = 10)),  
        legend.title = element_text(size = 18, face = "bold"), 
        legend.text = element_text(size = 16)) 

################################################################################
################################################################################
################################################################################


