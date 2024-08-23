library(tidyverse)
library(bnlearn)
library(mice)
library(survival)
library(mediation)
library(randomForest)
library(dplyr)
library(lubridate)
library(Benchmarking)
library(ggplot2)
library(gridExtra)
library(grid)
library(purrr)


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

dmu_data <- data.frame(
  Group = kidney_data$ETHCAT,
  ETHCAT = kidney_data$ETHCAT,
  WaitlistDuration = kidney_data$DAYSWAIT_CHRON_KI, # DEA Input 1
  QualityScore = kidney_data$KDPI, # DEA Input 2
  OutcomeScore = kidney_data$GTIME_KI, # DEA Output
  WL_days = kidney_data$DAYSWAIT_CHRON_KI,
  KDPI = kidney_data$KDPI,
  GTIME_KI = kidney_data$GTIME_KI,
  Year = kidney_data$TX_YEAR,
  DIABETES_DON = kidney_data$DIABETES_DON,
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

treatment_effect <- data.frame(
  Group = balanced_data$Group,
  Year = balanced_data$Year,
  WaitlistDuration = balanced_data$Waitlist_Score,
  QualityScore = balanced_data$Quality_Score,
  OutcomeScore = balanced_data$Outcome_Score
)


################################################################################
################################# DEA ANALYSIS #################################
################################################################################
# Compute DEA Efficiency Scores
define_io <- function(data) {
  list(
    inputs = as.matrix(data[, c("WaitlistDuration", "QualityScore")]),
    outputs = as.matrix(data[, "OutcomeScore"])
  )
}
compute_efficiency <- function(data) {
  io <- define_io(data)
  eff <- dea(X = io$inputs, Y = io$outputs, RTS = "vrs", ORIENTATION = "graph")$eff
  data$efficiency <- eff
  return(data)
}

treatment_effect$WaitlistDuration <- treatment_effect$WaitlistDuration + 
  abs(min(treatment_effect$WaitlistDuration))
treatment_effect$QualityScore <- treatment_effect$QualityScore + 
  abs(min(treatment_effect$QualityScore))
treatment_effect$OutcomeScore <- treatment_effect$OutcomeScore + 
  abs(min(treatment_effect$OutcomeScore))

treatment_effect <- compute_efficiency(treatment_effect)

################################################################################
############# Plotting Distributions of Efficiency Scores by Group #############
################################################################################
library(ggplot2)
library(ggridges)

color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4", 
                   "Hispanic" = "#6baed6", "White" = "#b3cde3")

treatment_effect_Group_Label <- treatment_effect %>% mutate(Group = case_when(
  Group == 1 ~ "Asian",
  Group == 2 ~ "Black",
  Group == 3 ~ "Hispanic",
  Group == 4 ~ "White",
  TRUE ~ as.character(Group) 
))

treatment_effect_Group_Label$Group <- factor(
  treatment_effect_Group_Label$Group, 
  levels = c("White", "Hispanic", "Black", "Asian")
)
ggplot(treatment_effect_Group_Label, aes(x = efficiency, y = Group, fill = Group)) +
  geom_density_ridges(
    aes(height = after_stat(density)), 
    stat = "density",
    scale = 1, 
    alpha = 0.7, 
    color = "black"
  ) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = "Efficiency",
    y = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 12)  
  )

ggsave("Ridgeline_efficiencies.png", width = 10, height = 6, dpi = 300)
################################################################################
################################################################################
################################################################################



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
############################### DATA IMBALANCE PLOT ############################
################################################################################
library(ggplot2)
library(dplyr)

treatment_effect <- treatment_effect %>% mutate(Group = case_when(
  Group == 1 ~ "Asian",
  Group == 2 ~ "Black",
  Group == 3 ~ "Hispanic",
  Group == 4 ~ "White"
))

treatment_effect$Group <- factor(
  treatment_effect$Group, levels = c("Asian", "Black", "Hispanic", "White")
  )

# Prepare the data
original_data <- dmu_data %>%
  count(Group) %>%
  mutate(Percentage = n / sum(n),
         Dataset = "Imbalanced Data")

resampled_data <- treatment_effect %>%
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
            vjust = -0.5, size = 5) +  
  scale_fill_manual(values = color_palette) +
  labs(x = "", y = "Relative Frequency") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),  
        axis.title = element_text(size = 16),  
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12)) 

ggsave("group_distribution_comparison.png", width = 10, height = 6, dpi = 300)

################################################################################
################################################################################
################################################################################




################################################################################
################################## MMD TESTING #################################
################################################################################
library(maotai)
library(kernlab)


# Uncomment if not running data prep section code
# treatment_effect <- read.csv("treatment_effect_data.csv")

set.seed(14)
groups <- c("1", "2", "3", "4")


n_comparisons <- choose(length(groups), 2)  # 6 comparisons for 4 groups

# Bonferroni-corrected significance level
alpha <- 0.05
alpha_corrected <- alpha / n_comparisons

results <- list()

# Loop through all combinations of group pairs
for(i in 1:(length(groups) - 1)) {
  for(j in (i + 1):length(groups)) {
    group_i <- groups[i]
    group_j <- groups[j]
    
    # Subset data for each group
    theta_i <- treatment_effect[treatment_effect$Group == group_i,]
    theta_j <- treatment_effect[treatment_effect$Group == group_j,]

    combined_data <- rbind(theta_i, theta_j)
    
    # Compute the Euclidean distance matrix
    dmat <- as.matrix(dist(combined_data))
    
    # Calculate the median distance to use as sigma
    sigma <- median(dmat[dmat > 0])
    
    # Compute the Gaussian kernel matrix
    kmat <- exp(-(dmat^2) / (2 * sigma^2))
    
    # Create labels
    lab <- c(theta_i$Group, theta_j$Group)
    
    # Perform the MMD test
    set.seed(123)
    mmd_result <- mmd2test(kmat, lab, mc.iter = 999)
    
    results[[paste0("Group_", group_i, "_vs_Group_", group_j)]] <- list(
      mmd_statistic = mmd_result$statistic,
      original_p_value = mmd_result$p.value,
      significant = mmd_result$p.value < alpha_corrected
    )
  }
}

results

################################################################################
################################################################################
################################################################################



################################################################################
############################ PRIORITY ANALYSIS #################################
################################################################################

# Uncomment if not running data prep section code
#balanced_data <- read.csv("balanced_data.csv") 

data <- dplyr::select(balanced_data, WL_days, PRA, ABO_A, ABO_B, ABO_AB,
                      ABO_O,ABO,ETH_Asian, ETH_Black, ETH_Hispanic, ETH_White,
                      AGE, GENDER, ON_DIALYSIS, ETHCAT, ABO, GTIME_KI, 
                      ETHCAT_DON, MED_COND_TRR, KDPI, GFR, AGE_DON, 
                      BMI_DON_CALC, CREAT_DON, COD_CAD_DON, ETHCAT_DON,
                      DIABETES_DON, HGT_CM_DON_CALC, WGT_KG_DON_CALC,
                      BMI_CALC, COLD_ISCH_KI, HLAMIS)

data$ABO_A <- as.factor(data$ABO_A)
data$ABO_B <- as.factor(data$ABO_B)
data$ABO_AB <- as.factor(data$ABO_AB)
data$ABO_O <- as.factor(data$ABO_O)
data$ETHCAT <- factor(data$ETHCAT, 
                      levels = c("White", "Asian", "Black", "Hispanic"))
data$ETH_Asian <- as.factor(data$ETH_Asian)
data$ETH_Black <- as.factor(data$ETH_Black)
data$ETH_Hispanic <- as.factor(data$ETH_Hispanic)
data$ETH_White <- as.factor(data$ETH_White)
data$HLAMIS <- as.factor(data$HLAMIS)
data$ABO <- factor(data$ABO, levels = c("O", "A", "B", "AB"))


summary(lm(WL_days ~ PRA + ABO + ON_DIALYSIS + AGE + GENDER, data = data))
summary(aov(PRA ~ ETH_Asian, data = data))
summary(aov(PRA ~ ETH_Black, data = data))
summary(aov(PRA ~ ETH_Hispanic, data = data))
summary(aov(AGE ~ ETH_Asian, data = data))
summary(aov(AGE ~ ETH_Black, data = data))
summary(aov(AGE ~ ETH_Hispanic, data = data))
chisq.test(data$ABO_A, data$ETH_Asian)
chisq.test(data$ABO_A, data$ETH_Black)
chisq.test(data$ABO_A, data$ETH_Hispanic)
chisq.test(data$ABO_B, data$ETH_Asian)
chisq.test(data$ABO_B, data$ETH_Black)
chisq.test(data$ABO_B, data$ETH_Hispanic)
chisq.test(data$ABO_AB, data$ETH_Asian)
chisq.test(data$ABO_AB, data$ETH_Black)
chisq.test(data$ABO_AB, data$ETH_Hispanic)


library(mma)
covs <- data.frame(data[, c(12,13)])
x <- data[, c(2,7, 12:14)]
y <- data.frame(data[,1])
x$ABO <- factor(x$ABO, levels = c("O", "A", "B", "AB"))
pred_ETH <- data$ETHCAT

x <- as.data.frame(x)

mma_ETH <- mma(x, y, pred_ETH, mediator = c(1,2,5), catmed = c(2), 
               catref = "O", para = F, n2=1000)

summary(mma_ETH)

################################################################################
# Asymptotics simulation
mma_sizes <- list()
sizes <- seq(100, 13000, 100)

for (i in 1:length(sizes)){
  ind <- sample(1:nrow(x),sizes[i], replace = F)
  x_subset <- x[ind,]
  x_subset <- as.data.frame(x_subset)
  y_subset <- as.data.frame(y[ind,])
  y_subset <- y[ind,]
  pred_ETH_subset <- pred_ETH[ind]
  mma_sizes[[i]] <- mma::mma(x_subset, y_subset, pred_ETH_subset, 
                             mediator = c(1,2,5), catmed = c(2), para = F, n2=2)
  print(i)
}

indirect_all_Asian <- c()
indirect_all_Black <- c()
indirect_all_Hispanic <- c()

direct_Asian <- c()
direct_Black <- c()
direct_Hispanic <- c()

total_Asian <- c()
total_Black <- c()
total_Hispanic <- c()

for (i in 1:length(sizes)){
  indirect_all_Asian[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["ie"]][1,1]
  indirect_all_Black[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["ie"]][2,1]
  indirect_all_Hispanic[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["ie"]][3,1]
  
  direct_Asian[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["de"]][1]
  direct_Black[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["de"]][2]
  direct_Hispanic[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["de"]][3]
  
  total_Asian[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["te"]][1]
  total_Black[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["te"]][2]
  total_Hispanic[i] <- mma_sizes[[i]][["a.binx"]][["estimation"]][["te"]][3]
}


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

data <- data.frame(
  Size = rep(sizes, each = 3),
  Ethnic_Group = rep(c("Asian", "Black", "Hispanic"), times = length(sizes)),
  Indirect = c(indirect_all_Asian, indirect_all_Black, indirect_all_Hispanic),
  Direct = c(direct_Asian, direct_Black, direct_Hispanic),
  Total = c(total_Asian, total_Black, total_Hispanic)
)

# Define the color palette and line types
color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4", 
                   "Hispanic" = "#6baed6", "White" = "#b3cde3")
line_types <- c("Asian" = "solid", "Black" = "dashed", "Hispanic" = "dotdash")

# Create a base plot with different line types for each ethnic group
base_plot <- ggplot(data, aes(x = Size, y = Indirect, color = Ethnic_Group, 
                              linetype = Ethnic_Group)) +
  scale_color_manual(values = color_palette, name = "Group") + 
  scale_linetype_manual(values = line_types, name = "Group") +  
  theme_minimal(base_size = 15) + 
  theme(legend.position = "top",
        legend.key.width = unit(2, "cm"),  
        legend.key.height = unit(0.6, "cm"),  
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16)) 

legend <- g_legend(base_plot + geom_line(size = 1.5))

# Plot for Indirect Effect
plot_indirect <- base_plot +
  geom_line(aes(y = Indirect), size = 1.5) +  
  labs(title = "Indirect Effect", y = "") +
  theme(axis.title.x = element_blank(), legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0)) 

# Plot for Direct Effect
plot_direct <- base_plot +
  geom_line(aes(y = Direct), size = 1.5) +  
  labs(title = "Direct Effect", y = "") +
  theme(axis.title.x = element_blank(), legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0))  

# Plot for Total Effect
plot_total <- base_plot +
  geom_line(aes(y = Total), size = 1.5) +  
  labs(title = "Total Effect", y = "", x = "Sample Size (n)") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0))  

# Arrange the plots with a common legend and a main title
grid.arrange(
  legend,
  arrangeGrob(plot_indirect, plot_direct, plot_total, ncol = 1),
  heights = c(1, 10)
)

png("Effects_Convergence.png", width = 10, height = 6, units = "in", res = 300)

grid.arrange(
  legend,
  arrangeGrob(plot_indirect, plot_direct, plot_total, ncol = 1),
  heights = c(1, 10) 
)

dev.off()

################################################################################
################################################################################
################################################################################



################################################################################
############################## ACCESS ANALYSIS #################################
################################################################################
# Uncomment if not running data prep section code
#balanced_data <- read.csv("balanced_data.csv") 

library(caret)
library(randomForest)

access_data <- dplyr::select(balanced_data, WL_days, PRA, AGE, GENDER, 
                             ON_DIALYSIS, ETHCAT, ABO, GTIME_KI, ETHCAT_DON, 
                             MED_COND_TRR, KDPI, GFR, AGE_DON, BMI_DON_CALC,
                             CREAT_DON, COD_CAD_DON, ETHCAT_DON, DIABETES_DON,
                             HGT_CM_DON_CALC, WGT_KG_DON_CALC, BMI_CALC, 
                             DISTANCE, NON_HRT_DON, HIST_HYPERTENS_DON,
                             HIST_DIABETES_DON, HEP_C_ANTI_DON)

access_data$HIST_DIABETES_DON <- as.factor(access_data$HIST_DIABETES_DON)
access_data <- access_data %>% mutate(HIST_DIABETES_DON = case_when(
  HIST_DIABETES_DON %in% c(2,3,4,5) ~ 1,
  TRUE ~ 0
))
access_data$HIST_DIABETES_DON <- as.factor(access_data$HIST_DIABETES_DON)
access_data$ETHCAT <- factor(access_data$ETHCAT, 
                             levels = c("White", "Asian", "Black", "Hispanic"))


set.seed(12)

model_formula <- KDPI ~ ETHCAT + PRA + AGE + HEP_C_ANTI_DON + NON_HRT_DON + 
  HIST_DIABETES_DON + AGE_DON + ETHCAT_DON + COD_CAD_DON + 
  CREAT_DON + HGT_CM_DON_CALC + WGT_KG_DON_CALC + BMI_CALC + 
  HIST_HYPERTENS_DON

# Set up 5-fold cross-validation
ctrl <- trainControl(method = "cv", number = 5)

# Train the model using cross-validation
cv_model <- train(model_formula, 
                  data = access_data, 
                  method = "rf", 
                  trControl = ctrl,
                  ntree = 500,
                  importance = TRUE)

# Get the mean R-squared across folds
mean_r_squared <- mean(cv_model$results$Rsquared)
print(paste("Mean R-squared across folds:", mean_r_squared))

# Feature importance
importanceScores <- varImp(cv_model)$importance
importanceDf <- data.frame(Feature = rownames(importanceScores), 
                           Importance = importanceScores$Overall)
importanceDf <- importanceDf[order(-importanceDf$Importance),]

# Plot feature importance
ggplot(importanceDf, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() +
  labs(x = "Feature", y = "Importance", 
       title = "Feature Importance (Increase in MSE when a Feature is Omitted)")

# Make predictions on the entire dataset
all_predictions <- predict(cv_model, newdata = access_data)

# Calculate overall R-squared
SSR <- sum((all_predictions - access_data$KDPI)^2)
SST <- sum((access_data$KDPI - mean(access_data$KDPI))^2)
overall_r_squared <- 1 - SSR/SST
print(paste("Overall R-squared:", overall_r_squared))

# Calculate other metrics
MAE <- mean(abs(all_predictions - access_data$KDPI))
MSE <- mean((all_predictions - access_data$KDPI)^2)
RMSE <- sqrt(MSE)

print(paste("MAE:", MAE))
print(paste("MSE:", MSE))
print(paste("RMSE:", RMSE))

# Create a copy of the entire dataset for counterfactual analysis
access_data_cf <- access_data

# Set all ethnicities to 'White' for the counterfactual scenario
access_data_cf$ETHCAT <- factor('White', 
                              levels = c("White", "Asian", "Black","Hispanic"))

# Predict KDPI scores for both actual and counterfactual scenarios
access_data$predicted_KDPI_actual <- predict(cv_model, newdata = access_data)
access_data$predicted_KDPI_cf <- predict(cv_model, newdata = access_data_cf)

# Compute the difference
access_data$difference <- access_data$predicted_KDPI_cf - access_data$predicted_KDPI_actual

# Summarize results
pred_data <- access_data %>% 
  group_by(ETHCAT) %>% 
  summarise(
    median_predicted_KDPI_Counterfactual = median(predicted_KDPI_cf),
    median_predicted_KDPI_Actual = median(predicted_KDPI_actual),
    median_difference = median(difference),
    var_difference = var(difference)
  )

print(pred_data)

# Kruskal-Wallis test
kruskal_result <- kruskal.test(difference ~ ETHCAT, data = access_data)
print(kruskal_result)

# Pairwise Wilcoxon tests
wilcox_result <- pairwise.wilcox.test(access_data$difference, access_data$ETHCAT,
                                      p.adjust.method = "BH")
print(wilcox_result)

################################################################################
################################################################################
################################################################################



################################################################################
############################## OUTCOME ANALYSIS ################################
################################################################################
# Uncomment if not running data prep section code
#balanced_data <- read.csv("balanced_data.csv") 

library(survival)
library(cmprsk)
library(dplyr)
library(purrr)
library(patchwork)
library(ggplot2)
library(tidyr)
library(dplyr)

outcome_data <- dplyr::select(balanced_data, 
                              GRF_FAIL_CAUSE_TY_KI, GSTATUS_KI, GTIME_KI, 
                              ETHCAT, AGE, GENDER, PRA, ON_DIALYSIS, ABO, 
                              KDPI, AGE_DON, BMI_DON_CALC, CREAT_DON, BMI)

outcome_data <- outcome_data %>%
  mutate(GRF_FAIL_CAUSE_TY_KI = case_when(
    GRF_FAIL_CAUSE_TY_KI == "Hyperacute Rejection" ~ 1,
    GRF_FAIL_CAUSE_TY_KI == "Acute Rejection" ~ 2,
    GRF_FAIL_CAUSE_TY_KI == "Primary Failure" ~ 3,
    GRF_FAIL_CAUSE_TY_KI == "Graft Thrombosis" ~ 4,
    GRF_FAIL_CAUSE_TY_KI == "Infection" ~ 5,
    GRF_FAIL_CAUSE_TY_KI == "Surgical Complications" ~ 6,
    GRF_FAIL_CAUSE_TY_KI == "Urological Complications" ~ 7,
    GRF_FAIL_CAUSE_TY_KI == "Recurrent Disease" ~ 8,
    GRF_FAIL_CAUSE_TY_KI == "Primary Non-Function (Graft Never Functioned Post-Transplant)" ~ 9,
    GRF_FAIL_CAUSE_TY_KI == "Chronic Rejection" ~ 10,
    GRF_FAIL_CAUSE_TY_KI == "BK (Polyoma) Virus" ~ 11,
    GRF_FAIL_CAUSE_TY_KI == "Primary Non-Function (Graft Never Functioned Post-Transplant)" ~ 12,
    GRF_FAIL_CAUSE_TY_KI == "Other" ~ 999,
    TRUE ~ NA_real_ 
  ))

# Recode graft failure causes
outcome_data <- outcome_data %>%
  mutate(GRF_FAIL_CAUSE = case_when(
    GRF_FAIL_CAUSE_TY_KI %in% c(1,2,10) ~ 1,
    !is.na(GRF_FAIL_CAUSE_TY_KI) ~ 2,
    TRUE ~ 0
  ))

# Create the status variable for competing risks
outcome_data <- outcome_data %>%
  mutate(CR_STATUS = case_when(
    GSTATUS_KI == 0 ~ 0,  # No graft failure
    GSTATUS_KI == 1 & GRF_FAIL_CAUSE == 1 ~ 1, # Graft failure due to rejection
    GSTATUS_KI == 1 & GRF_FAIL_CAUSE == 2 ~ 2  # Graft failure from other causes
  ))

outcome_data$ETHCAT <- relevel(factor(outcome_data$ETHCAT), ref = "White")

# Fit the competing risks model
cr_model <- cmprsk::crr(ftime = outcome_data$GTIME_KI,
                        fstatus = outcome_data$CR_STATUS,
                        cov1 = model.matrix(~ ETHCAT, 
                                            data = outcome_data)[,-1],
                        failcode = 1,  # 1 is the event of interest (rejection)
                        cencode = 0)   # 0 is censoring

summary(cr_model)




# Calculate cumulative incidence functions
cif <- cmprsk::cuminc(ftime = outcome_data$GTIME_KI,
                      fstatus = outcome_data$CR_STATUS,
                      group = outcome_data$ETHCAT)

# Function to extract data for each group and event type
extract_cif_data <- function(cif_obj, group, event) {
  data.frame(
    time = cif_obj[[paste(group, event, sep = " ")]]$time,
    est = cif_obj[[paste(group, event, sep = " ")]]$est,
    group = group,
    event = ifelse(event == "1", "Graft Rejection", "Competing Risk")
  )
}

# Create a list of all group and event combinations
groups <- c("White", "Asian", "Black", "Hispanic")
events <- c("1", "2")
combinations <- expand.grid(group = groups, event = events)

# Apply the function to all combinations and combine the results
cif_long <- purrr::pmap_dfr(combinations, ~extract_cif_data(cif, ..1, ..2))

cif_long <- cif_long %>%
  mutate(linetype = ifelse(event == "Graft Rejection", "solid", "dashed"))

y_limits <- range(cif_long$est)


cif_long$group <- factor(cif_long$group, 
                         levels = c("Asian", "Black", "Hispanic", "White"))

color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4", 
                   "Hispanic" = "#6baed6", "White" = "#b3cde3")

line_types <- c("Asian" = "solid", "Black" = "dashed", 
                "Hispanic" = "dotdash", "White" = "dotted")

# Create the plot for Graft Rejection
p1 <- ggplot(cif_long %>% filter(event == "Graft Rejection"), 
             aes(x = time, y = est, color = group, linetype = group)) +
  geom_step(size = 1.5) +  # Increase line thickness for better visibility
  scale_color_manual(values = color_palette, name = "Group") +
  scale_linetype_manual(values = line_types, name = "Group") +
  labs(title = "Graft Rejection", x = "Time (Days)", y = "Cumulative Incidence") +
  ylim(y_limits) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.key.width = unit(2, "cm"),  
    legend.key.height = unit(0.6, "cm"),  
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

# Create the plot for Competing Risk
p2 <- ggplot(cif_long %>% filter(event == "Competing Risk"), 
             aes(x = time, y = est, color = group, linetype = group)) +
  geom_step(size = 1.5) +  
  scale_color_manual(values = color_palette, name = "Group") +
  scale_linetype_manual(values = line_types, name = "Group") +
  labs(title = "Competing Risks", x = "Time (Days)", y = "Cumulative Incidence") +
  ylim(y_limits) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.key.width = unit(2, "cm"), 
    legend.key.height = unit(0.6, "cm"), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

# Combine the plots side by side with a shared legend at the bottom
combined_plot <- p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

# Display the combined plot
print(combined_plot)

ggsave("cumulative_incidence_plot.png", plot = combined_plot, width = 10, 
       height = 6, dpi = 300)

################################################################################
################################################################################
################################################################################




################################################################################
###################### SIMULATION FOR SUPPLEMENTARY MATERIALS ##################
################################################################################

library(dplyr)
library(purrr)
library(MASS)
library(ggplot2)
library(dplyr)

# Uncomment if not running data prep section code
# treatment_effect <- read.csv("treatment_effect_data.csv")

# Functions for computing DEA Efficiency Scores
define_io <- function(data) {
  list(
    inputs = as.matrix(data[, c("WaitlistDuration", "QualityScore")]),
    outputs = as.matrix(data[, "OutcomeScore"])
  )
}

compute_efficiency <- function(data) {
  io <- define_io(data)
  eff <- dea(X = io$inputs, Y = io$outputs, 
             RTS = "vrs", ORIENTATION = "graph")$eff
  data$efficiency <- eff
  return(data)
}


#### Computing U.S. Population proportions ####
# 2020 Census (https://data.census.gov/table?g=010XX00US) 
Asian_n <- 19886049
Black_n <- 41104200
Hispanic_n <- 62080044
White_n <- 204277273
Eth_total <- Asian_n + Black_n + Hispanic_n + White_n
Asian_prop <- Asian_n/Eth_total
Black_prop <- Black_n/Eth_total
Hispanic_prop <- Hispanic_n/Eth_total
White_prop <- White_n/Eth_total
desired_props <- c(Asian = Asian_prop, Black = Black_prop, Hispanic = Hispanic_prop, White = White_prop)

# Function to perform stratified sampling
stratified_sample <- function(data, group_var, desired_props, total_samples) {
  desired_counts <- floor(desired_props * total_samples)
  grouped_data <- group_split(data, !!sym(group_var))
  sampled_data <- map2_dfr(grouped_data, names(desired_counts), 
                           function(group_data, group_name) {
    slice_sample(group_data, n = desired_counts[group_name])})
  return(sampled_data)
}

treatment_effect_sim <- treatment_effect %>%
  mutate(Group = case_when(
    Group == "1" ~ "Asian",
    Group == "2" ~ "Black",
    Group == "3" ~ "Hispanic",
    Group == "4" ~ "White",
    TRUE ~ as.character(Group)
  ))

groups <- c("Asian", "Black", "Hispanic", "White")

# Computing covariance matrices for each group
cov_matrices <- list()

# Loop through each group and compute the covariance matrix
for (group in groups) {
  # Subset the data for the group
  group_data <- treatment_effect_sim %>%
    filter(Group == group) %>%
    dplyr::select(WaitlistDuration, QualityScore, OutcomeScore)
  cov_matrix <- cov(group_data)
  cov_matrices[[group]] <- cov_matrix
}


data_summary <- treatment_effect_sim %>% 
  group_by(Group) %>% 
  summarise(
    wait_mean = mean(WaitlistDuration),
    wait_median = median(WaitlistDuration),
    wait_sd = sd(WaitlistDuration),
    quality_mean = mean(QualityScore),
    quality_median = median(QualityScore),
    quality_sd = sd(QualityScore),
    outcome_mean = mean(OutcomeScore),
    outcome_median = median(OutcomeScore),
    outcome_sd = sd(OutcomeScore)
  )

#################################### SIMULATION ################################

# Generate synthetic data function
generate_synthetic_data <- function(n = 1000, summary_stats = data_summary, cov_matrices) {
  groups <- sample(c("Asian", "Black", "Hispanic", "White"), n, replace = TRUE,
                   prob = c(0.078, 0.289, 0.130, 0.502))  # Original proportions
  
  data <- lapply(unique(groups), function(group) {
    group_n <- sum(groups == group)
    stats <- summary_stats %>% filter(Group == group)
    means <- c(stats$wait_mean, stats$quality_mean, stats$outcome_mean)
    cov_matrix <- cov_matrices[[group]]
    synthetic_data <- MASS::mvrnorm(n = group_n, mu = means, Sigma = cov_matrix)
    synthetic_data <- pmax(synthetic_data, 0)
    tibble(
      Group = group,
      WaitlistDuration = synthetic_data[, 1],
      QualityScore = synthetic_data[, 2],
      OutcomeScore = synthetic_data[, 3]
    )
  })
  bind_rows(data)
}

# Simulation function
run_simulation <- function(n_iterations = 100) {
  results <- tibble()
  
  for (i in 1:n_iterations) {
    # Generate synthetic data
    synthetic_data <- generate_synthetic_data(
      n = 1000, summary_stats = data_summary, cov_matrices = cov_matrices
      )
    
    original_efficiency <- compute_efficiency(synthetic_data)
    
    # Determine the total sample size for stratified sampling
    group_sizes <- table(synthetic_data$Group)
    max_possible_samples <- min(group_sizes / desired_props[names(group_sizes)])
    total_samples <- floor(max_possible_samples)
    
    # Perform stratified sampling
    resampled_data <- stratified_sample(synthetic_data, "Group", 
                                        desired_props, total_samples)
    
    resampled_efficiency <- compute_efficiency(resampled_data)
    
    # Store results
    results <- bind_rows(results,
                         original_efficiency %>% 
                           mutate(Method = "Imbalanced") %>%
                           dplyr::select(Group, efficiency, Method),
                         resampled_efficiency %>% 
                           mutate(Method = "Resampled (U.S.)") %>%
                           dplyr::select(Group, efficiency, Method))
  }
  
  return(results)
}

# Run simulation
set.seed(123)
sim_results <- run_simulation()

# Analyze results
summary_stats <- sim_results %>%
  group_by(Group, Method) %>%
  summarise(
    Mean = mean(efficiency),
    SD = sd(efficiency),
    CI_lower = quantile(efficiency, 0.025),
    CI_upper = quantile(efficiency, 0.975)
  )

print(summary_stats)


ggplot(sim_results, aes(x = Group, y = efficiency, fill = Method)) +
  geom_violin(trim = T, alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), 
               outlier.size = 1, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    y = "Efficiency Score",
    x = "Ethnic Group"
  )

ggsave("efficiency_comparison.png", width = 10, height = 6)





############################### SENSITIVITY ANALYSIS ###########################

proportions_list <- list(
  c(Asian = 0.1, Black = 0.2, Hispanic = 0.3, White = 0.4), # Alternative 1
  c(Asian = 0.2, Black = 0.4, Hispanic = 0.3, White = 0.1), # Alternative 2
  c(Asian = 0.078, Black = 0.289, Hispanic = 0.130, White = 0.502) # Original proportions
)

# Generate synthetic data function (modified to use original proportions)
generate_synthetic_data <- function(n = 1000, summary_stats = data_summary, cov_matrices) {
  groups <- sample(c("Asian", "Black", "Hispanic", "White"), n, replace = TRUE,
                   prob = c(Asian = 0.061, Black = 0.126,
                            Hispanic = 0.190, White = 0.624))  # US proportions
  
  data <- lapply(unique(groups), function(group) {
    group_n <- sum(groups == group)
    stats <- summary_stats %>% filter(Group == group)
    means <- c(stats$wait_mean, stats$quality_mean, stats$outcome_mean)
    cov_matrix <- cov_matrices[[group]]
    synthetic_data <- MASS::mvrnorm(n = group_n, mu = means, Sigma = cov_matrix)
    synthetic_data <- pmax(synthetic_data, 0)
    tibble(
      Group = group,
      WaitlistDuration = synthetic_data[, 1],
      QualityScore = synthetic_data[, 2],
      OutcomeScore = synthetic_data[, 3]
    )
  })
  bind_rows(data)
}

run_sensitivity_analysis <- function(n_iterations = 10, proportions_list) {
  results <- tibble()
  
  us_proportions <- c(Asian = 0.061, Black = 0.126, 
                      Hispanic = 0.190, White = 0.624)
  
  for (proportions in proportions_list) {
    for (i in 1:n_iterations) {
      # Generate synthetic data using US proportions
      synthetic_data <- generate_synthetic_data(n = 1000, 
                                                summary_stats = data_summary, 
                                                cov_matrices = cov_matrices)
      
      # Compute efficiency for US proportions data
      us_efficiency <- compute_efficiency(synthetic_data)
      
      # Resample and compute efficiency for alternative proportions
      resampled_data <- stratified_sample(synthetic_data, "Group", 
                                          proportions, nrow(synthetic_data))
      resampled_efficiency <- compute_efficiency(resampled_data)
      
      # Store results
      results <- bind_rows(results,
                           us_efficiency %>% 
                             mutate(Method = "US Proportions", Proportions = "US") %>%
                             dplyr::select(Group, efficiency, Method, Proportions),
                           resampled_efficiency %>% 
                             mutate(Method = "Alternative", Proportions = paste(proportions, collapse = "-")) %>%
                             dplyr::select(Group, efficiency, Method, Proportions))
    }
  }
  return(results)
}

# Run the sensitivity analysis
set.seed(123)
sensitivity_results <- run_sensitivity_analysis(n_iterations = 100, proportions_list = proportions_list)

# Analyze results
sensitivity_summary <- sensitivity_results %>%
  group_by(Group, Method, Proportions) %>%
  summarise(
    Mean = mean(efficiency),
    SD = sd(efficiency),
    CI_lower = quantile(efficiency, 0.025),
    CI_upper = quantile(efficiency, 0.975), 
    n()
  )
print(sensitivity_summary, n = 24)

proportion_order <- c("US", "0.078-0.289-0.13-0.502", "0.1-0.2-0.3-0.4", "0.2-0.4-0.3-0.1")

proportions_labels <- c(
  "US" = "US Population",
  "0.078-0.289-0.13-0.502" = "Original Imbalanced Data",
  "0.1-0.2-0.3-0.4" = "Alternative 1",
  "0.2-0.4-0.3-0.1" = "Alternative 2"
)

# Calculate mean efficiency for each group and proportion
mean_efficiency <- sensitivity_results %>%
  group_by(Group, Proportions) %>%
  summarise(mean_efficiency = mean(efficiency), .groups = 'drop')

color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4", "Hispanic" = "#6baed6", "White" = "#b3cde3")

ggplot(mean_efficiency, aes(x = factor(Proportions, levels = proportion_order), y = mean_efficiency, color = Group, group = Group)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_discrete(labels = proportions_labels) +
  theme_minimal(base_size = 14) + 
  scale_color_manual(values = color_palette) + 
  theme(
    strip.text = element_text(size = 14, face = "bold"), 
    legend.position = "top",
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12), 
    plot.title = element_text(size = 16, face = "bold"), 
    plot.subtitle = element_text(size = 14) 
  ) +
  labs(
    y = "Mean Efficiency Score",
    x = "Proportion Scenarios"
  )
ggsave("sensitivity_efficiency_undersampling_stratified.png", width = 10, height = 6, dpi = 300)


################################## BIAS ANALYSIS ###############################
# Calculate mean efficiency for each group and proportion scenario
mean_efficiencies <- sensitivity_results %>%
  group_by(Group, Method, Proportions) %>%
  summarise(Mean_Efficiency = mean(efficiency), .groups = 'drop')

# Calculate bias relative to US proportions
bias_calculation <- mean_efficiencies %>%
  pivot_wider(names_from = Method, values_from = Mean_Efficiency) %>%
  group_by(Group) %>%
  mutate(
    US_Efficiency = `US Proportions`[Proportions == "US"],
    Bias = case_when(
      Proportions == "US" ~ 0,
      TRUE ~ Alternative - US_Efficiency
    )
  ) %>%
  dplyr::select(Group, Proportions, `US Proportions`, Alternative, Bias)

weighted_bias <- bias_calculation %>%
  filter(Proportions != "US") %>% 
  rowwise() %>%
  mutate(
    Group_Proportion = case_when(
      Group == "Asian" ~ as.numeric(word(Proportions, 1, sep = "-")),
      Group == "Black" ~ as.numeric(word(Proportions, 2, sep = "-")),
      Group == "Hispanic" ~ as.numeric(word(Proportions, 3, sep = "-")),
      Group == "White" ~ as.numeric(word(Proportions, 4, sep = "-"))
    )
  ) %>%
  group_by(Proportions) %>%
  summarise(Weighted_Bias = sum(Bias * Group_Proportion, na.rm = TRUE))

print(weighted_bias)
################################################################################
################################################################################
################################################################################
