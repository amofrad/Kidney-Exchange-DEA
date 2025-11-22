library(caret)
library(randomForest)
library(dplyr)


DEA_data <- read.csv("DEA_data.csv")


access_data <- DEA_data %>%
  dplyr::select(
    LKDPI, Access_Score, ETHCAT, eGFR, AGE_DON, HIST_CIG_DON,
    HLA_B_mismatches, HLA_DR_mismatches,
    SBP, GENDER, GENDER_DON, ABO_MAT,
    ABO_incompatible,
    B1, DB1, B2, DB2, DDR1, DR1, DDR2, DR2,
    WGT_KG_CALC, WGT_KG_DON_CALC,
    WL_days, ABO, ETHCAT_DON,
    BMI_DON_CALC, BMI_CALC, AGE,
    REGION, EDUCATION, CITIZENSHIP,
    PRI_PAYMENT_TCR_KI, PREV_KI_TX,
    DISTANCE, WORK_INCOME_TCR, Year
  ) %>%
  dplyr::mutate(
    Donor_Black = factor(ifelse(ETHCAT_DON == "Black", 1, 0)),
    Don_Rec_WR  = WGT_KG_DON_CALC / WGT_KG_CALC,
    Don_Rec_WR  = ifelse(Don_Rec_WR >= 0.9, 0.9, Don_Rec_WR)
  ) %>%
  dplyr::group_by(Year) %>%
  dplyr::mutate(
    LKDPI = LKDPI - mean(LKDPI)   # same as your original (no na.rm)
  ) %>%
  dplyr::ungroup()

access_data %>% group_by(ETHCAT) %>% summarise(mean(LKDPI), median(LKDPI))


#############################################################
# Training the model
#############################################################
# Model to Train
model_formula <- LKDPI ~ ETHCAT + Donor_Black + BMI_DON_CALC + eGFR + AGE_DON + HIST_CIG_DON + SBP + GENDER + GENDER_DON + #ABO_MAT + 
  HLA_B_mismatches + HLA_DR_mismatches + Don_Rec_WR + ABO_incompatible + Year + REGION + EDUCATION + CITIZENSHIP + PRI_PAYMENT_TCR_KI + PREV_KI_TX + WORK_INCOME_TCR

# Set up 5-fold cross-validation
ctrl <- trainControl(method = "cv", number = 5)

set.seed(12)

tuneGrid = expand.grid(
  mtry = seq(1, 15, 1),
  splitrule = "variance",
  min.node.size = seq(1, 10, 1)
)

cv_model <- train(
  model_formula,
  data = access_data,
  method = "ranger",
  trControl = ctrl,
  tuneGrid = tuneGrid,
  num.trees = 500,
  importance = "permutation"
)

cv_model

cv_model$finalModel$prediction.error

# Feature importance
importanceScores <- varImp(cv_model)$importance
importanceDf <- data.frame(Feature = rownames(importanceScores), 
                           Importance = importanceScores$Overall)
importanceDf <- importanceDf[order(-importanceDf$Importance),]

importanceDf

# Plot feature importance
ggplot(importanceDf, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() +
  labs(x = "Feature", y = "Importance", 
       title = "Feature Importance")

# Make predictions on the entire dataset
all_predictions <- predict(cv_model, newdata = access_data)

# Compute R-squared
SSR <- sum((all_predictions - access_data$LKDPI)^2)
SST <- sum((access_data$LKDPI - mean(access_data$LKDPI))^2)
overall_r_squared <- 1 - SSR/SST
print(paste("Overall R-squared:", overall_r_squared))
# Compute other metrics
(MAE <- mean(abs(all_predictions - access_data$LKDPI)))
(MSE <- mean((all_predictions - access_data$LKDPI)^2))
(RMSE <- sqrt(MSE))

# Create a copy of the dataset for counterfactual analysis
access_data_cf <- access_data

# Set all ethnicities to 'White' for the counterfactual scenario
access_data_cf$ETHCAT <- factor('White', 
                                levels = c("White", "Asian", "Black","Hispanic"))

# Predict LKDPI scores for both actual and counterfactual scenarios
access_data$predicted_LKDPI_actual <- predict(cv_model, newdata = access_data)
access_data$predicted_LKDPI_cf <- predict(cv_model, newdata = access_data_cf)

# Compute the difference
access_data$difference <- access_data$predicted_LKDPI_cf - access_data$predicted_LKDPI_actual

access_data %>% 
  group_by(ETHCAT) %>% 
  summarise(
    mean_LKDPI_Actual = formatC(mean(predicted_LKDPI_actual), format = "f", digits = 4),
    mean_LKDPI_CF = formatC(mean(predicted_LKDPI_cf), format = "f", digits = 4),
    mean_diff = formatC(mean(difference), format = "f", digits = 4),
    var_diff = var(difference)
  )

# What percentage of Patients saw an improvement (lower LKDPI) by only changing Ethnicity to White
access_data %>% 
  group_by(ETHCAT) %>% 
  summarise(mean(predicted_LKDPI_actual > predicted_LKDPI_cf))


# Kruskal-Wallis test
kruskal.test(difference ~ ETHCAT, data = access_data)$p.value

# Pairwise Wilcoxon tests
pairwise.wilcox.test(access_data$difference, access_data$ETHCAT,
                     p.adjust.method = "BH")



