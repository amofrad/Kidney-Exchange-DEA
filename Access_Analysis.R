library(caret)
library(randomForest)
library(dplyr)

balanced_data <- read.csv("balanced_data.csv")

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