library(mma)
library(lmtest)
require(nnet) 

DEA_data <- read.csv("DEA_data.csv")

priority_data <- DEA_data %>%
  dplyr::select(
    Priority_Score,
    Group,
    REGION, EDUCATION, CITIZENSHIP,
    PRI_PAYMENT_TCR_KI,
    DISTANCE, WORK_INCOME_TCR,
    ON_DIALYSIS, ABO,
    PRA, AGE, GENDER,
    PREV_KI_TX, MED_COND_TRR
  ) %>%
  dplyr::mutate(
    Group = stats::relevel(factor(Group), ref = "White"),
    dplyr::across(
      c(
        REGION,
        EDUCATION,
        CITIZENSHIP,
        PRI_PAYMENT_TCR_KI,
        WORK_INCOME_TCR,
        ON_DIALYSIS,
        ABO,
        GENDER,
        MED_COND_TRR,
        PREV_KI_TX
      ),
      ~ factor(.)
    )
  )


drop_all_levels <- function(data) {
  # Identify which columns are factors or characters
  col_classes <- sapply(data, class)
  categorical_cols <- names(data)[col_classes %in% c("factor", "character")]
  
  # For each categorical column, convert to factor and drop unused levels
  for (col in categorical_cols) {
    if (is.character(data[[col]])) {
      data[[col]] <- factor(data[[col]])
    }
    data[[col]] <- droplevels(data[[col]])
  }
  
  return(data)
}

priority_data <- drop_all_levels(priority_data)

# ################################################################################
# #                               Mediator/Covariate Testing                     #
# ################################################################################
# 
# ########## 1) Test for each potential Mediator (M) and Response (Y): ###########

candid_med_cov = names(priority_data[,which(!(colnames(priority_data) %in% c("Group", "Priority_Score")))])

# Assuming priority_data is already loaded
set.seed(2)

# Create full model
full_model <- glm(Priority_Score ~ ., data = priority_data)

# Get all predictor names
predictors <- candid_med_cov #names(coef(full_model))[-1]  # Exclude intercept

# Initialize results dataframe
table_M_Y <- data.frame(
  Term = character(0),
  LR_Chisq = numeric(0), 
  Df = numeric(0),
  `Pr(>Chisq)` = numeric(0)
)

# Full model log-likelihood
ll_full <- logLik(full_model)

# Get LR test for each predictor
for (pred in predictors) {
  # Create formula for reduced model
  reduced_formula <- paste("Priority_Score ~ . -", pred)
  
  # Fit reduced model
  reduced_model <- glm(as.formula(reduced_formula), data = priority_data)
  
  # Get log-likelihood for reduced model
  ll_reduced <- logLik(reduced_model)
  
  # Calculate LR statistic
  lr_stat <- -2 * (as.numeric(ll_reduced) - as.numeric(ll_full))
  
  # Degrees of freedom
  df <- attr(ll_full, "df") - attr(ll_reduced, "df")
  
  # Calculate p-value
  p_value <- 1 - pchisq(lr_stat, df = df)
  
  # Add to results dataframe
  table_M_Y <- rbind(table_M_Y, data.frame(
    Term = pred,
    LR_Chisq = lr_stat,
    Df = df,
    `Pr(>Chisq)` = p_value
  ))
}

# Add adjusted p-values and significance flags
table_M_Y$Adj_P <- p.adjust(table_M_Y$`Pr..Chisq.`, method = "BH")
table_M_Y$Significant <- ifelse(as.numeric(table_M_Y$Adj_P) <= 0.05, TRUE, FALSE)

# Display the table
table_M_Y







############# 2) Test for each potential Mediator (M) and Predictor (X): #######
mediator_names <- candid_med_cov
########################################################
#          Using Log likelihood ratio testing          #
########################################################

test_M_X_significance <- function(data, mediator_names, group_var = "Group") {
  # data: dataframe containing mediators and group variable
  # mediator_names: character vector of mediator column names
  # group_var: name of the grouping variable (default: "Group")
  
  
  results <- data.frame(
    mediator = mediator_names,
    p_value = NA,
    test_type = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(mediator_names)) {
    med_name <- mediator_names[i]
    
    # Skip if mediator has too many missing values
    if (sum(!is.na(data[[med_name]])) <= 1) {
      next
    }
    
    # Create formula
    formula_med <- as.formula(paste(med_name, "~", group_var))
    formula_null <- as.formula(paste(med_name, "~ 1"))

    # Determine mediator type and run appropriate test
    if (is.numeric(data[[med_name]])) {
      model_full <- tryCatch(
        glm(formula_med, data = data),
        error = function(e) NULL
      )
      model_null <- tryCatch(
        glm(formula_null, data = data),
        error = function(e) NULL
      )
      if (!is.null(model_full) && !is.null(model_null)) {
        results$p_value[i] <- lrtest(model_full,model_null)$`Pr(>Chisq)`[2]
        results$test_type[i] <- "Likelihood ratio test"
      }
      
    } else if (nlevels(as.factor(data[[med_name]])) == 2) {
      # Binary mediator: Likelihood ratio test with logistic regression
      model_full <- tryCatch(
        glm(formula_med, data = data, family = binomial),
        error = function(e) NULL
      )
      model_null <- tryCatch(
        glm(formula_null, data = data, family = binomial),
        error = function(e) NULL
      )
      
      
      if (!is.null(model_full) && !is.null(model_null)) {
        results$p_value[i] <- lrtest(model_full,model_null)$`Pr(>Chisq)`[2]

        results$test_type[i] <- "Likelihood ratio test"
      }
      
    } else {
      # Categorical mediator with >2 levels: Likelihood ratio test with multinomial
      model_full <- tryCatch(
        multinom(formula_med, data = data, trace = FALSE),
        error = function(e) NULL
      )
      model_null <- tryCatch(
        multinom(formula_null, data = data, trace = FALSE),
        error = function(e) NULL
      )
      if (!is.null(model_full) && !is.null(model_null)) {
        results$p_value[i] <- lrtest(model_full, model_null)$`Pr(>Chisq)`[2]
        results$test_type[i] <- "Multinomial likelihood ratio test"
      }
    }
  }
  return(results)
}

M_X_tests <- test_M_X_significance(priority_data, mediator_names, group_var = "Group")
M_X_tests$adj_p_value <- p.adjust(M_X_tests$p_value, method = "BH")

Adjsusted_pval_results <- data.frame(Candidate = candid_med_cov, M_Y = table_M_Y$Adj_P, M_X = M_X_tests$adj_p_value)

(Adjsusted_pval_results <- Adjsusted_pval_results %>% 
  mutate(Variable = case_when(
    M_Y < 0.05 & M_X < 0.05 ~ "Mediator",
    M_Y < 0.05 & M_X >= 0.05 ~ "Covariate",
    TRUE ~ "Discard"
  )))

(final_mediators <- (Adjsusted_pval_results[Adjsusted_pval_results$Variable == "Mediator",])$Candidate)
(final_covariates <- (Adjsusted_pval_results[Adjsusted_pval_results$Variable == "Covariate",])$Candidate)

################################################################################
######################### Conduct Mediation Analysis ###########################
################################################################################

x = priority_data[,which(colnames(priority_data) %in% final_mediators)]
pred = priority_data[,2]
y = priority_data[,1]

# Residualize the response (Priority_Score) on covariates
resid_model <- lm(as.formula(paste("Priority_Score ~", paste(final_covariates, collapse = "+"))), data = priority_data)
y_resid <- resid(resid_model)

# Fit Model
set.seed(4)
mma_model <- mma(x, y_resid, mediator = 1:ncol(x),  pred=pred, predref="White", alpha=0.05,alpha2=0.05, testtype = 1, n2 = 1000, para = F)
summary(mma_model)


boot_samples_Asian_indirect <- mma_model[["a.binx"]][["bootsresults"]][["ie"]][["pred"]][,1]
boot_samples_Black_indirect <- mma_model[["a.binx"]][["bootsresults"]][["ie"]][["pred.Black"]][,1]
boot_samples_Hispanic_indirect <- mma_model[["a.binx"]][["bootsresults"]][["ie"]][["pred.Hispanic"]][,1]

boot_samples_Asian_direct <- mma_model[["a.binx"]][["bootsresults"]][["de"]][,1]
boot_samples_Black_direct <- mma_model[["a.binx"]][["bootsresults"]][["de"]][,2]
boot_samples_Hispanic_direct <- mma_model[["a.binx"]][["bootsresults"]][["de"]][,3]

boot_samples_Asian_total <- mma_model[["a.binx"]][["bootsresults"]][["te"]][,1]
boot_samples_Black_total <- mma_model[["a.binx"]][["bootsresults"]][["te"]][,2]
boot_samples_Hispanic_total <- mma_model[["a.binx"]][["bootsresults"]][["te"]][,3]


# Extract statistics for all effects
calculate_stats <- function(boot_sample) {
  est_mean <- mean(boot_sample)
  est_se <- sd(boot_sample)
  ci_lower <- as.numeric(quantile(boot_sample, 0.025))
  ci_upper <- as.numeric(quantile(boot_sample, 0.975))
  z_val <- est_mean / est_se
  p_val <- 2 * (1 - pnorm(abs(z_val)))
  
  return(c(Estimate = est_mean, 
           SE = est_se, 
           CI_Lower = ci_lower, 
           CI_Upper = ci_upper, 
           P_Value = p_val))
}

# Apply function to each set of boot samples
Asian_indirect_stats <- calculate_stats(boot_samples_Asian_indirect)
Black_indirect_stats <- calculate_stats(boot_samples_Black_indirect)
Hispanic_indirect_stats <- calculate_stats(boot_samples_Hispanic_indirect)

Asian_direct_stats <- calculate_stats(boot_samples_Asian_direct)
Black_direct_stats <- calculate_stats(boot_samples_Black_direct)
Hispanic_direct_stats <- calculate_stats(boot_samples_Hispanic_direct)

Asian_total_stats <- calculate_stats(boot_samples_Asian_total)
Black_total_stats <- calculate_stats(boot_samples_Black_total)
Hispanic_total_stats <- calculate_stats(boot_samples_Hispanic_total)


results_df <- data.frame(
  Race = rep(c("Asian", "Black", "Hispanic"), 3),
  Effect_Type = c(rep("Indirect", 3), rep("Direct", 3), rep("Total", 3)),
  Estimate = c(
    Asian_indirect_stats["Estimate"], Black_indirect_stats["Estimate"], Hispanic_indirect_stats["Estimate"],
    Asian_direct_stats["Estimate"], Black_direct_stats["Estimate"], Hispanic_direct_stats["Estimate"],
    Asian_total_stats["Estimate"], Black_total_stats["Estimate"], Hispanic_total_stats["Estimate"]
  ),
  SE = c(
    Asian_indirect_stats["SE"], Black_indirect_stats["SE"], Hispanic_indirect_stats["SE"],
    Asian_direct_stats["SE"], Black_direct_stats["SE"], Hispanic_direct_stats["SE"],
    Asian_total_stats["SE"], Black_total_stats["SE"], Hispanic_total_stats["SE"]
  ),
  CI_Lower = c(
    Asian_indirect_stats["CI_Lower"], Black_indirect_stats["CI_Lower"], Hispanic_indirect_stats["CI_Lower"],
    Asian_direct_stats["CI_Lower"], Black_direct_stats["CI_Lower"], Hispanic_direct_stats["CI_Lower"],
    Asian_total_stats["CI_Lower"], Black_total_stats["CI_Lower"], Hispanic_total_stats["CI_Lower"]
  ),
  CI_Upper = c(
    Asian_indirect_stats["CI_Upper"], Black_indirect_stats["CI_Upper"], Hispanic_indirect_stats["CI_Upper"],
    Asian_direct_stats["CI_Upper"], Black_direct_stats["CI_Upper"], Hispanic_direct_stats["CI_Upper"],
    Asian_total_stats["CI_Upper"], Black_total_stats["CI_Upper"], Hispanic_total_stats["CI_Upper"]
  ),
  P_Value = c(
    Asian_indirect_stats["P_Value"], Black_indirect_stats["P_Value"], Hispanic_indirect_stats["P_Value"],
    Asian_direct_stats["P_Value"], Black_direct_stats["P_Value"], Hispanic_direct_stats["P_Value"],
    Asian_total_stats["P_Value"], Black_total_stats["P_Value"], Hispanic_total_stats["P_Value"]
  )
)

results_df$CI_95 <- paste0("[", round(results_df$CI_Lower, 3), ", ", round(results_df$CI_Upper, 3), "]")

(final_table <- results_df[, c("Race", "Effect_Type", "Estimate", "SE", "CI_95", "P_Value")])
