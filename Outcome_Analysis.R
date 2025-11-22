library(survival)
library(cmprsk)
library(dplyr)
library(riskRegression)
library(prodlim)
library(ggplot2)
library(tidyr)
library(patchwork)

DEA_data <- read.csv("DEA_data.csv")

# Build outcome_data
outcome_data <- DEA_data %>%
  dplyr::select(
    GTIME_KI, GRF_FAIL_CAUSE_TY_KI, GSTATUS_KI,
    ETHCAT, AGE, GENDER, PRA, ON_DIALYSIS, ABO,
    LKDPI, AGE_DON, BMI_DON_CALC, BMI_CALC,
    REGION, EDUCATION, CITIZENSHIP,
    PRI_PAYMENT_TCR_KI, PREV_KI_TX,
    DISTANCE, WORK_INCOME_TRR,
    TRTREJ1Y_KI, TRTREJ6M_KI, FUNC_STAT_TRF, Year
  ) %>%
  # Center GTIME_KI within year and round
  dplyr::group_by(Year) %>%
  dplyr::mutate(
    GTIME_KI = round(GTIME_KI - mean(GTIME_KI))
  ) %>%
  dplyr::ungroup()

# Ensure GTIME_KI is strictly positive
min_time <- min(outcome_data$GTIME_KI)
if (min_time < 0) {
  outcome_data <- outcome_data %>%
    dplyr::mutate(GTIME_KI = GTIME_KI - min_time + 1)
}

# Recode GRF_FAIL_CAUSE_TY_KI to numeric codes
outcome_data <- outcome_data %>%
  dplyr::mutate(
    GRF_FAIL_CAUSE_TY_KI = dplyr::case_when(
      GRF_FAIL_CAUSE_TY_KI == "Hyperacute Rejection" ~ 1,
      GRF_FAIL_CAUSE_TY_KI == "Acute Rejection" ~ 2,
      GRF_FAIL_CAUSE_TY_KI == "Primary Failure" ~ 3,
      GRF_FAIL_CAUSE_TY_KI == "Graft Thrombosis" ~ 4,
      GRF_FAIL_CAUSE_TY_KI == "Infection" ~ 5,
      GRF_FAIL_CAUSE_TY_KI == "Surgical Complications" ~ 6,
      GRF_FAIL_CAUSE_TY_KI == "Urological Complications" ~ 7,
      GRF_FAIL_CAUSE_TY_KI == "Recurrent Disease" ~ 8,
      GRF_FAIL_CAUSE_TY_KI ==
        "Primary Non-Function (Graft Never Functioned Post-Transplant)" ~ 9,
      GRF_FAIL_CAUSE_TY_KI == "Chronic Rejection" ~ 10,
      GRF_FAIL_CAUSE_TY_KI == "BK (Polyoma) Virus" ~ 11,
      GRF_FAIL_CAUSE_TY_KI == "Primary Non-Function (Graft Never Functioned Post-Transplant)" ~ 12,
      GRF_FAIL_CAUSE_TY_KI == "Other" ~ 999,
      TRUE ~ NA_real_
    )
  ) %>%
  # Collapse causes into rejection vs other
  dplyr::mutate(
    GRF_FAIL_CAUSE = dplyr::case_when(
      GRF_FAIL_CAUSE_TY_KI %in% c(1, 2, 10) ~ 1,
      !is.na(GRF_FAIL_CAUSE_TY_KI)          ~ 2,
      TRUE                                  ~ 2
    )
  ) %>%
  # Competing risks status
  dplyr::mutate(
    CR_STATUS = dplyr::case_when(
      GSTATUS_KI == 0                      ~ 0L,  # No graft failure
      GSTATUS_KI == 1 & GRF_FAIL_CAUSE == 1 ~ 1L, # Failure due to rejection
      GSTATUS_KI == 1 & GRF_FAIL_CAUSE == 2 ~ 2L  # Failure from other causes
    )
  ) %>%
  # Factor conversions and reference level for ETHCAT
  dplyr::mutate(
    REGION        = factor(REGION),
    EDUCATION     = factor(EDUCATION),
    CITIZENSHIP   = factor(CITIZENSHIP),
    TRTREJ1Y_KI   = factor(TRTREJ1Y_KI),
    WORK_INCOME_TRR = factor(WORK_INCOME_TRR),
    ETHCAT        = stats::relevel(factor(ETHCAT), ref = "White")
  )

# Sum-to-zero contrasts for ETHCAT (4 levels)
contrasts(outcome_data$ETHCAT) <- contr.sum(4)
colnames(contrasts(outcome_data$ETHCAT)) <- c("Asian", "Black", "Hispanic")

# Fit the competing risks model
cs_model <- CSC(
  formula = Hist(GTIME_KI, CR_STATUS) ~ ETHCAT + TRTREJ1Y_KI + WORK_INCOME_TRR,
  data    = outcome_data
)


# Results (as summarized in Table 6)
(cs_model)


# Function to calculate corrected CIF for a given cause
calculate_corrected_cif <- function(cs_model, cause_num, outcome_data) {
  
  # Coefficients from the specific cause model
  cox_model <- cs_model[["models"]][[paste0("Cause ", cause_num)]]
  coefs <- coef(cox_model)
  
  # Extract the ethnicity coefficients
  asian_coef <- coefs["ETHCATAsian"]
  black_coef <- coefs["ETHCATBlack"] 
  hispanic_coef <- coefs["ETHCATHispanic"]
  white_coef <- -(asian_coef + black_coef + hispanic_coef)
  
  # Get modal values and coefficients for covariates
  modal_trtrej <- names(sort(table(outcome_data$TRTREJ1Y_KI), decreasing = TRUE))[1]
  modal_work <- names(sort(table(outcome_data$WORK_INCOME_TRR), decreasing = TRUE))[1]
  trtrej_coefs <- coefs[grep("TRTREJ1Y_KI", names(coefs))]
  work_coefs <- coefs[grep("WORK_INCOME_TRR", names(coefs))]
  
  # Calculate baseline linear predictor
  baseline_lp <- 0
  if(modal_trtrej != "N" && length(trtrej_coefs) > 0) {
    trtrej_coef_name <- paste0("TRTREJ1Y_KI", modal_trtrej)
    if(trtrej_coef_name %in% names(trtrej_coefs)) {
      baseline_lp <- baseline_lp + trtrej_coefs[trtrej_coef_name]
    }
  }
  if(modal_work != "N" && length(work_coefs) > 0) {
    work_coef_name <- paste0("WORK_INCOME_TRR", modal_work)
    if(work_coef_name %in% names(work_coefs)) {
      baseline_lp <- baseline_lp + work_coefs[work_coef_name]
    }
  }
  
  # Calculate hazard ratios for each ethnicity
  hr_asian <- exp(baseline_lp + asian_coef)
  hr_black <- exp(baseline_lp + black_coef)
  hr_hispanic <- exp(baseline_lp + hispanic_coef)
  hr_white <- exp(baseline_lp + white_coef)
  
  # Create survival object for this specific cause
  surv_obj <- Surv(outcome_data$GTIME_KI, outcome_data$CR_STATUS == cause_num)
  
  # Get baseline survival
  base_fit <- survfit(surv_obj ~ 1, data = outcome_data)
  
  # Create time points
  time_points <- seq(50, 4000, by = 50)
  
  # Get baseline survival at these time points
  baseline_surv <- summary(base_fit, times = time_points, extend = TRUE)
  
  if(length(baseline_surv$time) < length(time_points)) {

    baseline_surv_extended <- approx(x = baseline_surv$time, 
                                     y = baseline_surv$surv,
                                     xout = time_points, 
                                     rule = 2)$y
  } else {
    baseline_surv_extended <- baseline_surv$surv
  }
  
  # Calculate cumulative incidence and confidence intervals for each group
  calculate_group_cif <- function(ethnicity_name, group_name) {
    
    # Create design matrix for this ethnicity with modal covariates
    if(ethnicity_name == "Asian") {
      eth_vec <- c(1, 0, 0)  # Asian contrast
    } else if(ethnicity_name == "Black") {
      eth_vec <- c(0, 1, 0)  # Black contrast  
    } else if(ethnicity_name == "Hispanic") {
      eth_vec <- c(0, 0, 1)  # Hispanic contrast
    } else { # White
      eth_vec <- c(-1, -1, -1)  # White contrast (sum-to-zero)
    }
    
    # Build design vector
    design_vec <- rep(0, length(coefs))
    names(design_vec) <- names(coefs)
    
    # Set ethnicity contrasts
    design_vec["ETHCATAsian"] <- eth_vec[1]
    design_vec["ETHCATBlack"] <- eth_vec[2] 
    design_vec["ETHCATHispanic"] <- eth_vec[3]

    if(modal_trtrej != "N" && paste0("TRTREJ1Y_KI", modal_trtrej) %in% names(design_vec)) {
      design_vec[paste0("TRTREJ1Y_KI", modal_trtrej)] <- 1
    }
    if(modal_work != "N" && paste0("WORK_INCOME_TRR", modal_work) %in% names(design_vec)) {
      design_vec[paste0("WORK_INCOME_TRR", modal_work)] <- 1
    }
    
    # Calculate linear predictor and its variance
    lp <- sum(design_vec * coefs)
    lp_var <- as.numeric(t(design_vec) %*% vcov(cox_model) %*% design_vec)
    lp_se <- sqrt(lp_var)
    
    # Calculate HR and confidence interval for HR
    hr <- exp(lp)
    hr_ci_lower <- exp(lp - 1.96 * lp_se)
    hr_ci_upper <- exp(lp + 1.96 * lp_se)
    
    # Calculate CIF and confidence intervals
    cif <- 1 - (baseline_surv_extended ^ hr)
    ci_lower <- 1 - (baseline_surv_extended ^ hr_ci_upper)  # Higher HR gives lower survival, higher CIF
    ci_upper <- 1 - (baseline_surv_extended ^ hr_ci_lower)  # Lower HR gives higher survival, lower CIF
    
    data.frame(
      time = time_points,
      cif = cif,
      ci_lower = pmax(0, ci_lower), 
      ci_upper = pmin(1, ci_upper), 
      group = group_name
    )
  }
  
  asian_data <- calculate_group_cif("Asian", "Asian")
  black_data <- calculate_group_cif("Black", "Black") 
  hispanic_data <- calculate_group_cif("Hispanic", "Hispanic")
  white_data <- calculate_group_cif("White", "White")
  

  combined_data <- rbind(asian_data, black_data, hispanic_data, white_data)

  result <- list(
    adj = combined_data
  )
  
  return(result)
}

adjcif_rejection_corrected <- calculate_corrected_cif(cs_model, 1, outcome_data)
adjcif_competing_corrected <- calculate_corrected_cif(cs_model, 2, outcome_data)

color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4",
                   "Hispanic" = "#6baed6", "White" = "#b3cdf8")

line_types <- c(
  "Asian" = "solid",         # Solid line (1)
  "Black" = "dashed",        # Dashed line (2)
  "Hispanic" = "dotdash",    # Dot-dash line (4)
  "White" = "11"             # Dotted line with shorter gaps
)

line_widths <- c(
  "Asian" = 1.2,
  "Black" = 1.2,
  "Hispanic" = 1.2,
  "White" = 1.5  
)


ethnicity_order <- c("Asian", "Black", "Hispanic", "White")

plot_cif <- function(adjcif, event_label, conf_int = TRUE) {
  # Extract the data from the adjcif object
  plot_data <- adjcif$adj

  plot_data$group <- factor(plot_data$group, levels = ethnicity_order)

  p <- ggplot(plot_data, aes(x = time, y = cif, color = group, linetype = group)) +
    geom_line(aes(linewidth = group)) +
    ylim(0, 0.4)+
    {if(conf_int) geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = group), 
                              alpha = 0.15, color = NA)} +
    scale_color_manual(values = color_palette, breaks = ethnicity_order) +
    scale_linetype_manual(values = line_types, breaks = ethnicity_order) +
    scale_fill_manual(values = color_palette, breaks = ethnicity_order) +
    scale_linewidth_manual(values = line_widths, breaks = ethnicity_order) +
    labs(title = event_label, 
         x = "Time (Days)", 
         y = "Cumulative Incidence") +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      legend.position = "bottom",
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13)
    ) +
    guides(
      color = guide_legend(title = "Ethnicity", override.aes = list(linewidth = line_widths)),
      linetype = guide_legend(title = "Ethnicity"),
      fill = guide_legend(title = "Ethnicity"),
      linewidth = "none"
    )
  
  return(p)
}

# Create the plots
p1 <- plot_cif(adjcif_rejection_corrected, "Graft Rejection", conf_int = TRUE)
p2 <- plot_cif(adjcif_competing_corrected, "Competing Risks", conf_int = TRUE)

# Combine plots
p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")


cat("Hazard Ratios for Rejection (Cause 1):\n")
cox_model_1 <- cs_model[["models"]][["Cause 1"]]
coefs_1 <- coef(cox_model_1)
hrs_1 <- exp(c(
  Asian = coefs_1["ETHCATAsian"],
  Black = coefs_1["ETHCATBlack"],
  Hispanic = coefs_1["ETHCATHispanic"],
  White = -(coefs_1["ETHCATAsian"] + coefs_1["ETHCATBlack"] + coefs_1["ETHCATHispanic"])
))
print(round(hrs_1, 3))



#########################################################################
# Likelihood Ratio Test for Variable Selection (Supplementary Material) #
#########################################################################
# Function to perform likelihood ratio test for a specific variable
test_variable_lrt <- function(data, full_formula, variable_to_test) {
  # Create reduced formula by removing the variable
  formula_text <- as.character(full_formula)
  lhs <- formula_text[2]  
  rhs <- formula_text[3]  

  terms <- unlist(strsplit(rhs, "\\+"))
  terms <- trimws(terms)
  new_terms <- terms[!grepl(variable_to_test, terms)]
  new_rhs <- paste(new_terms, collapse = " + ")
  
  # Create reduced formula
  reduced_formula <- as.formula(paste(lhs, "~", new_rhs))
  
  # Fit models
  cat("Fitting full model...\n")
  full_model <- CSC(full_formula, data = data)
  
  cat("Fitting reduced model without", variable_to_test, "...\n")
  reduced_model <- CSC(reduced_formula, data = data)
  
  results <- list()
  
  # Test for each cause
  for (cause in c(1, 2)) {
    cause_name <- paste0("Cause ", cause)
    cat("Testing", cause_name, "...\n")
    
    # Extract models and log-likelihoods
    full_cause <- full_model$models[[cause_name]]
    reduced_cause <- reduced_model$models[[cause_name]]
    ll_full <- logLik(full_cause)
    ll_reduced <- logLik(reduced_cause)
    
    # Calculate LRT
    lrt_stat <- 2 * (as.numeric(ll_full) - as.numeric(ll_reduced))
    df <- attr(ll_full, "df") - attr(ll_reduced, "df")
    p_value <- pchisq(lrt_stat, df = df, lower.tail = FALSE)
    
    results[[cause_name]] <- list(
      LRT = lrt_stat,
      df = df,
      p_value = p_value
    )
  }
  
  return(list(
    variable = variable_to_test,
    results = results
  ))
}

# Test variables
variables_to_test <- c("ETHCAT", "REGION", "EDUCATION", "CITIZENSHIP", "TRTREJ1Y_KI", "WORK_INCOME_TRR")
full_formula <- Hist(GTIME_KI, CR_STATUS) ~ ETHCAT + REGION + EDUCATION + CITIZENSHIP + TRTREJ1Y_KI + WORK_INCOME_TRR

lrt_results <- list()
for (var in variables_to_test) {
  cat("\n==== Testing", var, "====\n")
  result <- tryCatch({
    test_variable_lrt(outcome_data, full_formula, var)
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    lrt_results[[var]] <- result
  }
}

cat("\n========= LIKELIHOOD RATIO TEST RESULTS =========\n")
for (var in names(lrt_results)) {
  cat("\nVariable:", var, "\n")
  
  for (cause in c("Cause 1", "Cause 2")) {
    cat("  ", cause, ":\n")
    res <- lrt_results[[var]]$results[[cause]]
    cat("    LRT:", round(res$LRT, 3), 
        "\n    df:", res$df, 
        "\n    p-value:", format.pval(res$p_value, digits = 4), "\n")
  }
}

# Create summary table
summary_table <- data.frame(
  Variable = character(),
  Cause1_LRT = numeric(),
  Cause1_df = numeric(),
  Cause1_pvalue = numeric(),
  Cause2_LRT = numeric(),
  Cause2_df = numeric(),
  Cause2_pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (var in names(lrt_results)) {
  res <- lrt_results[[var]]
  row <- data.frame(
    Variable = var,
    Cause1_LRT = res$results$`Cause 1`$LRT,
    Cause1_df = res$results$`Cause 1`$df,
    Cause1_pvalue = res$results$`Cause 1`$p_value,
    Cause2_LRT = res$results$`Cause 2`$LRT,
    Cause2_df = res$results$`Cause 2`$df,
    Cause2_pvalue = res$results$`Cause 2`$p_value,
    stringsAsFactors = FALSE
  )
  summary_table <- rbind(summary_table, row)
}

# Sort by Cause 1 p-value (rejection is primary interest)
summary_table <- summary_table[order(summary_table$Cause1_pvalue), ]


summary_table %>% 
  mutate(
    Sig_Cause1 = Cause1_pvalue < 0.05,
    Sig_Cause2 = Cause2_pvalue < 0.05,
    Any_Significant = Sig_Cause1 | Sig_Cause2,
    
    Marginal_Cause1 = Cause1_pvalue >= 0.05 & Cause1_pvalue < 0.10,
    Marginal_Cause2 = Cause2_pvalue >= 0.05 & Cause2_pvalue < 0.10,
    Any_Marginal = Marginal_Cause1 | Marginal_Cause2,
    
    # Decision recommendation
    Decision = case_when(
      Variable == "ETHCAT" ~ "Keep (Primary Variable)",
      Any_Significant ~ "Keep (Significant)",
      Any_Marginal ~ "Consider (Marginal)",
      TRUE ~ "Remove (Not Significant)"
    ),
    
    Cause1_pvalue = format.pval(Cause1_pvalue, digits = 3),
    Cause2_pvalue = format.pval(Cause2_pvalue, digits = 3)
  ) %>% select(
    Variable, Cause1_LRT, Cause1_pvalue, Sig_Cause1,
    Cause2_LRT, Cause2_pvalue, Sig_Cause2, Decision
  )

