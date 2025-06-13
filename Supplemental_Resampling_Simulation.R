library(dplyr)
library(purrr)
library(MASS)
library(ggplot2)
library(Benchmarking)

dmu_data <- read.csv("dmu_data.csv")
Final_Data <- read.csv("Final_Data.csv")

treatment_effect <- data.frame(
  Group = Final_Data$Group,
  Year = Final_Data$Year,
  WaitlistDuration = Final_Data$WaitlistDuration,
  LKDPI = Final_Data$LKDPI,
  GraftLifespan = Final_Data$GTIME_KI,
  EDUCATION = as.factor(Final_Data$EDUCATION),
  REGION = as.factor(Final_Data$REGION),
  CITIZENSHIP = as.factor(Final_Data$CITIZENSHIP),
  AGE = Final_Data$AGE,
  GENDER = Final_Data$GENDER
)

# ESRD Prevalence Counts by Year (https://usrds-adr.niddk.nih.gov/2024/reference-tables, Table B.1)
esrd_counts_by_year <- list(
  "2010" = c(Asian = 25048, Black = 183948, Hispanic = 94650, White = 276777),
  "2011" = c(Asian = 26676, Black = 190914, Hispanic = 100440, White = 282207),
  "2012" = c(Asian = 28345, Black = 197394, Hispanic = 105730, White = 289956),
  "2013" = c(Asian = 30148, Black = 204333, Hispanic = 111778, White = 298822),
  "2014" = c(Asian = 31956, Black = 211315, Hispanic = 117768, White = 308436),
  "2015" = c(Asian = 33800, Black = 217563, Hispanic = 124270, White = 318856),
  "2016" = c(Asian = 35825, Black = 222328, Hispanic = 130416, White = 328987),
  "2017" = c(Asian = 37810, Black = 226188, Hispanic = 136111, White = 338207),
  "2018" = c(Asian = 40088, Black = 230048, Hispanic = 142358, White = 348471),
  "2019" = c(Asian = 42642, Black = 235844, Hispanic = 149856, White = 356921)
)

# Function to convert counts to proportions
get_proportions <- function(counts) {
  total <- sum(counts)
  round(counts / total, 5)
}

# Compute covariance matrices for each group
groups <- c("Asian", "Black", "Hispanic", "White")
cov_matrices <- list()

for (group in groups) {
  group_data <- treatment_effect %>%
    filter(Group == group) %>%
    dplyr::select(WaitlistDuration, LKDPI, GraftLifespan)
  cov_matrix <- cov(group_data)
  cov_matrices[[group]] <- cov_matrix
}

# Data summary for simulation
data_summary <- treatment_effect %>% 
  group_by(Group) %>% 
  summarise(
    wait_mean = mean(WaitlistDuration),
    wait_median = median(WaitlistDuration),
    wait_sd = sd(WaitlistDuration),
    lkdpi_mean = mean(LKDPI),  # Changed from quality to LKDPI
    lkdpi_median = median(LKDPI),
    lkdpi_sd = sd(LKDPI),
    graft_mean = mean(GraftLifespan),  # Changed from outcome to graft
    graft_median = median(GraftLifespan),
    graft_sd = sd(GraftLifespan)
  )

# Reference Frontier Mapping (RFM) function
create_reference_set <- function(data, ref_prop = 0.1, gamma = 0.05) {
  shift_scores_positive <- function(data) {
    data$WaitlistDuration <- data$WaitlistDuration + abs(min(data$WaitlistDuration, na.rm = TRUE))
    data$LKDPI <- data$LKDPI + abs(min(data$LKDPI, na.rm = TRUE))
    data$GraftLifespan <- data$GraftLifespan + abs(min(data$GraftLifespan, na.rm = TRUE))
    return(data)
  }
  data <- shift_scores_positive(data)
  get_diverse_reference <- function(data) {
    bind_rows(
      data %>% filter(WaitlistDuration <= quantile(WaitlistDuration, gamma)),
      data %>% filter(LKDPI <= quantile(LKDPI, gamma)),
      data %>% filter(GraftLifespan >= quantile(GraftLifespan, 1 - gamma))
    ) %>% distinct()
  }
  total_n <- nrow(data)
  target_n_ref <- ceiling(ref_prop * total_n)
  diverse_reference <- get_diverse_reference(data)
  n_diverse <- nrow(diverse_reference)
  if (n_diverse > target_n_ref) {
    # Score observations based on frontier relevance
    frontier_scores <- diverse_reference %>%
      mutate(
        priority_score = 1 / (WaitlistDuration - min(WaitlistDuration) + 1),
        quality_score = 1 / (LKDPI - min(LKDPI) + 1),
        outcome_score = (GraftLifespan - min(GraftLifespan)) / (max(GraftLifespan) - min(GraftLifespan)),
        frontier_value = priority_score * quality_score * outcome_score
      )
    reference_data <- frontier_scores %>%
      arrange(desc(frontier_value)) %>%
      slice_head(n = target_n_ref) %>%
      dplyr::select(-priority_score, -quality_score, -outcome_score, -frontier_value)
  } else {
    n_remaining <- target_n_ref - n_diverse
    remaining_pool <- anti_join(data, diverse_reference)
    
    if(nrow(remaining_pool) > 0) {
      additional_reference <- remaining_pool %>% 
        slice_sample(n = min(n_remaining, nrow(remaining_pool)))
      reference_data <- bind_rows(diverse_reference, additional_reference)
    } else {
      reference_data <- diverse_reference
    }
  }
  
  evaluation_data <- anti_join(data, reference_data)
  
  return(list(reference = reference_data, evaluation = evaluation_data))
}

# Conditional DEA with RFM
compute_conditional_efficiency_rfm <- function(data, exog_vars = c("REGION", "EDUCATION", "CITIZENSHIP"), 
                                               bandwidth = NULL) {
  # Use optimal bandwidth if not specified
  if(is.null(bandwidth)) {
    r <- length(exog_vars)
    n <- nrow(data)
    bandwidth <- n^(-1/(r+4))  # Optimal bandwidth from Badin et al. (2012)
  }
  
  # Create reference and evaluation sets using RFM
  rfm_sets <- create_reference_set(data)
  reference_data <- rfm_sets$reference
  evaluation_data <- rfm_sets$evaluation
  if(nrow(reference_data) < 3 || nrow(evaluation_data) == 0) {
    return(data.frame(Group = character(0), efficiency = numeric(0)))
  }
  # Split exogenous variables
  split_exog_vars <- function(data, exog_vars) {
    continuous_vars <- exog_vars[sapply(data[, exog_vars], is.numeric)]
    categorical_vars <- exog_vars[!(exog_vars %in% continuous_vars)]
    list(continuous = continuous_vars, categorical = categorical_vars)
  }
  # Product kernel function
  product_kernel <- function(ref_data, eval_row, cont_vars, cat_vars, h) {
    # Continuous: 
    if (length(cont_vars) > 0) {
      Zc_ref <- as.matrix(ref_data[, cont_vars, drop = FALSE])
      Zc_eval <- as.numeric(eval_row[cont_vars])
      cont_dists <- sqrt(rowSums((Zc_ref - matrix(Zc_eval, nrow = nrow(Zc_ref), 
                                                  ncol = length(Zc_eval), byrow = TRUE))^2))
      Kc <- exp(- (cont_dists^2) / (2 * h^2))
    } else {
      Kc <- rep(1, nrow(ref_data))
    }
    # Categorical: 
    if (length(cat_vars) > 0) {
      Zd_ref <- ref_data[, cat_vars, drop = FALSE]
      Zd_eval <- eval_row[cat_vars]
      
      match_matrix <- mapply(function(col_ref, val_eval) {
        match_vec <- as.character(col_ref) == as.character(val_eval)
        smoothed_vec <- ifelse(match_vec, 1, exp(-1 / h^2))
        return(smoothed_vec)
      }, Zd_ref, Zd_eval, SIMPLIFY = FALSE)
      
      Kd <- Reduce(`*`, match_matrix)
    } else {
      Kd <- rep(1, nrow(ref_data))
    }
    
    return(Kc * Kd)
  }
  # Define variables
  split_vars <- split_exog_vars(reference_data, exog_vars)
  cont_vars <- split_vars$continuous
  cat_vars <- split_vars$categorical
  
  # Inputs and outputs
  X_ref <- as.matrix(reference_data[, c("WaitlistDuration", "LKDPI")])
  Y_ref <- as.matrix(reference_data[, "GraftLifespan"])
  
  X_eval <- as.matrix(evaluation_data[, c("WaitlistDuration", "LKDPI")])
  Y_eval <- as.matrix(evaluation_data[, "GraftLifespan"])
  
  # Conditional DEA computation
  eff_scores <- numeric(nrow(X_eval))
  
  for (i in 1:nrow(X_eval)) {
    tryCatch({
      eval_row <- evaluation_data[i, ]
      weights <- product_kernel(reference_data, eval_row, cont_vars, cat_vars, bandwidth)
      selected_idx <- which(weights > 1e-5)
      
      if (length(selected_idx) >= 2) {
        X_local <- X_ref[selected_idx, , drop = FALSE]
        Y_local <- Y_ref[selected_idx, , drop = FALSE]
        
        eff_scores[i] <- Benchmarking::dea(
          X = X_eval[i, , drop = FALSE],
          Y = Y_eval[i, , drop = FALSE],
          XREF = X_local,
          YREF = Y_local,
          RTS = "vrs",
          ORIENTATION = "graph"
        )$eff
      } else {
        eff_scores[i] <- NA
      }
    }, error = function(e) {
      eff_scores[i] <- NA
    })
  }
  
  # Add efficiency scores to evaluation data
  evaluation_data$efficiency <- eff_scores
  evaluation_data <- evaluation_data %>%
    group_by(Year) %>%
    mutate(relative_efficiency = efficiency - mean(efficiency, na.rm = TRUE)) %>%
    ungroup()
  
  return(evaluation_data)
}

# Year-specific stratified sampling function
stratified_sample_yearly <- function(data, group_var, esrd_counts_by_year) {
  years <- intersect(names(esrd_counts_by_year), unique(data$Year) |> as.character())
  
  balanced_data_list <- lapply(years, function(year) {
    year_data <- data %>% filter(Year == as.integer(year))
    if (nrow(year_data) == 0) return(NULL)
    
    desired_props <- get_proportions(esrd_counts_by_year[[year]])
    group_sizes <- table(year_data[[group_var]])
    
    # Calculate maximum possible balanced sample size
    max_possible_samples <- min(group_sizes / desired_props[names(group_sizes)])
    total_samples <- floor(max_possible_samples)
    
    if(total_samples < 10) return(NULL) 
    
    # Sample from each group
    desired_counts <- round(desired_props[names(group_sizes)] * total_samples)
    sampled_data <- map2_dfr(split(year_data, year_data[[group_var]]), 
                             names(desired_counts), function(group_data, group_name) {
                               slice_sample(group_data, n = min(nrow(group_data), desired_counts[group_name]))
                             })
    
    return(sampled_data)
  })
  
  return(bind_rows(balanced_data_list))
}

# Generate synthetic data function
generate_synthetic_data <- function(n = 1000, summary_stats = data_summary, 
                                    cov_matrices, use_original_props = TRUE) {
  if(use_original_props) {
    # Use original imbalanced proportions
    groups <- sample(c("Asian", "Black", "Hispanic", "White"), n, replace = TRUE,
                     prob = as.vector(table(dmu_data$Group)/nrow(dmu_data)))
  } else {
    # Use average ESRD proportions across years
    avg_props <- get_proportions(rowMeans(sapply(esrd_counts_by_year, function(x) x)))
    groups <- sample(c("Asian", "Black", "Hispanic", "White"), n, replace = TRUE,
                     prob = avg_props[c("Asian", "Black", "Hispanic", "White")])
  }
  
  # Add year information (randomly assign years 2010-2019)
  years <- sample(2010:2019, n, replace = TRUE)
  
  data <- lapply(unique(groups), function(group) {
    group_indices <- which(groups == group)
    group_n <- length(group_indices)
    stats <- summary_stats %>% filter(Group == group)
    means <- c(stats$wait_mean, stats$lkdpi_mean, stats$graft_mean)  
    cov_matrix <- cov_matrices[[group]]
    
    synthetic_data <- MASS::mvrnorm(n = group_n, mu = means, Sigma = cov_matrix)
    synthetic_data <- pmax(synthetic_data, 0)  
    
    tibble(
      Group = group,
      Year = years[group_indices],
      WaitlistDuration = synthetic_data[, 1],
      LKDPI = synthetic_data[, 2],  
      GraftLifespan = synthetic_data[, 3], 
      # Add dummy exogenous variables
      REGION = sample(1:11, group_n, replace = TRUE),
      EDUCATION = sample(1:6, group_n, replace = TRUE),
      CITIZENSHIP = sample(1:3, group_n, replace = TRUE)
    )
  })
  
  bind_rows(data)
}

run_simulation_rfm <- function(n_iterations = 100) {
  results <- tibble()
  
  for (i in 1:n_iterations) {
    cat(sprintf("Running iteration %d/%d...\n", i, n_iterations))
    
    # Generate synthetic data with original proportions
    synthetic_data <- generate_synthetic_data(n = 1000, summary_stats = data_summary, 
                                              cov_matrices = cov_matrices, 
                                              use_original_props = TRUE)
    
    # Compute efficiency for original (imbalanced) data using RFM
    original_efficiency <- compute_conditional_efficiency_rfm(synthetic_data)
    
    # Perform year-specific stratified resampling to match ESRD proportions
    resampled_data <- stratified_sample_yearly(synthetic_data, "Group", esrd_counts_by_year)
    
    if(nrow(resampled_data) > 50) { 
      # Compute efficiency for resampled data using RFM
      resampled_efficiency <- compute_conditional_efficiency_rfm(resampled_data)

      results <- bind_rows(results,
                           original_efficiency %>%
                             mutate(Method = "Imbalanced", Iteration = i) %>%
                             dplyr::select(Group, relative_efficiency, Method, Iteration),
                           resampled_efficiency %>%
                             mutate(Method = "Resampled (ESRD)", Iteration = i) %>%
                             dplyr::select(Group, relative_efficiency, Method, Iteration))
    }
  }
  
  return(results)
}

# Run the simulation
set.seed(123)
cat("Starting RFM-based simulation with LKDPI and yearly ESRD resampling...\n")
sim_results <- run_simulation_rfm(n_iterations = 100) 

# Analyze results
summary_stats <- sim_results %>%
  group_by(Group, Method) %>%
  summarise(
    Mean = mean(relative_efficiency, na.rm = TRUE),
    SD = sd(relative_efficiency, na.rm = TRUE),
    CI_lower = quantile(relative_efficiency, 0.025, na.rm = TRUE),
    CI_upper = quantile(relative_efficiency, 0.975, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# Sensitivity analysis with different proportion scenarios
sensitivity_proportions <- list(
  "Alternative 1" = c(Asian = 0.1, Black = 0.2, Hispanic = 0.3, White = 0.4),
  "Alternative 2" = c(Asian = 0.2, Black = 0.4, Hispanic = 0.3, White = 0.1)
)


# Sensitivity analysis with different proportion scenarios
run_sensitivity_analysis_rfm <- function(n_iterations = 100, alternative_proportions) {
  results <- tibble()
  
  # Average ESRD proportions across years as baseline
  avg_esrd_props <- get_proportions(rowMeans(sapply(esrd_counts_by_year, function(x) x)))
  baseline_props <- avg_esrd_props[c("Asian", "Black", "Hispanic", "White")]
  
  # Test each proportion scenario
  all_scenarios <- c(list("ESRD Average" = baseline_props), alternative_proportions)
  
  for (scenario_name in names(all_scenarios)) {
    cat(sprintf("Testing scenario: %s\n", scenario_name))
    scenario_props <- all_scenarios[[scenario_name]]
    
    for (i in 1:n_iterations) {
      if(i %% 10 == 0) cat(sprintf("  Iteration %d/%d\n", i, n_iterations))
      
      # Generate synthetic data using baseline ESRD proportions
      synthetic_data <- generate_synthetic_data(n = 1000, summary_stats = data_summary, 
                                                cov_matrices = cov_matrices, 
                                                use_original_props = FALSE)  # Use ESRD props
      
      # Compute efficiency for baseline ESRD proportions
      if(scenario_name == "ESRD Average") {
        scenario_efficiency <- compute_conditional_efficiency_rfm(synthetic_data)
        method_label <- "ESRD Average"
      } else {
        # Resample according to alternative proportions
        group_sizes <- table(synthetic_data$Group)
        max_possible_samples <- min(group_sizes / scenario_props[names(group_sizes)])
        total_samples <- floor(max_possible_samples)
        
        if(total_samples > 50) {
          # Perform stratified sampling for alternative scenario
          desired_counts <- round(scenario_props[names(group_sizes)] * total_samples)
          resampled_data <- map2_dfr(split(synthetic_data, synthetic_data$Group), 
                                     names(desired_counts), function(group_data, group_name) {
                                       if(group_name %in% names(desired_counts)) {
                                         slice_sample(group_data, n = min(nrow(group_data), desired_counts[group_name]))
                                       } else {
                                         return(tibble())
                                       }
                                     })
          
          scenario_efficiency <- compute_conditional_efficiency_rfm(resampled_data)
          method_label <- scenario_name
        } else {
          next  
        }
      }

      if(nrow(scenario_efficiency) > 0) {
        results <- bind_rows(results,
                             scenario_efficiency %>%
                               mutate(Method = method_label, 
                                      Scenario = scenario_name,
                                      Iteration = i) %>%
                               dplyr::select(Group, relative_efficiency, Method, Scenario, Iteration))
      }
    }
  }
  
  return(results)
}

# Define alternative proportion scenarios
sensitivity_proportions <- list(
  "Original Imbalanced" = table(dmu_data$Group)/nrow(dmu_data),
  "Alternative 1" = c(Asian = 0.1, Black = 0.2, Hispanic = 0.3, White = 0.4),
  "Alternative 2" = c(Asian = 0.2, Black = 0.4, Hispanic = 0.3, White = 0.1)
)

# Run sensitivity analysis
set.seed(123)
cat("Starting sensitivity analysis...\n")
sensitivity_results <- run_sensitivity_analysis_rfm(n_iterations = 100, 
                                                    alternative_proportions = sensitivity_proportions)

 # Calculate mean efficiency for visualization
mean_efficiency_sensitivity <- sensitivity_results %>%
  group_by(Group, Scenario) %>%
  summarise(mean_efficiency = mean(relative_efficiency, na.rm = TRUE), .groups = 'drop')

# Create proportion labels
scenario_order <- c("ESRD Average", "Original Imbalanced", "Alternative 1", "Alternative 2")
scenario_labels <- c(
  "ESRD Average" = "ESRD Average",
  "Original Imbalanced" = "Original Imbalanced Data", 
  "Alternative 1" = "Alternative 1",
  "Alternative 2" = "Alternative 2"
)

# Sensitivity analysis plot
color_palette <- c("Asian" = "#08306b", "Black" = "#1f78b4", 
                   "Hispanic" = "#6baed6", "White" = "#b3cde3")
ggplot(mean_efficiency_sensitivity, 
                           aes(x = factor(Scenario, levels = scenario_order), 
                               y = mean_efficiency, color = Group, group = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 3.5) +
  scale_x_discrete(labels = scenario_labels) +
  scale_color_manual(values = color_palette) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "Mean Relative Efficiency Score",
    x = "Proportion Scenarios",
    #title = "Sensitivity Analysis: RFM-Based Efficiency Across Proportion Scenarios",
    #subtitle = "LKDPI-based analysis with conditional DEA and Reference Frontier Mapping",
    color = "Ethnic Group"
  )




# Summary statistics for sensitivity analysis
sensitivity_summary <- sensitivity_results %>%
  group_by(Group, Scenario) %>%
  summarise(
    Mean = mean(relative_efficiency, na.rm = TRUE),
    SD = sd(relative_efficiency, na.rm = TRUE),
    CI_lower = quantile(relative_efficiency, 0.025, na.rm = TRUE),
    CI_upper = quantile(relative_efficiency, 0.975, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Scenario, Group)

print("Sensitivity Analysis Summary:")
print(sensitivity_summary)

sensitivity_table <- sensitivity_summary %>%
  dplyr::select(Scenario, Group, Mean) %>%
  pivot_wider(names_from = Group, values_from = Mean) %>%
  arrange(factor(Scenario, levels = scenario_order))

print("Mean Relative Efficiency Scores Across Different Proportion Scenarios:")
print(sensitivity_table)

# Bias calculation relative to ESRD Average
bias_calculation <- sensitivity_summary %>%
  dplyr::select(Group, Scenario, Mean) %>%
  pivot_wider(names_from = Scenario, values_from = Mean) %>%
  mutate(
    `Original Imbalanced Bias` = `Original Imbalanced` - `ESRD Average`,
    `Alternative 1 Bias` = `Alternative 1` - `ESRD Average`,
    `Alternative 2 Bias` = `Alternative 2` - `ESRD Average`
  ) %>%
  dplyr::select(Group, contains("Bias"))

print("Bias in Efficiency Scores Relative to ESRD Average:")
print(bias_calculation)

# Calculate weighted bias for each scenario
weighted_bias_calculation <- bias_calculation %>%
  pivot_longer(cols = contains("Bias"), names_to = "Scenario", values_to = "Bias") %>%
  mutate(
    Scenario = str_remove(Scenario, " Bias"),
    Group_Proportion = case_when(
      Scenario == "Original Imbalanced" & Group == "Asian" ~ 0.078,
      Scenario == "Original Imbalanced" & Group == "Black" ~ 0.289,
      Scenario == "Original Imbalanced" & Group == "Hispanic" ~ 0.130,
      Scenario == "Original Imbalanced" & Group == "White" ~ 0.502,
      
      Scenario == "Alternative 1" & Group == "Asian" ~ 0.1,
      Scenario == "Alternative 1" & Group == "Black" ~ 0.2,
      Scenario == "Alternative 1" & Group == "Hispanic" ~ 0.3,
      Scenario == "Alternative 1" & Group == "White" ~ 0.4,
      
      Scenario == "Alternative 2" & Group == "Asian" ~ 0.2,
      Scenario == "Alternative 2" & Group == "Black" ~ 0.4,
      Scenario == "Alternative 2" & Group == "Hispanic" ~ 0.3,
      Scenario == "Alternative 2" & Group == "White" ~ 0.1,
      
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(Scenario) %>%
  summarise(Weighted_Bias = sum(Bias * Group_Proportion, na.rm = TRUE), .groups = 'drop')

print("Weighted Bias by Scenario:")
print(weighted_bias_calculation)

bias_summary_table <- bias_calculation %>%
  pivot_longer(cols = contains("Bias"), names_to = "Scenario_Type", values_to = "Bias") %>%
  mutate(Scenario_Type = str_remove(Scenario_Type, " Bias")) %>%
  pivot_wider(names_from = Scenario_Type, values_from = Bias) %>%
  bind_rows(
    tibble(
      Group = "Weighted Bias",
      `Original Imbalanced` = weighted_bias_calculation$Weighted_Bias[1],
      `Alternative 1` = weighted_bias_calculation$Weighted_Bias[2],
      `Alternative 2` = weighted_bias_calculation$Weighted_Bias[3]
    )
  )

print("Comprehensive Bias Analysis Table:")
print(bias_summary_table)
