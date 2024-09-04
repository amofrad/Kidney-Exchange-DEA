library(maotai)
library(kernlab)
library(dplyr)


balanced_data <- read.csv("balanced_data.csv")

treatment_effect <- data.frame(
  Group = balanced_data$Group,
  Year = balanced_data$Year,
  WaitlistDuration = balanced_data$Waitlist_Score,
  QualityScore = balanced_data$Quality_Score,
  OutcomeScore = balanced_data$Outcome_Score
)
treatment_effect$WaitlistDuration <- treatment_effect$WaitlistDuration + abs(min(treatment_effect$WaitlistDuration))
treatment_effect$QualityScore <- treatment_effect$QualityScore + abs(min(treatment_effect$QualityScore))
treatment_effect$OutcomeScore <- treatment_effect$OutcomeScore + abs(min(treatment_effect$OutcomeScore))



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