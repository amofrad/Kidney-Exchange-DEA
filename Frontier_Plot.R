library(Benchmarking)

Final_Data <- read.csv("Final_Data.csv")

treatment_effect <- data.frame(
  Group = Final_Data$Group,
  Year = Final_Data$Year,
  PriorityScore = Final_Data$WaitlistDuration,
  QualityScore = Final_Data$LKDPI,
  OutcomeScore = Final_Data$GTIME_KI,
  EDUCATION = as.factor(Final_Data$EDUCATION),
  REGION = as.factor(Final_Data$REGION),
  CITIZENSHIP = as.factor(Final_Data$CITIZENSHIP),
  AGE = Final_Data$AGE,
  GENDER = Final_Data$GENDER
)

treatment_effect_2010_2016 <- treatment_effect[treatment_effect$Year <=2016,]

################################################################################
###### Plotting Conditional DEA Production Possibility Frontier (Fig. 1) #######
################################################################################

# Function to get local reference set for a given evaluation point
get_local_reference <- function(eval_row, reference_data, cont_vars, cat_vars, bandwidth = nrow(reference_data)^(-1/7)) {
  product_kernel <- function(ref_data, eval_row, cont_vars, cat_vars, h) {
    # Continuous: 
    if (length(cont_vars) > 0) {
      Zc_ref <- as.matrix(ref_data[, cont_vars, drop = FALSE])
      Zc_eval <- as.numeric(eval_row[cont_vars])
      cont_dists <- sqrt(rowSums((Zc_ref - matrix(Zc_eval, nrow = nrow(Zc_ref), ncol = length(Zc_eval), byrow = TRUE))^2))
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
  
  weights <- product_kernel(reference_data, eval_row, cont_vars, cat_vars, bandwidth)
  selected_idx <- which(weights > 1e-5)
  
  return(list(
    indices = selected_idx,
    weights = weights[selected_idx],
    reference_subset = reference_data[selected_idx, ]
  ))
}

# Function to create conditional frontier points for a local reference set
create_conditional_frontier <- function(local_ref_x, local_ref_y, xlim, ylim, steps = 100) {
  if(length(local_ref_x) < 2) return(NULL)
  
  # Perform DEA on the local reference set to get the frontier
  dea_result <- dea(X = matrix(local_ref_x, ncol = 1),
                    Y = matrix(local_ref_y, ncol = 1),
                    RTS = "vrs",
                    ORIENTATION = "graph")
  
  # Get efficient points (those on the frontier)
  efficient_points <- which(dea_result$eff >= 0.999)  # Account for numerical precision
  
  if(length(efficient_points) < 2) return(NULL)
  
  # Sort efficient points by x-coordinate for frontier plotting
  frontier_x <- local_ref_x[efficient_points]
  frontier_y <- local_ref_y[efficient_points]
  sorted_order <- order(frontier_x)
  frontier_x <- frontier_x[sorted_order]
  frontier_y <- frontier_y[sorted_order]
  extended_x <- frontier_x
  extended_y <- frontier_y
  leftmost_x <- frontier_x[1]
  rightmost_x <- frontier_x[length(frontier_x)]
  rightmost_y <- frontier_y[length(frontier_y)]
  
  if(leftmost_x > xlim[1] + 0.05) {
    extension_length <- min(0.1, leftmost_x - xlim[1])
    extended_x <- c(leftmost_x - extension_length, extended_x)
    extended_y <- c(frontier_y[1], extended_y)
  }
  
  if(rightmost_x < xlim[2] - 0.05) {
    extension_length <- min(0.1, xlim[2] - rightmost_x)
    extended_x <- c(extended_x, rightmost_x + extension_length)
    extended_y <- c(extended_y, rightmost_y)
  }
  
  return(list(
    x = extended_x,
    y = extended_y
  ))
}


calculate_conditional_hyperbolic_path <- function(X, Y, local_ref_x, local_ref_y, steps = 100) {
  # Get efficiency
  dea_result <- dea(X = matrix(X, ncol = 1),
                    Y = matrix(Y, ncol = 1),
                    XREF = matrix(local_ref_x, ncol = 1),
                    YREF = matrix(local_ref_y, ncol = 1),
                    RTS = "vrs",
                    ORIENTATION = "graph")
  
  eff <- dea_result$eff
  
  # Calculate hyperbolic path
  t <- seq(1, eff, length.out = steps)
  path_X <- X * t
  path_Y <- Y / t
  
  return(list(X = path_X, Y = path_Y, eff = eff))
}

# Prepare the data - Convert group numbers to labels
DEA_Result_2010_2016 <- DEA_Result_2010_2016 %>%
  mutate(Group = case_when(
    Group == "1" ~ "Asian",
    Group == "2" ~ "Black", 
    Group == "3" ~ "Hispanic",
    Group == "4" ~ "White",
    TRUE ~ as.character(Group)
  ))
DEA_Result_2010_2016$Group <- as.factor(DEA_Result_2010_2016$Group)

color_palette <- c(
  "Asian" = "#08306b",
  "Black" = "#1f78b4", 
  "Hispanic" = "#6baed6",
  "White" = "#b3cde3"
)


min_max_scale <- function(data) {
  scale_column <- function(x) {
    rng <- range(x, na.rm = TRUE)
    (x - rng[1]) / (rng[2] - rng[1])
  }
  
  data$PriorityScore <- scale_column(data$PriorityScore)
  data$QualityScore  <- scale_column(data$QualityScore)
  data$OutcomeScore  <- scale_column(data$OutcomeScore)
  
  return(data)
}

shift_scores_positive <- function(data) {
  data$PriorityScore <- data$PriorityScore + abs(min(data$PriorityScore, na.rm = TRUE))
  data$QualityScore  <- data$QualityScore  + abs(min(data$QualityScore, na.rm = TRUE))
  data$OutcomeScore  <- data$OutcomeScore  + abs(min(data$OutcomeScore, na.rm = TRUE))
  return(data)
}


gamma <- 0.05
ref_prop <- 0.10
treatment_effect_scaled <- min_max_scale(treatment_effect_2010_2016)

# Create reference set
get_diverse_reference <- function(data) {
  bind_rows(
    data %>% filter(PriorityScore <= quantile(PriorityScore, gamma)),
    data %>% filter(QualityScore <= quantile(QualityScore, gamma)),
    data %>% filter(OutcomeScore >= quantile(OutcomeScore, 1 - gamma))
  ) %>% distinct()
}

target_n_ref <- ceiling(ref_prop * nrow(treatment_effect_scaled))

diverse_reference <- get_diverse_reference(treatment_effect_scaled)
n_diverse <- nrow(diverse_reference)

if (n_diverse > target_n_ref) {
  # Score each observation based on how extreme it is
  frontier_scores <- diverse_reference %>%
    mutate(
      # Higher scores = more valuable for frontier construction
      priority_score = 1 / (PriorityScore - min(PriorityScore) + 1),  # Favor low priority scores
      quality_score = 1 / (QualityScore - min(QualityScore) + 1),     # Favor low quality scores  
      outcome_score = (OutcomeScore - min(OutcomeScore)) / (max(OutcomeScore) - min(OutcomeScore)), # Favor high outcome scores
      # Combined frontier value
      frontier_value = priority_score * quality_score * outcome_score
    )
  
  reference_data <- frontier_scores %>%
    arrange(desc(frontier_value)) %>%  
    slice_head(n = target_n_ref) %>%
    dplyr::select(-priority_score, -quality_score, -outcome_score, -frontier_value)
  
  suppressMessages(evaluation_data <- anti_join(treatment_effect_scaled, reference_data))
  
} else {
  n_remaining <- target_n_ref - n_diverse
  remaining_pool <- anti_join(treatment_effect, diverse_reference)
  
  additional_reference <- remaining_pool %>% slice_sample(n = n_remaining)
  reference_data <- bind_rows(diverse_reference, additional_reference)
  evaluation_data <- anti_join(treatment_effect_scaled, reference_data)
}

# Merge efficiency scores from DEA_Result_2010_2016 to the scaled evaluation data
evaluation_data$efficiency <- DEA_Result_2010_2016$efficiency[match(
  paste(evaluation_data$Group, evaluation_data$Year, evaluation_data$AGE),
  paste(DEA_Result_2010_2016$Group, DEA_Result_2010_2016$Year, DEA_Result_2010_2016$AGE)
)]

# Define variable splits
cont_vars <- c()  
cat_vars <- c("REGION", "EDUCATION", "CITIZENSHIP")

eval_ref_ids <- sapply(seq_len(nrow(evaluation_data)), function(i) {
  lr <- get_local_reference(evaluation_data[i, ], reference_data, cont_vars, cat_vars)
  if (length(lr$indices) >= 10) {
    paste(sort(lr$indices), collapse = "_")
  } else {
    NA
  }
})

ref_counts <- table(eval_ref_ids, useNA = "no")
common_id  <- names(ref_counts)[which.max(ref_counts)]
all_inds <- which(eval_ref_ids == common_id)
set.seed(2)
if(length(all_inds) >= 100) {
  chosen_inds <- sample(all_inds, 100)
} else {
  chosen_inds <- all_inds
}

treatment_effect_subset <- evaluation_data[chosen_inds, ]

n_hyper <- ceiling(0.10 * nrow(treatment_effect_subset))

ord        <- order(treatment_effect_subset$efficiency, decreasing = TRUE)
hyperbolic_indices <- ord[1:n_hyper]


par(mfrow = c(1, 2), mar = c(5, 5, 1, 2) + 0.1, oma = c(0, 0, 4, 0))
xlim_common <- c(0, 1)
ylim_common <- c(0, 1)

# Plot 1: Priority Score vs Outcome Score
plot(treatment_effect_subset$PriorityScore, 
     treatment_effect_subset$OutcomeScore,
     col = color_palette[treatment_effect_subset$Group],
     pch = 19, cex = 1.,
     xlim = xlim_common, ylim = ylim_common,
     xlab = "Waitlist Duration", ylab = "Graft Lifespan",
     main = "",
     cex.lab = 1.5, cex.axis = 1.2,
     font.lab = 2, font.axis = 1,
     bty        = "l" )

# Add conditional frontiers and hyperbolic path lines for selected points
frontiers_plotted <- list()

for(i in hyperbolic_indices) {
  eval_row <- treatment_effect_subset[i, ]
  
  # Get local reference set for this evaluation point
  local_ref <- get_local_reference(eval_row, reference_data, cont_vars, cat_vars)
  
  print(length(local_ref$indices))
  
  if(length(local_ref$indices) >= 3) {  # Need at least 3 points for a meaningful frontier
    local_ref_data <- local_ref$reference_subset
    
    # Create a unique identifier for this local reference set
    ref_id <- paste(sort(local_ref$indices), collapse = "_")
    
    # Only plot the frontier once per unique local reference set
    if(!(ref_id %in% names(frontiers_plotted))) {
      frontier <- create_conditional_frontier(local_ref_data$PriorityScore, 
                                              local_ref_data$OutcomeScore,
                                              xlim_common, ylim_common)
      
      if(!is.null(frontier)) {
        # Plot the conditional frontier
        lines(frontier$x, frontier$y, col = "blue", lwd = 3, lty = 1)
        frontiers_plotted[[ref_id]] <- TRUE
      }
    }
    
    # Calculate conditional hyperbolic path
    path <- calculate_conditional_hyperbolic_path(
      treatment_effect_subset$PriorityScore[i],
      treatment_effect_subset$OutcomeScore[i],
      local_ref_data$PriorityScore,
      local_ref_data$OutcomeScore
    )
    
    # Draw hyperbolic path
    lines(path$X, path$Y, col = "darkgray", lty = 2, lwd = 2)
    points(path$X[length(path$X)], path$Y[length(path$Y)], pch = 8, col = "red", cex = 2)
  }
}

# Plot 2: Quality Score vs Outcome Score
plot(treatment_effect_subset$QualityScore, 
     treatment_effect_subset$OutcomeScore,
     col = color_palette[treatment_effect_subset$Group],
     pch = 19, cex = 1.,
     xlim = xlim_common, ylim = ylim_common,
     xlab = "LKDPI", ylab = "Graft Lifespan",
     main = "",
     cex.lab = 1.5, cex.axis = 1.2,
     font.lab = 2, font.axis = 1,
     bty        = "l" )


frontiers_plotted <- list()  # Reset for second plot

for(i in hyperbolic_indices) {
  eval_row <- treatment_effect_subset[i, ]
  
  # Get local reference set for this evaluation point
  local_ref <- get_local_reference(eval_row, reference_data, cont_vars, cat_vars)
  
  print(length(local_ref$indices))
  
  if(length(local_ref$indices) >= 3) {  # Need at least 3 points for a meaningful frontier
    local_ref_data <- local_ref$reference_subset

    ref_id <- paste(sort(local_ref$indices), collapse = "_")
    if(!(ref_id %in% names(frontiers_plotted))) {
      frontier <- create_conditional_frontier(local_ref_data$QualityScore, 
                                              local_ref_data$OutcomeScore,
                                              xlim_common, ylim_common)
      
      if(!is.null(frontier)) {
        # Plot the conditional frontier
        lines(frontier$x, frontier$y, col = "blue", lwd = 3, lty = 1)
        frontiers_plotted[[ref_id]] <- TRUE
      }
    }
    
    # Calculate conditional hyperbolic path
    path <- calculate_conditional_hyperbolic_path(
      treatment_effect_subset$QualityScore[i],
      treatment_effect_subset$OutcomeScore[i],
      local_ref_data$QualityScore,
      local_ref_data$OutcomeScore
    )
    # Plot hyperbolic path 
    lines(path$X, path$Y, col = "darkgray", lty = 2, lwd = 2)
    points(path$X[length(path$X)], path$Y[length(path$Y)], pch = 8, col = "red", cex = 2)
  }
}
# Add legend
par(fig = c(0, 1, 0, 1), oma = c(8, 0, 0, 0), mar = c(0, 0, .5, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", 
       legend = c("Asian", "Black", "Hispanic", "White", "Projection"),
       col = c(color_palette, "red"),
       pch = c(rep(19, 4), 8),
       title = "Group", horiz = TRUE, xpd = TRUE, inset = c(0, 0), 
       cex = 1.15, 
       text.font = 2, 
       bty = "n", 
       x.intersp  = 0.3,
       pt.cex = 1) 

