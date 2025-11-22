library(dplyr)
library(tidyr)
library(ggplot2)
library(Benchmarking)
library(tidyverse)


set.seed(1480)

DEA_data <- read.csv("DEA_data.csv")

frontier_plot_data <- data.frame(
  Group = DEA_data$Group,
  Year = DEA_data$Year,
  PriorityScore = DEA_data$WaitlistDuration,
  QualityScore = DEA_data$LKDPI,
  OutcomeScore = DEA_data$GTIME_KI,
  EDUCATION = as.factor(DEA_data$EDUCATION),
  REGION = as.factor(DEA_data$REGION),
  CITIZENSHIP = as.factor(DEA_data$CITIZENSHIP),
  AGE = DEA_data$AGE,
  GENDER = DEA_data$GENDER
)

frontier_plot_data_2010_2019 <- frontier_plot_data[frontier_plot_data$Year <= 2019,]


DEA_Result_2010_2019 = read.csv("DEA_Result_2010_2019.csv")
# Prepare the data - Convert group numbers to labels
DEA_Result_2010_2019 <- DEA_Result_2010_2019 %>%
  mutate(Group = case_when(
    Group == "1" ~ "Asian",
    Group == "2" ~ "Black", 
    Group == "3" ~ "Hispanic",
    Group == "4" ~ "White",
    TRUE ~ as.character(Group)
  ))
DEA_Result_2010_2019$Group <- as.factor(DEA_Result_2010_2019$Group)

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
frontier_plot_data_scaled <- min_max_scale(frontier_plot_data_2010_2019)

# Create reference set
get_diverse_reference <- function(data) {
  bind_rows(
    data %>% filter(PriorityScore <= quantile(PriorityScore, gamma)),
    data %>% filter(QualityScore <= quantile(QualityScore, gamma)),
    data %>% filter(OutcomeScore >= quantile(OutcomeScore, 1 - gamma))
  ) %>% distinct()
}

target_n_ref <- ceiling(ref_prop * nrow(frontier_plot_data_scaled))

diverse_reference <- get_diverse_reference(frontier_plot_data_scaled)
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
  
  suppressMessages(evaluation_data <- anti_join(frontier_plot_data_scaled, reference_data))
  
} else {
  n_remaining <- target_n_ref - n_diverse
  remaining_pool <- anti_join(frontier_plot_data, diverse_reference)
  
  additional_reference <- remaining_pool %>% slice_sample(n = n_remaining)
  reference_data <- bind_rows(diverse_reference, additional_reference)
  evaluation_data <- anti_join(frontier_plot_data_scaled, reference_data)
}

# Merge efficiency scores from DEA_Result_2010_2019 to the scaled evaluation data
evaluation_data$efficiency <- DEA_Result_2010_2019$efficiency[match(
  paste(evaluation_data$Group, evaluation_data$Year, evaluation_data$AGE),
  paste(DEA_Result_2010_2019$Group, DEA_Result_2010_2019$Year, DEA_Result_2010_2019$AGE)
)]

# Define variable splits
cont_vars <- c()  
cat_vars <- c("REGION", "EDUCATION", "CITIZENSHIP")


# --- Choose one focal patient from evaluation_data ---

focal_idx <- sample(seq_len(nrow(evaluation_data)), 1)


focal      <- evaluation_data[focal_idx, ]


# Get local reference set for this focal patient
local_ref <- get_local_reference(
  eval_row       = focal,
  reference_data = reference_data,
  cont_vars      = cont_vars,
  cat_vars       = cat_vars
)
local_ref_data <- local_ref$reference_subset

# Conditional frontiers + hyperbolic paths
frontier_wait <- create_conditional_frontier(
  local_ref_x = local_ref_data$PriorityScore,
  local_ref_y = local_ref_data$OutcomeScore,
  xlim        = c(0, 1), ylim = c(0, 1)
)

path_wait <- calculate_conditional_hyperbolic_path(
  X           = focal$PriorityScore,
  Y           = focal$OutcomeScore,
  local_ref_x = local_ref_data$PriorityScore,
  local_ref_y = local_ref_data$OutcomeScore
)

frontier_qual <- create_conditional_frontier(
  local_ref_x = local_ref_data$QualityScore,
  local_ref_y = local_ref_data$OutcomeScore,
  xlim        = c(0, 1), ylim = c(0, 1)
)

path_qual <- calculate_conditional_hyperbolic_path(
  X           = focal$QualityScore,
  Y           = focal$OutcomeScore,
  local_ref_x = local_ref_data$QualityScore,
  local_ref_y = local_ref_data$OutcomeScore
)




# -----------------------------
# Build plotting dataframe
# -----------------------------

# All patients
df_all <- frontier_plot_data_scaled %>%
  select(PriorityScore, QualityScore, OutcomeScore) %>%
  mutate(type = "All patients")

# Local reference set points
df_ref <- local_ref_data %>%
  transmute(
    PriorityScore, QualityScore, OutcomeScore,
    type = "Local reference set",
    Group = Group
  )

# Focal patient (recycled for both panels)
df_focal <- tibble(
  PriorityScore = focal$PriorityScore,
  QualityScore  = focal$QualityScore,
  OutcomeScore  = focal$OutcomeScore,
  Group         = focal$Group,
  type          = "Focal patient"
)

# Hyperbolic path (two panels)
df_path_wait <- tibble(
  X = path_wait$X,
  Y = path_wait$Y,
  type = "Hyperbolic path",
  panel = "Waitlist vs Graft Lifespan"
)

df_path_qual <- tibble(
  X = path_qual$X,
  Y = path_qual$Y,
  type = "Hyperbolic path",
  panel = "LKDPI vs Graft Lifespan"
)

# Projection point
df_proj <- bind_rows(
  tibble(
    X = tail(path_wait$X, 1),
    Y = tail(path_wait$Y, 1),
    type = "Projection",
    panel = "Waitlist vs Graft Lifespan"
  ),
  tibble(
    X = tail(path_qual$X, 1),
    Y = tail(path_qual$Y, 1),
    type = "Projection",
    panel = "LKDPI vs Graft Lifespan"
  )
)

# Frontier curves
df_frontier_wait <- tibble(
  X = frontier_wait$x,
  Y = frontier_wait$y,
  type = "Conditional frontier",
  panel = "Waitlist vs Graft Lifespan"
)

df_frontier_qual <- tibble(
  X = frontier_qual$x,
  Y = frontier_qual$y,
  type = "Conditional frontier",
  panel = "LKDPI vs Graft Lifespan"
)

# -----------------------------
# Assemble data for each panel
# -----------------------------
df_panel_wait <- tibble(
  X = c(df_all$PriorityScore),
  Y = c(df_all$OutcomeScore),
  type = "All patients",
  panel = "Waitlist vs Graft Lifespan"
)

df_panel_qual <- tibble(
  X = c(df_all$QualityScore),
  Y = c(df_all$OutcomeScore),
  type = "All patients",
  panel = "LKDPI vs Graft Lifespan"
)

# Combined all-patient frame
df_all_long <- bind_rows(df_panel_wait, df_panel_qual)

# Reference set long format
df_ref_long <- bind_rows(
  tibble(
    X = df_ref$PriorityScore,
    Y = df_ref$OutcomeScore,
    Group = df_ref$Group,
    type = "Local reference set",
    panel = "Waitlist vs Graft Lifespan"
  ),
  tibble(
    X = df_ref$QualityScore,
    Y = df_ref$OutcomeScore,
    Group = df_ref$Group,
    type = "Local reference set",
    panel = "LKDPI vs Graft Lifespan"
  )
)

# Focal patient long format
df_focal_long <- bind_rows(
  tibble(
    X = df_focal$PriorityScore,
    Y = df_focal$OutcomeScore,
    Group = df_focal$Group,
    type = "Focal patient",
    panel = "Waitlist vs Graft Lifespan"
  ),
  tibble(
    X = df_focal$QualityScore,
    Y = df_focal$OutcomeScore,
    Group = df_focal$Group,
    type = "Focal patient",
    panel = "LKDPI vs Graft Lifespan"
  )
)




# -----------------------------
# Combine for plotting
# -----------------------------
df_frontier <- bind_rows(df_frontier_wait, df_frontier_qual)
df_path     <- bind_rows(df_path_wait, df_path_qual)
df_proj     <- df_proj

# -----------------------------
# Aesthetics
# -----------------------------
color_palette <- c(
  "Asian"    = "#08306b",
  "Black"    = "#1f78b4",
  "Hispanic" = "#6baed6",
  "White"    = "#b1dff9"
)

panel_levels <- c("Waitlist vs Graft Lifespan", "LKDPI vs Graft Lifespan")

df_all_long    <- df_all_long    %>% mutate(panel = factor(panel, levels = panel_levels))
df_ref_long    <- df_ref_long    %>% mutate(panel = factor(panel, levels = panel_levels))
df_focal_long  <- df_focal_long  %>% mutate(panel = factor(panel, levels = panel_levels))
df_frontier    <- df_frontier    %>% mutate(panel = factor(panel, levels = panel_levels))
df_path        <- df_path        %>% mutate(panel = factor(panel, levels = panel_levels))
df_proj        <- df_proj        %>% mutate(panel = factor(panel, levels = panel_levels))


# -----------------------------
# Final Plot
# -----------------------------
p <- ggplot() +
  ## ----------------------------------------------------- ##
  ## 1. Geoms
  ## ----------------------------------------------------- ##
  # All patients
  geom_point(
    data  = df_all_long,
    aes(X, Y, color = "All patients"),
    alpha = 0.4, size = 1
  ) +
  
  # Local reference set (colored by Group, outlined in black)
  geom_point(
    data  = df_ref_long,
    aes(X, Y, fill = Group),
    shape = 21, size = 4, stroke = 0.6,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
  # Focal patient (diamond)
  geom_point(
    data  = df_focal_long,
    aes(X, Y, fill = Group,
        color = "Focal patient",
        shape = "Focal patient"),
    size  = 5, stroke = 1.2
  ) +
  
  # Conditional frontier (line)
  geom_line(
    data = df_frontier,
    aes(X, Y,
        color    = "Conditional frontier",
        linetype = "Conditional frontier"),
    linewidth = 2
  ) +
  
  # Frontier points 
  geom_point(
    data  = df_frontier,
    aes(X, Y,
        color = "Frontier point",
        shape = "Frontier point"),
    stroke = 1.4, size = 4
  ) +
  
  # Hyperbolic path (dashed)
  geom_path(
    data = df_path,
    aes(X, Y,
        color    = "Hyperbolic path",
        linetype = "Hyperbolic path"),
    linewidth = 1
  ) +
  
  # Projection (star)
  geom_point(
    data  = df_proj,
    aes(X, Y,
        color = "Projection",
        shape = "Projection"),
    size = 5, stroke = 1.2
  ) +
  
  # Facets for the two panels
  facet_wrap(~ panel, scales = "free_x") +
  
  ## ----------------------------------------------------- ##
  ## 2. Scales
  ## ----------------------------------------------------- ##
  # Row 1: ethnic groups (fill legend)
  scale_fill_manual(
    name   = "Ethnic group",
    values = color_palette,
    breaks = c("Asian", "Black", "Hispanic", "White")
  ) +
  
  # Row 2: components (color / linetype / shape all share *another* legend)
  scale_color_manual(
    name   = "Components",
    values = c(
      #"All patients"         = "grey70",
      "Conditional frontier" = "royalblue3",
      "Hyperbolic path"      = "black",
      "Projection"           = "firebrick3",
      "Frontier point"       = "brown",
      "Focal patient"        = "brown"
    )
  ) +
  scale_linetype_manual(
    name   = "Components",
    values = c(
      "Conditional frontier" = "solid",
      "Hyperbolic path"      = "22"
    )
  ) +
  scale_shape_manual(
    name   = "Components",
    values = c(
      "Focal patient"   = 23,  # diamond
      "Frontier point"  = 16,   # X
      "Projection"      = 8    # star
    )
  ) +
  
  ## ----------------------------------------------------- ##
  ## 3. Guides and theme
  ## ----------------------------------------------------- ##
  guides(
    # First row: ethnic groups
    fill = guide_legend(
      order = 1,
      nrow  = 1,
      override.aes = list(
        shape = 21,
        size  = 6,
        colour = "black",stroke = 1
      )
    ),
    # Second row: components (merged color/shape/linetype)
    color = guide_legend(
      order = 2,
      nrow  = 1,
      override.aes = list(
        linetype = c("solid", "blank", "blank", "22", "blank"),
        shape    = c( NA, 23, 16, NA, 8),
        size     = c(8, 6, 6, 8, 6),
        fill     = NA
      )
    ),
    linetype = "none",  # handled via color guide above
    shape    = "none"   # handled via color guide above
  ) +
  
  labs(
    x = "",
    y = "Graft Lifespan"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "top",
    legend.box       = "vertical",   # groups row above components row
    legend.title     = element_blank(),
    strip.text       = element_text(size = 16, face = "bold"),
    axis.title.x     = element_text(size = 12),
    axis.title.y     = element_text(size = 12),
    legend.key.size  = unit(0.9, "lines"),
    legend.text      = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

p
