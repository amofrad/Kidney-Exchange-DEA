library(Benchmarking)
library(dplyr)
library(ggplot2)
library(ggridges)

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
######### Plotting Distributions of Efficiency Scores by Group (Fig. 4) ########
################################################################################


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
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14, face = "bold"), 
    axis.title.x = element_text(size = 16, face = "bold")  
  )


################################################################################
############ Plotting Production Possibility Set + Frontier (Fig. 1) ###########
################################################################################
# Convert group numbers to labels
treatment_effect <- treatment_effect %>%
  mutate(Group = case_when(
    Group == "1" ~ "Asian",
    Group == "2" ~ "Black",
    Group == "3" ~ "Hispanic",
    Group == "4" ~ "White",
    TRUE ~ as.character(Group)
  ))
treatment_effect$Group <- as.factor(treatment_effect$Group)

# Define color palette
color_palette <- c(
  "Asian" = "#08306b",
  "Black" = "#1f78b4",
  "Hispanic" = "#6baed6",
  "White" = "#b3cde3"
)

# Function to calculate points along the hyperbolic path
calculate_hyperbolic_path <- function(X, Y, eff, steps = 100) {
  t <- seq(1, eff, length.out = steps)
  path_X <- X * t
  path_Y <- Y / t
  return(list(X = path_X, Y = path_Y))
}

# Set seed and select a subset of 50 observations
set.seed(14)
subset_indices <- sample(which(treatment_effect$efficiency <= 1), 200)
treatment_effect_subset <- treatment_effect[subset_indices, ]

# Select 10% of the subset for hyperbolic path lines
hyperbolic_indices <- sample(1:nrow(treatment_effect_subset), size = ceiling(0.1 * nrow(treatment_effect_subset)))

# Save plot as PNG
png("frontier.png", width = 1200, height = 600)

# Adjust margins to accommodate legend at the bottom
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1, oma = c(0, 0, 4, 0))  # Increased oma to 4 for more space

# Set axis limits to [0, 1] for both axes
xlim_common <- c(0, 1)
ylim_common <- c(0, 1)

# DEA plot for Waitlist Duration vs Outcome Score
dea.plot(treatment_effect_subset$WaitlistDuration, 
         treatment_effect_subset$OutcomeScore,
         RTS = "vrs",
         ORIENTATION = "graph",
         main = "",
         xlab = "Waitlist Duration",
         ylab = "Graft Lifespan", 
         pch = NA, 
         xlim = xlim_common,
         ylim = ylim_common,
         cex.lab = 2,  # Further increase label size
         cex.axis = 1.5, # Further increase axis tick label size
         font.lab = 2,   # Make label text bold
         font.axis = 2, lwd = 3)  # Make axis text bold

# Add points colored by ethnicity
points(treatment_effect_subset$WaitlistDuration, 
       treatment_effect_subset$OutcomeScore, 
       col = color_palette[treatment_effect_subset$Group],
       pch = 19, cex = 2)  # Increase point size

# Perform DEA
dea_result_WL <- dea(X = matrix(treatment_effect_subset$WaitlistDuration, ncol = 1),
                     Y = matrix(treatment_effect_subset$OutcomeScore, ncol = 1),
                     RTS = "vrs",
                     ORIENTATION = "graph", FAST = FALSE)

# Add hyperbolic path lines for 10% of the sampled points
for(i in hyperbolic_indices) {
  path <- calculate_hyperbolic_path(treatment_effect_subset$WaitlistDuration[i],
                                    treatment_effect_subset$OutcomeScore[i],
                                    dea_result_WL$eff[i])
  
  lines(path$X, path$Y, col = "darkgray", lty = 2, lwd = 2)
  points(path$X[length(path$X)], path$Y[length(path$Y)], pch = 8, col = "red", cex = 2)
}

# DEA plot for KDPI vs Outcome Score
dea.plot(treatment_effect_subset$QualityScore, 
         treatment_effect_subset$OutcomeScore,
         RTS = "vrs",
         ORIENTATION = "graph",
         main = "",
         xlab = "KDPI",
         ylab = "Graft Lifespan", 
         pch = NA, 
         xlim = xlim_common,
         ylim = ylim_common,
         cex.lab = 2,  # Further increase label size
         cex.axis = 1.5, # Further increase axis tick label size
         font.lab = 2,   # Make label text bold
         font.axis = 2, lwd = 3)  # Make axis text bold

# Add points colored by ethnicity
points(treatment_effect_subset$QualityScore, 
       treatment_effect_subset$OutcomeScore, 
       col = color_palette[treatment_effect_subset$Group],
       pch = 19, cex = 2)  # Increase point size

# Perform DEA
dea_result_KDPI <- dea(X = matrix(treatment_effect_subset$QualityScore, ncol = 1),
                       Y = matrix(treatment_effect_subset$OutcomeScore, ncol = 1),
                       RTS = "vrs",
                       ORIENTATION = "graph", FAST = FALSE)

# Add hyperbolic path lines for 10% of the sampled points
for(i in hyperbolic_indices) {
  path <- calculate_hyperbolic_path(treatment_effect_subset$QualityScore[i],
                                    treatment_effect_subset$OutcomeScore[i],
                                    dea_result_KDPI$eff[i])
  
  lines(path$X, path$Y, col = "darkgray", lty = 2, lwd = 2)
  points(path$X[length(path$X)], path$Y[length(path$Y)], pch = 8, col = "red", cex = 2)
}

# Add a legend under both plots
par(fig = c(0, 1, 0, 1), oma = c(2, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", 
       legend = c("Asian", "Black", "Hispanic", "White", "Projection"),
       col = c(color_palette, "red"),
       pch = c(rep(19, 4), 8),
       title = "Group", horiz = TRUE, xpd = TRUE, inset = c(0, -0.01),  # Adjusted inset for proper spacing
       cex = 2, # Further increase legend text size
       text.font = 2, # Bold legend text
       bty = "n")  # Remove the legend box outline

# Close the PNG device
dev.off()



