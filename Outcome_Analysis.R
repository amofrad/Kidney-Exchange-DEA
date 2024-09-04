library(survival)
library(cmprsk)
library(dplyr)
library(purrr)
library(patchwork)
library(ggplot2)
library(tidyr)
library(dplyr)

balanced_data <- read.csv("balanced_data.csv")


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



################################################################################
################### CUMULATIVE INCIDENCE FUNCTION PLOT (Fig. 3) ################
################################################################################
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
