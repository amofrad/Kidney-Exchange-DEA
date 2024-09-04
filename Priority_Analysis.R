
library(mma)
library(dplyr)

balanced_data <- read.csv("balanced_data.csv")

data <- dplyr::select(balanced_data, WL_days, PRA, ABO_A, ABO_B, ABO_AB,
                      ABO_O,ABO,ETH_Asian, ETH_Black, ETH_Hispanic, ETH_White,
                      AGE, GENDER, ON_DIALYSIS, ETHCAT, ABO, GTIME_KI, 
                      ETHCAT_DON, MED_COND_TRR, KDPI, GFR, AGE_DON, 
                      BMI_DON_CALC, CREAT_DON, COD_CAD_DON, ETHCAT_DON,
                      DIABETES_DON, HGT_CM_DON_CALC, WGT_KG_DON_CALC,
                      BMI_CALC, COLD_ISCH_KI, HLAMIS)

data$ABO_A <- as.factor(data$ABO_A)
data$ABO_B <- as.factor(data$ABO_B)
data$ABO_AB <- as.factor(data$ABO_AB)
data$ABO_O <- as.factor(data$ABO_O)
data$ETHCAT <- factor(data$ETHCAT, 
                      levels = c("White", "Asian", "Black", "Hispanic"))
data$ETH_Asian <- as.factor(data$ETH_Asian)
data$ETH_Black <- as.factor(data$ETH_Black)
data$ETH_Hispanic <- as.factor(data$ETH_Hispanic)
data$ETH_White <- as.factor(data$ETH_White)
data$HLAMIS <- as.factor(data$HLAMIS)
data$ABO <- factor(data$ABO, levels = c("O", "A", "B", "AB"))


summary(lm(WL_days ~ PRA + ABO + ON_DIALYSIS + AGE + GENDER, data = data))
summary(aov(PRA ~ ETH_Asian, data = data))
summary(aov(PRA ~ ETH_Black, data = data))
summary(aov(PRA ~ ETH_Hispanic, data = data))
summary(aov(AGE ~ ETH_Asian, data = data))
summary(aov(AGE ~ ETH_Black, data = data))
summary(aov(AGE ~ ETH_Hispanic, data = data))
chisq.test(data$ABO_A, data$ETH_Asian)
chisq.test(data$ABO_A, data$ETH_Black)
chisq.test(data$ABO_A, data$ETH_Hispanic)
chisq.test(data$ABO_B, data$ETH_Asian)
chisq.test(data$ABO_B, data$ETH_Black)
chisq.test(data$ABO_B, data$ETH_Hispanic)
chisq.test(data$ABO_AB, data$ETH_Asian)
chisq.test(data$ABO_AB, data$ETH_Black)
chisq.test(data$ABO_AB, data$ETH_Hispanic)



covs <- data.frame(data[, c(12,13)])
x <- data[, c(2,7, 12:14)]
y <- data.frame(data[,1])
x$ABO <- factor(x$ABO, levels = c("O", "A", "B", "AB"))
pred_ETH <- data$ETHCAT

x <- as.data.frame(x)

mma_ETH <- mma(x, y, pred_ETH, mediator = c(1,2,5), catmed = c(2), 
               catref = "O", para = F, n2=1000)

summary(mma_ETH)