library(maotai)
library(kernlab)

# Load data result from RFM_Mapping.R
DEA_Result_2010_2013 = read.csv("DEA_Result_2010_2013.csv")
DEA_Result_2010_2016 = read.csv("DEA_Result_2010_2016.csv")
DEA_Result_2010_2019 = read.csv("DEA_Result_2010_2019.csv")


################################################################################
############################## Pairwise Tests ##################################
################################################################################

pairwise_mmd <- function(data){
  groups <- unique(data$Group)
  
  results <- list()
  # Loop through all combinations of group pairs
  for(i in 1:(length(groups) - 1)) {
    for(j in (i + 1):length(groups)) {
      group_i <- groups[i]
      group_j <- groups[j]
      
      # Subset data for each group
      theta_i <- (data[data$Group == group_i,])
      theta_j <- (data[data$Group == group_j,])
      
      theta_i_mat <- matrix(theta_i$Relative_Eff, ncol = 1)
      theta_j_mat <- matrix(theta_j$Relative_Eff, ncol = 1)
      
      combined_data <- rbind(theta_i_mat, theta_j_mat)
      
      # Compute the Euclidean distance matrix
      dmat <- as.matrix(dist(combined_data))
      
      sigma <- median(dmat)
      
      laplace <- laplacedot(sigma = sigma)
      kmat <- kernelMatrix(laplace, combined_data)
      
      # Create labels
      lab <- c(theta_i$Group, theta_j$Group)

      set.seed(496)
      suppressWarnings(mmd_result <- mmd2test(kmat, lab, mc.iter = 500, method = "u"))

      cat("Group ", group_i, " vs Group ", group_j, ": p=", mmd_result$p.value, "\n")
      
      results[[paste0("Group_", group_i, "_vs_Group_", group_j)]] <- list(
        mmd_statistic = mmd_result$statistic,
        pvalue = mmd_result$p.value
      )
    }
  }
  return(results)
}


pairs_2010_2013 <- pairwise_mmd(DEA_Result_2010_2013)
pairs_2010_2016 <- pairwise_mmd(DEA_Result_2010_2016)
pairs_2010_2019 <- pairwise_mmd(DEA_Result_2010_2019)

MMD_2010_2013 <- sapply(pairs_2010_2013, function(x) x$mmd_statistic)
MMD_2010_2016 <- sapply(pairs_2010_2016, function(x) x$mmd_statistic)
MMD_2010_2019 <- sapply(pairs_2010_2019, function(x) x$mmd_statistic)

# Extract p-values for each
p_2010_2013 <- sapply(pairs_2010_2013, function(x) x$pvalue)
p_2010_2016 <- sapply(pairs_2010_2016, function(x) x$pvalue)
p_2010_2019 <- sapply(pairs_2010_2019, function(x) x$pvalue)

# Adjust p-values using Benjamini-Hochberg
adj_2010_2013 <- p.adjust(p_2010_2013, method = "BH")
adj_2010_2016 <- p.adjust(p_2010_2016, method = "BH")
adj_2010_2019 <- p.adjust(p_2010_2019, method = "BH")

(summary_table_pairwise <- data.frame(
  Stat_2010_2013 = round(MMD_2010_2013, 8),
  Stat_2010_2016 = round(MMD_2010_2016, 8),
  Stat_2010_2019 = round(MMD_2010_2019, 8),
  P_2010_2013 = round(p_2010_2013, 4),
  Adj_2010_2013 = round(adj_2010_2013, 4),
  P_2010_2016 = round(p_2010_2016, 4),
  Adj_2010_2016 = round(adj_2010_2016, 4),
  P_2010_2019 = round(p_2010_2019, 4),
  Adj_2010_2019 = round(adj_2010_2019, 4)
))


################################################################################
############################## Group-vs-Rest Tests #############################
################################################################################
group_vs_rest_mmd <- function(data){
  groups <- unique(data$Group)
  results <- list()
  
  # Loop through all combinations of group pairs
  for(i in seq_along(groups)) {
    group_i <- groups[i]
    
    
    # Subset data for each group
    theta_i <- (data[data$Group == group_i,])
    theta_rest <- (data[data$Group != group_i,])
    
    theta_i_mat <- matrix(theta_i$Relative_Eff, ncol = 1)
    theta_rest_mat <- matrix(theta_rest$Relative_Eff, ncol = 1)
    
    combined_data <- rbind(theta_i_mat, theta_rest_mat)
  
    # Compute the Euclidean distance matrix
    dmat <- as.matrix(dist(combined_data))
    
    sigma <- median((dmat))

    laplace <- laplacedot(sigma = sigma)
    kmat <- kernelMatrix(laplace, combined_data)
    
    # Create labels
    lab <- c(theta_i$Group, rep("Rest", length(theta_rest_mat)))
    set.seed(480) # for Laplace
    
    
    suppressWarnings(mmd_result <- mmd2test(kmat, lab, mc.iter = 500, method = "u"))

    cat("Group", group_i, " vs Rest: p=", mmd_result$p.value, "\n")
    
    results[[paste0("Group_", group_i, "_vs_Rest")]] <- list(
      mmd_statistic = mmd_result$statistic,
      pvalue = mmd_result$p.value
    )
  }
  
  return(results)
}



group_rest_2010_2013 <- group_vs_rest_mmd(DEA_Result_2010_2013)
group_rest_2010_2016 <- group_vs_rest_mmd(DEA_Result_2010_2016)
group_rest_2010_2019 <- group_vs_rest_mmd(DEA_Result_2010_2019)


group_MMD_2010_2013 <- sapply(group_rest_2010_2013, function(x) x$mmd_statistic)
group_MMD_2010_2016 <- sapply(group_rest_2010_2016, function(x) x$mmd_statistic)
group_MMD_2010_2019 <- sapply(group_rest_2010_2019, function(x) x$mmd_statistic)

# Extract p-values for each
p_rest_2010_2013 <- sapply(group_rest_2010_2013, function(x) x$pvalue)
p_rest_2010_2016 <- sapply(group_rest_2010_2016, function(x) x$pvalue)
p_rest_2010_2019 <- sapply(group_rest_2010_2019, function(x) x$pvalue)

# Adjust p-values using Benjamini-Hochberg
adj_rest_2010_2013 <- p.adjust(p_rest_2010_2013, method = "BH")
adj_rest_2010_2016 <- p.adjust(p_rest_2010_2016, method = "BH")
adj_rest_2010_2019 <- p.adjust(p_rest_2010_2019, method = "BH")

(summary_table <- data.frame(
  Stat_2010_2013 = round(group_MMD_2010_2013, 8),
  Stat_2010_2016 = round(group_MMD_2010_2016, 8),
  Stat_2010_2019 = round(group_MMD_2010_2019, 8),
  P_2010_2013 = round(p_rest_2010_2013, 4),
  Adj_2010_2013 = round(adj_rest_2010_2013, 4),
  P_2010_2016 = round(p_rest_2010_2016, 4),
  Adj_2010_2016 = round(adj_rest_2010_2016, 4),
  P_2010_2019 = round(p_rest_2010_2019, 4),
  Adj_2010_2019 = round(adj_rest_2010_2019, 4)
))
