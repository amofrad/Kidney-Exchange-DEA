library(dplyr)
library(purrr)
library(lpSolveAPI)
library(tidyverse)
library(dplyr)

# Load data result from RFM_Mapping.R
DEA_Result_2010_2013 = read.csv("DEA_Result_2010_2013.csv")
DEA_Result_2010_2016 = read.csv("DEA_Result_2010_2016.csv")
DEA_Result_2010_2019 = read.csv("DEA_Result_2010_2019.csv")

run_conformal_simulation <- function(df, n_splits = 100) {
  # Data preparation
  df$Group <- factor(df$Group)
  df_dummies <- model.matrix(~ Group - 1, df) %>%
    as.data.frame() %>%
    bind_cols(df %>% dplyr::select(-Group)) %>%
    rename(
      Group_Asian    = GroupAsian,
      Group_Black    = GroupBlack,
      Group_Hispanic = GroupHispanic,
      Group_White    = GroupWhite
    )
  
  features <- c("PriorityScore", "QualityScore", "OutcomeScore", "Relative_Eff",
                "Group_Asian", "Group_Black", "Group_Hispanic", "Group_White", "Year")
  dataSub <- df_dummies %>% dplyr::select(all_of(features))
  
  # Prepare X and Y
  X <- dataSub %>% dplyr::select(-Relative_Eff)
  Y <- dataSub$Relative_Eff
  
  # Extract protected group indicators
  phiFn <- function(X) {
    X %>% dplyr::select(starts_with("Group_")) %>% as.matrix()
  }
  
  # Single-simulation function
  run_single_simulation <- function(X, Y, seed) {
    set.seed(seed)
    
    # Split data into training, calibration, and test sets
    n <- nrow(X)
    test_prop <- 0.35
    cal_prop <- 0.5
    
    # Get group column names
    group_cols <- grep("^Group_", colnames(X), value = TRUE)
    
    # Determine group label for each row
    group_labels <- apply(X[, group_cols], 1, function(row) which(row == 1))
    group_names <- gsub("Group_", "", group_cols)
    group_factor <- factor(group_labels, labels = group_names)
    
    calib_idx <- c()
    test_idx <- c()
    
    for (g in levels(group_factor)) {
      group_rows <- which(group_factor == g)
      n_g <- length(group_rows)
      nTest_g <- round(test_prop * n_g)
      nCal_g <- round(cal_prop * n_g)
      
      shuffled <- sample(group_rows)
      test_idx <- c(test_idx, shuffled[1:nTest_g])
      calib_idx <- c(calib_idx, shuffled[(nTest_g + 1):(nTest_g + nCal_g)])
    }
    
    # Training is everything else
    idx_all <- seq_len(n)
    train_idx <- setdiff(idx_all, union(test_idx, calib_idx))
    
    XTrain <- X[train_idx, ]; YTrain <- Y[train_idx]
    XCalib <- X[calib_idx, ]; YCalib <- Y[calib_idx]
    XTest <- X[test_idx, ]; YTest <- Y[test_idx]
    
    # Print split information
    cat("Split sizes - Train:", length(train_idx),
        "Calib:", length(calib_idx),
        "Test:", length(test_idx), "\n")
    
    group_cols <- grep("^Group_", colnames(X), value = TRUE)
    for (g in group_cols) {
      train_prop <- mean(X[train_idx, g])
      calib_prop <- mean(X[calib_idx, g])
      test_prop <- mean(X[test_idx, g])
      cat(sprintf("%s - Train: %.2f, Calib: %.2f, Test: %.2f\n",
                  g, train_prop, calib_prop, test_prop))
    }
    cat("======================================================\n")
    
    # Fit regression model
    reg <- lm(YTrain ~ PriorityScore + QualityScore + OutcomeScore + 
                Group_Asian + Group_Black + Group_Hispanic + Group_White + 
                as.factor(Year), data = XTrain)
    
    scoreFn <- function(x, y) abs(y - predict(reg, newdata = x))
    alpha <- 0.05
    
    # Calibration scores & phi
    scores_calib <- scoreFn(XCalib, YCalib)
    phi_calib <- phiFn(XCalib)
    svd_res <- svd(phi_calib, nu = min(dim(phi_calib)), nv = min(dim(phi_calib)))
    tol <- 1e-5
    r <- sum(svd_res$d > tol)
    
    if(r < length(svd_res$d)) {
      Tmat <- svd_res$v[, 1:r, drop = FALSE]
      phiFn_reduced <- function(x) phiFn(x) %*% Tmat
      Phi_calib_reduced <- phi_calib %*% Tmat
    } else {
      phiFn_reduced <- phiFn
      Phi_calib_reduced <- phi_calib
    }
    
    # LP solver functions
    get_eta <- function(S, x) {
      # Rebuild full-vector and design-matrix
      phi_x <- phiFn_reduced(x)                 # 1×d
      S_full <- c(scores_calib, S)               # (n_calib+1)
      Phi_full <- rbind(Phi_calib_reduced, phi_x)  # (n_calib+1)×d
      n_full <- length(S_full)
      
      # Create empty LP with n_full variables
      m <- make.lp(0, n_full)
      
      # Add equality constraints (one per column of Phi_full)
      for (j in seq_len(ncol(Phi_full))) {
        add.constraint(m,
                       xt = Phi_full[, j],
                       type = "=",
                       rhs = 0)
      }
      
      # Set exactly −alpha <= eta_i <= 1−alpha
      set.bounds(m,
                 lower = rep(-alpha, n_full),
                 upper = rep(1-alpha, n_full),
                 columns = 1:n_full)
      
      # Objective
      set.objfn(m, -S_full)
      
      # Solve
      if (solve(m) != 0) stop("LP failed")
      sol <- get.variables(m)
      
      # Equality-constraints:
      viol <- t(Phi_full) %*% sol
      stopifnot(max(abs(viol)) < 1e-6)
      
      # box-constraints:
      stopifnot(all(sol >= -alpha - 1e-8), all(sol <= 1-alpha + 1e-8))
      
      sol
    }
    
    find_S_star <- function(x, randomize = TRUE) {
      U <- if (randomize) runif(1, -alpha, 1-alpha) else 1-alpha
      f <- function(S) {
        eta <- get_eta(S, x)
        eta[length(eta)] - U
      }
      Smin <- 0
      Smax <- max(scores_calib)
      fmin <- f(Smin)
      fmax <- f(Smax)
      
      # If eta_{n+1}(S) never crosses U in [0, max(scores_calib)], just return max(scores)
      if (fmin*fmax > 0) {
        return(Smax)
      }
      
      # Otherwise bracket is valid
      uniroot(f, c(Smin, Smax))$root
    }
    
    # Verify coverage & build intervals
    verify_coverage <- function(Xt, Yt) {
      m <- nrow(Xt)
      covs <- logical(m)
      preds <- vector("list", m)
      for(i in seq_len(m)) {
        xi <- Xt[i, , drop = FALSE]
        yi <- Yt[i]
        S0 <- scoreFn(xi, yi)
        eta0 <- get_eta(S0, xi)
        U <- runif(1, -alpha, 1-alpha)
        covs[i] <- (eta0[length(eta0)] < U)
        mu <- predict(reg, xi)
        Sstar <- find_S_star(xi)
        preds[[i]] <- c(mu - Sstar, mu + Sstar)
      }
      list(coverage = mean(covs),
           predictions = preds)
    }
    
    overall <- verify_coverage(XTest, YTest)
    group_cols <- paste0("Group_", levels(df$Group))
    group_res <- map(group_cols, ~{
      sel <- XTest[[.x]] == 1
      if(any(sel)) verify_coverage(XTest[sel, ], YTest[sel])
      else list(coverage = NA, predictions = list())
    }) %>% set_names(group_cols)
    
    cat(sprintf(
      "Split coverage:\n  Asian:    %.4f\n  Black:    %.4f\n  Hispanic: %.4f\n  White:    %.4f\n",
      group_res$Group_Asian$coverage,
      group_res$Group_Black$coverage,
      group_res$Group_Hispanic$coverage,
      group_res$Group_White$coverage
    ))
    
    list(
      overall_coverage = overall$coverage,
      overall_predictions = overall$predictions,
      group_results = group_res
    )
  }
  
  # Main simulation loop
  overall_cov <- numeric(n_splits)
  overall_preds <- list()
  group_cols <- paste0("Group_", levels(df$Group))
  group_res <- map(group_cols, ~list(coverages = numeric(), predictions = list())) %>% 
    set_names(group_cols)
  
  for (f in seq(n_splits)) {
    cat("\n======================================================\n")
    cat("======================================================\n")
    cat("Processing Split", f, "of", n_splits, "...\n")
    res <- run_single_simulation(X, Y, seed = f)
    overall_cov[f] <- res$overall_coverage
    overall_preds <- c(overall_preds, res$overall_predictions)
    cat("Mean Splits Coverage:\n")
    for(g in group_cols) {
      grp <- res$group_results[[g]]
      if(!is.na(grp$coverage)) {
        group_res[[g]]$coverages <- c(group_res[[g]]$coverages, grp$coverage)
        group_res[[g]]$predictions <- c(group_res[[g]]$predictions, grp$predictions)
      }
      cat(sprintf("  %-9s %.4f\n", gsub("Group_", "", g), mean(group_res[[g]]$coverages, na.rm = TRUE)))
    }
  }
  list(
    overall_coverages = overall_cov,
    overall_predictions = overall_preds,
    group_results = group_res
  )
}

conformal_Result_2010_2013 <- run_conformal_simulation(df = DEA_Result_2010_2013, n_splits = 100)
conformal_Result_2010_2016 <- run_conformal_simulation(df = DEA_Result_2010_2016, n_splits = 100)
conformal_Result_2010_2019 <- run_conformal_simulation(df = DEA_Result_2010_2019, n_splits = 100)



display_results <- function(sim_res){
  cat("\nAggregate Results:\n")
  cat(sprintf("Overall - Mean Coverage: %.4f, Std: %.4f\n",
              mean(sim_res$overall_coverages),
              sd(sim_res$overall_coverages)))
  for(g in names(sim_res$group_results)) {
    cvec <- sim_res$group_results[[g]]$coverages
    cat(sprintf("%s - Mean Coverage: %.4f, Std: %.4f\n",
                g, mean(cvec, na.rm=TRUE), sd(cvec, na.rm=TRUE)))
  }
  
  cat("\nPrediction Intervals:\n")
  overall_means <- map_dbl(sim_res$overall_predictions, mean)
  overall_lo    <- map_dbl(sim_res$overall_predictions, ~.[1])
  overall_hi    <- map_dbl(sim_res$overall_predictions, ~.[2])
  cat(sprintf("Overall - Mean: %.4f, 95%% CI: [%.4f, %.4f]\n",
              mean(overall_means),
              quantile(overall_lo, 0.025),
              quantile(overall_hi, 0.975)))
  for(g in names(sim_res$group_results)) {
    preds <- sim_res$group_results[[g]]$predictions
    if(length(preds)>0) {
      mvec <- map_dbl(preds, mean)
      lo   <- map_dbl(preds, ~.[1])
      hi   <- map_dbl(preds, ~.[2])
      cat(sprintf("%s - Mean: %.4f, 95%% CI: [%.4f, %.4f]\n",
                  g, mean(mvec), quantile(lo,0.025), quantile(hi,0.975)))
    }
  }
}

display_results(conformal_Result_2010_2013)
display_results(conformal_Result_2010_2016)
display_results(conformal_Result_2010_2019)
