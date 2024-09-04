# Fairness in Kidney Exchange Programs: A DEA Approach

This repository contains the code used for the analysis in the paper "A Data Envelopment Analysis Approach for Assessing Fairness in Resource Allocation: Application to Kidney Exchange Programs".

## Code Files and Corresponding Paper Sections

### Data_Preprocessing.R
This script performs exploratory data analysis (Section 3.1), generating summary statistics for key fairness criteria across ethnic groups. It also implements data resampling and debiasing techniques (Section 3.2), addressing data imbalance through stratified sampling and applying a debiasing procedure to account for potential confounding variables.

### Priority_Analysis.R
This code conducts mediation analysis for the Priority Fairness criteria (Section 3.3). It decomposes the total effect of ethnicity into direct and indirect effects to evaluate the direct effect of ethnicity on waitlist duration.

### Access_Analysis.R
This script performs counterfactual analysis for KDPI scores (Section 3.3) to evaluate the influence of ethnicity in estimating KDPI scores. It employs a random forest algorithm and conducts variable importance analysis using out-of-bag permutation importance method.

### Outcome_Analysis.R
This code implements competing risks analysis for the Outcome Fairness criteria (Section 3.3). It evaluates the effect of ethnicity on graft lifespan by calculating hazard ratios for graft rejection while accounting for other causes of graft failure.

### DEA_Analysis.R
This script contains the Data Envelopment Analysis (DEA) model used in the study (methodology described in Section 2.2, results in Section 4.1). It implements the hyperbolic graph efficiency measure, structuring the model with waitlist duration and KDPI score as inputs, and graft lifespan as the output. The script calculates efficiency scores for each patient and analyzes efficiency distributions across ethnic groups.

### GroupConditional_Conformal_Predictions.py
This code performs uncertainty quantification on the DEA efficiency scores using the conformal prediction method developed by Gibbs et al. (2023) (methodology in Section 2.3, results in Section 4.1). It generates prediction intervals to quantify uncertainty for the efficiency scores, providing group-conditional coverage guarantees.

### MMD_Testing.R
This script conducts pairwise hypothesis testing for DEA efficiency scores between ethnic groups using the Maximum Mean Discrepancy (MMD) test (Section 4.2). It implements a kernel-based approach to assess the statistical significance of observed differences in efficiency score distributions across ethnic groups.
