# Fairness in Kidney Exchange Programs: A DEA Approach

This repository contains the R code used for the analysis in the paper "A Data Envelopment Analysis Approach for Assessing Fairness in Resource Allocation: Application to Kidney Exchange Programs".


**Paper:** [https://arxiv.org/abs/2410.02799](https://arxiv.org/abs/2410.02799)


## Data

The empirical analysis utilizes data from the United Network for Organ Sharing (UNOS), specifically the kidney/pancreas transplant registry and living donor database. Due to data use agreements and patient privacy protections, the original datasets cannot be distributed with the code. To support reproducibility, we provide a version of the processed data with sensitive information removed. These files are available for download in `data_files.zip`, included in this repository.

## Code Files and Corresponding Paper Sections

### Data_Preprocessing.R
This script performs exploratory data analysis (Sections 3.1-3.2), generating summary statistics for key fairness criteria across ethnic groups. It computes LKDPI scores by matching transplant observations with donor characteristics, implements stratified resampling to match ESRD prevalence statistics, and applies year-specific centering of fairness measures. The script produces Tables 1-3 and Figure 2, and generates `DEA_data.csv` for subsequent analyses.

### Priority_Analysis.R
This code conducts mediation analysis for the Priority Fairness criteria (Section 3.3). It employs likelihood ratio testing for mediator selection and decomposes the total effect of ethnicity on waitlist duration into direct and indirect pathways using bootstrap confidence intervals. The analysis produces Table 4 summarizing mediation effects by ethnic group.

### Access_Analysis.R
This script performs counterfactual analysis for LKDPI scores (Section 3.3) to evaluate the influence of recipient ethnicity on organ quality allocation. It employs a random forest model with 5-fold cross-validation and conducts variable importance analysis using permutation-based methods. The analysis generates Table 5 documenting LKDPI score disparities across ethnic groups.

### Outcome_Analysis.R
This code implements competing risks analysis for the Outcome Fairness criteria (Section 3.3). It applies Fine-Gray subdistribution hazard models to evaluate ethnic disparities in graft rejection risk while accounting for competing causes of graft failure. The script uses sum-to-zero contrasts for ethnicity and generates Table 6 and Figure 3.

### Frontier_Plot.R
This script creates the illustrative visualization of the conditional DEA methodology (Section 2.2). It implements kernel-based localization procedures, constructs conditional efficiency frontiers for selected evaluation points, and visualizes hyperbolic efficiency trajectories. The script produces Figure 1 demonstrating the conditional DEA framework.

### RFM_Mapping.R
This script implements the core conditional DEA analysis with Reference Frontier Mapping (Sections 2.3 and 4.1). It handles data splitting to separate frontier construction from efficiency evaluation, implements kernel-based conditional production sets, and computes efficiency scores using targeted sampling strategies. The script generates efficiency datasets for three time periods and produces Figure 4 showing efficiency score distributions.

### Conformal_Predictions.R
This code performs uncertainty quantification on the conditional DEA efficiency scores using group-conditional conformal prediction (Sections 2.3 and 4.1). It employs linear programming formulation to optimize conformal scores while maintaining group-conditional coverage guarantees through stratified train-calibration-test splitting. The analysis conducts 100 independent splits and generates Table 7 with prediction intervals.

### MMD_test.R
This script conducts distributional hypothesis testing for DEA efficiency scores using Maximum Mean Discrepancy testing (Section 4.2). It employs Laplacian kernel-based two-sample tests to evaluate whether ethnic groups exhibit significantly different efficiency score distributions. The script performs both pairwise and group-versus-rest comparisons with Benjamini-Hochberg correction, generating Table 8.

### Supplemental_Resampling_Simulation.R
This code implements the simulation study described in Section S2 of the Supplementary Material. It generates synthetic datasets using multivariate normal distributions parameterized by empirical data characteristics and evaluates the impact of demographic resampling on conditional DEA efficiency scores. The script conducts sensitivity analysis across multiple demographic scenarios and produces Supplementary Tables 1-4 and Figure 1.

## Execution Order

Run the scripts in the following sequence:

1. `Data_Preprocessing.R` (creates `DEA_data.csv`)
2. Individual fairness analyses (can be run in any order):
   - `Priority_Analysis.R`
   - `Access_Analysis.R` 
   - `Outcome_Analysis.R`
3. `RFM_Mapping.R` (creates DEA result files)
4. `Frontier_Plot.R`, `Conformal_Predictions.R`, `MMD_test.R` (can be run in any order)
5. `Supplemental_Resampling_Simulation.R` (for supplementary analysis)
