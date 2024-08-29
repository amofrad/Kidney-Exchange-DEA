# -*- coding: utf-8 -*-

!pip install conditionalconformal # Install the 'conditionalconformal' package for conditional coverage prediction

import pandas as pd
import numpy as np
import cvxpy as cp
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.linear_model import LinearRegression
from scipy.optimize import linprog
from sklearn.metrics.pairwise import pairwise_kernels
from random import sample
import random
from sklearn.model_selection import KFold
from conditionalconformal import CondConf



# Load the data
df = pd.read_csv('treatment_effect_data.csv')

# Map the "Group" column to categorical names
mapping = {1: "Asian", 2: "Black", 3: "Hispanic", 4: "White"}
df['Group'] = df['Group'].map(mapping)

# Create dummy variables for each group
df_dummies = pd.get_dummies(df, columns=['Group'], prefix='Group', drop_first=False)

# Convert all columns that start with 'Group_' to integers
df_dummies.loc[:, df_dummies.columns.str.startswith('Group_')] = df_dummies.loc[:, df_dummies.columns.str.startswith('Group_')].astype(int)


# Select features for the model
features = ['WaitlistDuration', 'QualityScore', 'OutcomeScore', 'efficiency', 'Group_Asian', 'Group_Black', 'Group_Hispanic', 'Group_White']
dataSub = df_dummies[features]

# Prepare X and Y for modeling
X = dataSub.drop(['efficiency'], axis=1)
Y = dataSub['efficiency']


def phiFn(X):
    protectedFeatures = [3, 4, 5, 6]  # Indices for Group Indicators
    if X.ndim == 1:
        X = X.reshape(1, -1)
    phi = X[:, protectedFeatures]
    return phi

def run_single_simulation(X, Y, seed):
    np.random.seed(seed)

    nTest = round(X.shape[0] * 0.3)
    nCalib = round(X.shape[0] * 0.3)

    # Sample indices for the test set
    fullTrainingSet = np.random.choice(len(Y), len(Y) - nTest, replace=False)
    testPoints = list(set(range(0, len(Y))) - set(fullTrainingSet))

    # Split the full training set into calibration and training sets
    calibPoints = fullTrainingSet[:nCalib]
    trainPoints = fullTrainingSet[nCalib:]

    # Convert to NumPy arrays
    XTest = X.iloc[testPoints, :].values
    YTest = Y.iloc[testPoints].values
    XCalib = X.iloc[calibPoints, :].values
    YCalib = Y.iloc[calibPoints].values
    XTrain = X.iloc[trainPoints, :].values
    YTrain = Y.iloc[trainPoints].values

    reg = LinearRegression().fit(XTrain, YTrain)
    scoreFn = lambda x, y: np.abs(y - reg.predict(x))

    alpha = 0.05

    condCovProgram = CondConf(scoreFn, Phi_fn=phiFn, infinite_params={})
    condCovProgram.setup_problem(x_calib=XCalib, y_calib=YCalib)

    def custom_get_dual_solution(condCovProgram, S, x, alpha):
        S = np.concatenate([condCovProgram.scores_calib, [S]])
        Phi = np.vstack((condCovProgram.phi_calib, condCovProgram.Phi_fn(x)))
        zeros = np.zeros((Phi.shape[1],))

        n = len(S)
        bounds = np.column_stack((np.repeat(-alpha, n), np.repeat(1-alpha, n)))

        res = linprog(-1 * S, A_eq=Phi.T, b_eq=zeros, bounds=bounds, method='highs-ds', options={'presolve': False})
        return res.x

    def custom_verify_coverage(condCovProgram, x_test, y_test, alpha, randomize=True):
        coverages = []
        prediction_sets = []

        for x, y in zip(x_test, y_test):
            x = x.reshape(1, -1)
            true_score = condCovProgram.score_fn(x, y)

            eta = custom_get_dual_solution(condCovProgram, true_score.item(), x, alpha)[-1]
            U = np.random.uniform(-alpha, 1-alpha) if randomize else 1-alpha

            covered = (eta < U)
            coverages.append(covered)

            prediction = reg.predict(x)[0]
            S_star = condCovProgram.predict(alpha, x, lambda s, x: [prediction - s, prediction + s], exact=True, randomize=True)
            S_star = np.concatenate(S_star)
            prediction_sets.append((S_star[1], S_star[0]))  # Lower bound, Upper bound

        return np.mean(coverages), prediction_sets

    overall_coverage, overall_predictions = custom_verify_coverage(condCovProgram, XTest, YTest, alpha, randomize=True)

    group_results = {}
    group_columns = ["Group_Asian", "Group_Black", "Group_Hispanic", "Group_White"]
    for i, group in enumerate(group_columns):
        group_mask = XTest[:, 3+i] == 1
        group_coverage, group_predictions = custom_verify_coverage(condCovProgram, XTest[group_mask], YTest[group_mask], alpha, randomize=True)
        group_results[group] = (group_coverage, group_predictions)

    return overall_coverage, overall_predictions, group_results

def run_kfold_simulation(X, Y, n_splits=5):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    overall_coverages = []
    overall_predictions = []
    group_results = {group: {"coverages": [], "predictions": []} for group in ["Group_Asian", "Group_Black", "Group_Hispanic", "Group_White"]}

    for fold, (train_index, test_index) in enumerate(kf.split(X), 1):
        print(f"Running fold {fold}/{n_splits}")
        np.random.seed(fold)

        # Split data into training, calibration, and test sets
        X_train_full, X_test = X.iloc[train_index], X.iloc[test_index]
        Y_train_full, Y_test = Y.iloc[train_index], Y.iloc[test_index]

        nCalib = round(X_train_full.shape[0] * 0.3)
        calibPoints = np.random.choice(X_train_full.index, nCalib, replace=False)
        trainPoints = list(set(X_train_full.index) - set(calibPoints))

        XCalib = X.loc[calibPoints].values
        YCalib = Y.loc[calibPoints].values
        XTrain = X.loc[trainPoints].values
        YTrain = Y.loc[trainPoints].values

        # Fit the model and compute scores
        reg = LinearRegression().fit(XTrain, YTrain)
        scoreFn = lambda x, y: np.abs(y - reg.predict(x))

        alpha = 0.05

        condCovProgram = CondConf(scoreFn, Phi_fn=phiFn, infinite_params={})
        condCovProgram.setup_problem(x_calib=XCalib, y_calib=YCalib)

        overall_cov, overall_preds, group_res = run_single_simulation(X, Y, fold)
        overall_coverages.append(overall_cov)
        overall_predictions.extend(overall_preds)

        for group, (coverage, predictions) in group_res.items():
            group_results[group]["coverages"].append(coverage)
            group_results[group]["predictions"].extend(predictions)

    return overall_coverages, overall_predictions, group_results

# Run the k-fold simulations
n_splits = 10
overall_coverages, overall_predictions, group_results = run_kfold_simulation(X, Y, n_splits)

# Print aggregate results
print("\nAggregate Results:")
print(f"Overall - Mean Coverage: {np.mean(overall_coverages):.4f}, Std: {np.std(overall_coverages):.4f}")
for group, results in group_results.items():
    mean_coverage = np.mean(results["coverages"])
    std_coverage = np.std(results["coverages"])
    print(f"{group} - Mean Coverage: {mean_coverage:.4f}, Std: {std_coverage:.4f}")

# Calculate and print prediction intervals
print("\nPrediction Intervals:")
print(f"Overall - Mean: {np.mean([np.mean(p) for p in overall_predictions]):.4f}, "
      f"95% CI: [{np.percentile([p[0] for p in overall_predictions], 2.5):.4f}, "
      f"{np.percentile([p[1] for p in overall_predictions], 97.5):.4f}]")

for group, results in group_results.items():
    predictions = results["predictions"]
    mean_pred = np.mean([np.mean(p) for p in predictions])
    ci_lower = np.percentile([p[0] for p in predictions], 2.5)
    ci_upper = np.percentile([p[1] for p in predictions], 97.5)
    print(f"{group} - Mean: {mean_pred:.4f}, 95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")
