from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.exceptions import NotFittedError
from datetime import datetime
from sklearn.preprocessing import MinMaxScaler

import os
import urllib

import pandas as pd
from sklearn.metrics import roc_auc_score, precision_score, recall_score, average_precision_score
from scipy.stats.stats import pearsonr  
from scipy.stats.stats import spearmanr

import sys

import random

from sklearn.linear_model import LogisticRegression


def mse(y_true, y_pred):
    return K.mean((y_true-y_pred)**2)

### helper functions

# function for calculating calibration
def calc_oe(y_true, y_scores):
    return(np.sum(y_true)/np.sum(y_scores))

# function for calculating Brier Score
def calc_bs(y_true, y_scores):
    return(np.sqrt(np.mean((y_true-y_scores)**2)))

# function for calculating bootstrap CIs for AUC, calibration, Brier Score
def boot(y_scores, y_true, nboot=100):

    
    auc_actual = roc_auc_score(y_true, y_scores)
    pr_actual = average_precision_score(y_true, y_scores)
    oe_actual = calc_oe(y_true, y_scores)
    bs_actual = calc_bs(y_true, y_scores)
    
    
    auc = np.zeros(nboot)
    pr = np.zeros(nboot)
    oe = np.zeros(nboot)
    bs = np.zeros(nboot)
    
    
    for i in range(nboot):
        ind = np.random.choice(range(len(y_scores)), len(y_scores))
        scores_boot = y_scores[ind]
        outcomes_boot = y_true[ind]
        auc[i] = roc_auc_score(outcomes_boot, scores_boot)
        pr[i] = average_precision_score(outcomes_boot, scores_boot)
        oe[i] = calc_oe(outcomes_boot, scores_boot)
        bs[i] = calc_bs(outcomes_boot, scores_boot)

    return(auc_actual, np.percentile(auc, 2.5), np.percentile(auc, 97.5),
           oe_actual, np.percentile(oe, 2.5), np.percentile(oe, 97.5),
           bs_actual, np.percentile(bs, 2.5), np.percentile(bs, 97.5),
           pr_actual, np.percentile(pr, 2.5), np.percentile(pr, 97.5))

def boot_cor(y_scores, pred_brca, nboot=100):
    cor_actual = pearsonr(y_scores, pred_brca)[0]
    cor = np.zeros(nboot)

    for i in range(nboot):
        ind = np.random.choice(range(len(y_scores)), len(y_scores))
        scores_boot = y_scores[ind]
        cor[i] = pearsonr(scores_boot, pred_brca[ind])[0]

    return(cor_actual, np.percentile(cor, 2.5), np.percentile(cor, 97.5))


pred_all = pd.read_csv('../results/pred_lr_noise.csv')
pred_brca = pred_all.iloc[:, 1]
pred_fcnn = pd.read_csv('../results/pred_noise80000_50_0.csv').iloc[:, 1]
pred_cnn =  pd.read_csv('../results/pred_cnn_noise800000_15_0.csv').iloc[:, 1]
pred_lr = pred_all.iloc[:, 2]

nn_input = pd.read_csv('../nn_input/nn_input_percent10.csv', low_memory=False)
test_size = 87353
ind_y = nn_input.columns.get_loc("BC.10")
y_test = np.array(nn_input.iloc[0:test_size, ind_y])



nboot = 1000
oe = np.zeros((4, nboot))
auc = np.zeros((4, nboot))
bs = np.zeros((4, nboot))
corr = np.zeros((4, nboot))
prauc = np.zeros((4, nboot))




np.random.seed(0)
for i in range(nboot):
    print(i)
    ind = np.random.choice(range(len(y_test)), len(y_test))

    if (i==0):
        ind = np.arange(0, len(y_test), 1) 

    outcomes_boot = y_test[ind]

    auc[0, i] = roc_auc_score(outcomes_boot, pred_fcnn[ind])
    oe[0, i] = calc_oe(outcomes_boot, pred_fcnn[ind])
    bs[0, i] = calc_bs(outcomes_boot, pred_fcnn[ind])
    corr[0, i] = pearsonr(pred_fcnn[ind], pred_brca[ind])[0]
    prauc[0, i] = average_precision_score(outcomes_boot, pred_fcnn[ind])

    auc[1, i] = roc_auc_score(outcomes_boot, pred_cnn[ind])
    oe[1, i] = calc_oe(outcomes_boot, pred_cnn[ind])
    bs[1, i] = calc_bs(outcomes_boot, pred_cnn[ind])
    corr[1, i] = pearsonr(pred_cnn[ind], pred_brca[ind])[0]
    prauc[1, i] = average_precision_score(outcomes_boot, pred_cnn[ind])

    auc[2, i] = roc_auc_score(outcomes_boot, pred_brca[ind])
    oe[2, i] = calc_oe(outcomes_boot, pred_brca[ind])
    bs[2, i] = calc_bs(outcomes_boot, pred_brca[ind])
    corr[2, i] = 1
    prauc[2, i] = average_precision_score(outcomes_boot, pred_brca[ind])

    auc[3, i] = roc_auc_score(outcomes_boot, pred_lr[ind])
    oe[3, i] = calc_oe(outcomes_boot, pred_lr[ind])
    bs[3, i] = calc_bs(outcomes_boot, pred_lr[ind])
    corr[3, i] = pearsonr(pred_lr[ind], pred_brca[ind])[0]
    prauc[3, i] = average_precision_score(outcomes_boot, pred_lr[ind])


auc_df = pd.DataFrame(auc)
bs_df = pd.DataFrame(bs)
oe_df = pd.DataFrame(oe)
corr_df = pd.DataFrame(corr)
prauc_df = pd.DataFrame(prauc)

auc_df.to_csv('../results/auc_noise_boot.csv')
oe_df.to_csv('../results/oe_noise_boot.csv')
bs_df.to_csv('../results/bs_noise_boot.csv')
corr_df.to_csv('../results/corr_noise_boot.csv')
prauc_df.to_csv('../results/prauc_noise_boot.csv')



