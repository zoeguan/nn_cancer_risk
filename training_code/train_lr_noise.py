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
from sklearn.metrics import roc_auc_score
from scipy.stats.stats import pearsonr  
from scipy.stats.stats import spearmanr

import sys

import random

from sklearn.linear_model import LogisticRegression

os.environ['KERAS_BACKEND'] = 'theano'

from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten, Reshape
from keras.optimizers import adam, RMSprop
from keras.regularizers import l2, l1
from keras.callbacks import EarlyStopping, ModelCheckpoint

from graph_convolution import GraphConv
from graph_convolution_pro import GraphConvPro

import keras.backend as K


### helper functions

# function for calculating calibration
def calc_oe(y_true, y_scores):
    return(np.sum(y_true)/np.sum(y_scores))

# function for calculating Brier Score
def calc_bs(y_true, y_scores):
    return(np.sqrt(np.mean((y_true-y_scores)**2)))

# function for calculating bootstrap CIs for AUC, calibration, Brier Score
def boot(y_scores, y_true, nboot=200):
    
    auc_actual = roc_auc_score(y_true, y_scores)
    oe_actual = calc_oe(y_true, y_scores)
    bs_actual = calc_bs(y_true, y_scores)
    
    
    auc = np.zeros(nboot)
    oe = np.zeros(nboot)
    bs = np.zeros(nboot)
    
    
    for i in range(nboot):
        ind = np.random.choice(range(len(y_scores)), len(y_scores))
        scores_boot = y_scores[ind]
        outcomes_boot = y_true[ind]
        auc[i] = roc_auc_score(outcomes_boot, scores_boot)
        oe[i] = calc_oe(outcomes_boot, scores_boot)
        bs[i] = calc_bs(outcomes_boot, scores_boot)

    return(auc_actual, np.percentile(auc, 2.5), np.percentile(auc, 97.5),
           oe_actual, np.percentile(oe, 2.5), np.percentile(oe, 97.5),
           bs_actual, np.percentile(bs, 2.5), np.percentile(bs, 97.5))

def boot_cor(y_scores, pred_brca, nboot=200):
    cor_actual = pearsonr(y_scores, pred_brca)[0]
    cor = np.zeros(nboot)

    for i in range(nboot):
        ind = np.random.choice(range(len(y_scores)), len(y_scores))
        scores_boot = y_scores[ind]
        cor[i] = pearsonr(scores_boot, pred_brca[ind])[0]

    return(cor_actual, np.percentile(cor, 2.5), np.percentile(cor, 97.5))

### train FCNN

nn_input = pd.read_csv('../nn_input/nn_input_noise.csv', low_memory=False)

train_size = 800000
test_size = 87353
start_x = 2
end_x = 198-14
ind_y = 202+1-14
ind_brcapro = 200+1-14
pred_brca = np.array(nn_input.iloc[0:test_size, ind_brcapro])

scaler = MinMaxScaler().fit(nn_input.iloc[test_size:test_size+train_size, start_x:end_x])
scaler = MinMaxScaler().fit(nn_input.iloc[test_size:test_size+train_size, start_x:end_x])
X_train = np.array(scaler.transform(nn_input.iloc[test_size:test_size+train_size, start_x:end_x]))
y_train = np.array(nn_input.iloc[test_size:test_size+train_size, ind_y])

X_test = np.array(scaler.transform(nn_input.iloc[0:test_size, start_x:end_x]))
y_test = np.array(nn_input.iloc[0:test_size, ind_y])


LR = LogisticRegression(C=100000)
LR.fit(X_train, y_train)

predictions_LR_list = list(LR.predict_proba(X_test))
pred_lr = np.array([p[1] for p in predictions_LR_list])


perf = boot(pred_lr, y_test)





results ={'auc': perf[0:3],
         'oe': perf[3:6],
         'bs': perf[6:9]}

filename = "../results/res_lr_noise.csv"

results_df = pd.DataFrame(data=results)
results_df.to_csv(filename)

predictions = {"lr": pred_lr,
              "brcapro": pred_brca}

pred_df = pd.DataFrame(data=predictions)
pred_df.to_csv('../results/pred_lr_noise.csv')



