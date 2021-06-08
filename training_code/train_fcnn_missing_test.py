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
import theano.ifelse

import sys

import random

from sklearn.linear_model import LogisticRegression

os.environ['KERAS_BACKEND'] = 'theano'

from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten, Reshape
from keras.optimizers import adam, RMSprop
from keras.regularizers import l2, l1
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.layers import Dropout

from graph_convolution import GraphConv
from graph_convolution_pro import GraphConvPro

import keras.backend as K

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
def boot(y_scores, y_true, nboot=100, thresh=0.0334):

    
    auc_actual = roc_auc_score(y_true, y_scores)
    p_actual = precision_score(y_true, y_scores >= thresh)
    r_actual = recall_score(y_true, y_scores >= thresh)
    pr_actual = average_precision_score(y_true, y_scores)
    oe_actual = calc_oe(y_true, y_scores)
    bs_actual = calc_bs(y_true, y_scores)
    
    
    auc = np.zeros(nboot)
    pr = np.zeros(nboot)
    p = np.zeros(nboot)
    r = np.zeros(nboot)
    oe = np.zeros(nboot)
    bs = np.zeros(nboot)
    
    
    for i in range(nboot):
        ind = np.random.choice(range(len(y_scores)), len(y_scores))
        scores_boot = y_scores[ind]
        outcomes_boot = y_true[ind]
        auc[i] = roc_auc_score(outcomes_boot, scores_boot)
        pr[i] = average_precision_score(outcomes_boot, scores_boot)
        p[i] = precision_score(outcomes_boot, scores_boot >= thresh)
        r[i] = recall_score(outcomes_boot, scores_boot >= thresh)
        oe[i] = calc_oe(outcomes_boot, scores_boot)
        bs[i] = calc_bs(outcomes_boot, scores_boot)

    return(auc_actual, np.percentile(auc, 2.5), np.percentile(auc, 97.5),
           oe_actual, np.percentile(oe, 2.5), np.percentile(oe, 97.5),
           bs_actual, np.percentile(bs, 2.5), np.percentile(bs, 97.5),
           pr_actual, np.percentile(pr, 2.5), np.percentile(pr, 97.5),
           p_actual, np.percentile(p, 2.5), np.percentile(p, 97.5),
           r_actual, np.percentile(r, 2.5), np.percentile(r, 97.5))

def boot_cor(y_scores, pred_brca, nboot=100):
    cor_actual = pearsonr(y_scores, pred_brca)[0]
    cor = np.zeros(nboot)

    for i in range(nboot):
        ind = np.random.choice(range(len(y_scores)), len(y_scores))
        scores_boot = y_scores[ind]
        cor[i] = pearsonr(scores_boot, pred_brca[ind])[0]

    return(cor_actual, np.percentile(cor, 2.5), np.percentile(cor, 97.5))

### train FCNN

max_train = 25000
epochs = 10
seed_value = 0

if len(sys.argv)>1:
    max_train = int(sys.argv[1])
    epochs = int(sys.argv[2])
    seed_value = int(sys.argv[3])
       

nn_input = pd.read_csv('../nn_input/nn_input.csv', low_memory=False)


train_size = 800000
test_size = 87353
start_x = 2
end_x = nn_input.columns.get_loc("FDR")
ind_y = nn_input.columns.get_loc("BC.10")
ind_brcapro = nn_input.columns.get_loc("brcapro10")
pred_brca = np.array(nn_input.iloc[0:test_size, ind_brcapro])

scaler = MinMaxScaler().fit(nn_input.iloc[test_size:test_size+train_size, start_x:end_x])
scaler = MinMaxScaler().fit(nn_input.iloc[test_size:test_size+train_size, start_x:end_x])
X_train = np.array(scaler.transform(nn_input.iloc[test_size:test_size+train_size, start_x:end_x]))
y_train = np.array(nn_input.iloc[test_size:test_size+train_size, ind_y])




os.environ['PYTHONHASHSEED'] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)

fcnn = Sequential()
fcnn.add(Dense(30, input_shape=(X_train.shape[1],), activation='elu'))
fcnn.add(Dropout(0.2))
fcnn.add(Dense(10, activation='elu'))
fcnn.add(Dense(1, activation='sigmoid') )
fcnn.compile(loss='mean_squared_error', optimizer='Adam')
history = fcnn.fit(X_train[0:max_train], y_train[0:max_train], 
                   epochs=epochs,
                   batch_size=512, verbose=0)


prop_list = [0.05, 0.1, 0.3]

suffix_list = ["missingages" + str(x) for x in prop_list] + ["missingrels" + str(x) for x in prop_list] + ["missingunaff" + str(x) for x in prop_list]

for i in range(len(suffix_list)):
	suffix = suffix_list[i]
	nn_test = pd.read_csv('../nn_input/test_input_' + suffix + '.csv', low_memory=False)
	X_test = np.array(scaler.transform(nn_test.iloc[0:test_size, start_x:end_x]))
	y_test = np.array(nn_test.iloc[0:test_size, ind_y])

	pred_fcnn = fcnn.predict(X_test).flatten()
	perf = boot(pred_fcnn, y_test)

	cor = boot_cor(pred_fcnn, pred_brca)
	results ={'auc': perf[0:3],
         'oe': perf[3:6],
         'bs': perf[6:9],
         'pr_auc': perf[9:12],
         'prec': perf[12:15],
         'rec': perf[15:18],
         'corr': cor}

	filename = "../results/res_test_" + suffix + "_" + str(max_train) + "_" + str(epochs) + "_" + str(seed_value) + ".csv"
	results_df = pd.DataFrame(data=results)
	results_df.to_csv(filename)


#filename2 = "../results/pred_" + suffix + "_" + str(max_train) + "_" + str(epochs) + "_" + str(seed_value) + ".csv"


#if max_train==800000:
#	predictions = {"fcnn": pred_fcnn}
#
#	pred_df = pd.DataFrame(data=predictions)
#	pred_df.to_csv(filename2)
