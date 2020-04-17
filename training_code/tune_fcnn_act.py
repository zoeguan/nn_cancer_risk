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

def mse(y_true, y_pred):
    return K.mean((y_true-y_pred)**2)

### helper functions

# function for calculating calibration
def calc_oe(y_true, y_scores):
    return(np.sum(y_true)/np.sum(y_scores))

# function for calculating Brier Score
def calc_bs(y_true, y_scores):
    return(np.mean((y_true-y_scores)**2))

### tune architecture

nn_input = pd.read_csv('../nn_input/nn_input.csv', low_memory=False)

train_size = 500000
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

max_train = 450000


seed_value = 0
os.environ['PYTHONHASHSEED'] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)


act_fns = ['sigmoid', 'relu', 'elu']

bs = list()
auc = list()

for i in range(len(act_fns)):
    print(i)
    fcnn = Sequential()
    fcnn.add(Dense(30, input_shape=(X_train.shape[1],), activation=act_fns[i]))
    fcnn.add(Dense(10, activation=act_fns[i]) )
    fcnn.add(Dense(1, activation='sigmoid') )
    fcnn.compile(loss='mean_squared_error', optimizer='Adam')
    history = fcnn.fit(X_train[0:max_train], y_train[0:max_train], 
                   epochs=20,
                   batch_size=512, verbose=0)
    pred_fcnn = fcnn.predict(X_train[max_train:train_size]).flatten()
    auc.append(roc_auc_score(y_train[max_train:train_size], pred_fcnn))
    bs.append(calc_bs(y_train[max_train:train_size], pred_fcnn))


perf = {"activation": act_fns,
        "auc": auc,
        "bs": bs}

perf_df = pd.DataFrame(data=perf)
perf_df.to_csv('../tuning_results/tune_fcnn_act.csv')
