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
from graph_convolution import GraphConv
from graph_convolution_pro import GraphConvPro
from keras.layers import Dropout

import theano.ifelse

from keras.callbacks import EarlyStopping, ModelCheckpoint

import keras.backend as K

def mse(y_true, y_pred):
    return K.mean((y_true-y_pred)**2)

# function for calculating calibration
def calc_oe(y_true, y_scores):
    return(np.sum(y_true)/np.sum(y_scores))

# function for calculating Brier Score
def calc_bs(y_true, y_scores):
    return(np.mean((y_true-y_scores)**2))


cnn_input = pd.read_csv('../nn_input/cnn_input.csv', low_memory=False)
graph_mat = pd.read_csv('../nn_input/graph_mat.csv')
graph_mat = np.array(graph_mat)
num_neighbors = graph_mat.shape[1]




train_size = 500000
test_size = 87353
start_x = 2
end_x = 198-14
ind_y = 210-14
features = 7

y_train = np.array(cnn_input.iloc[test_size:test_size+train_size, ind_y])
y_test = np.array(cnn_input.iloc[0:test_size, ind_y])

scaler2 = MinMaxScaler().fit(cnn_input.iloc[test_size:test_size+train_size, start_x:(end_x+features)])
X_train_cnn = np.array(scaler2.transform(cnn_input.iloc[test_size:test_size+train_size, start_x:(end_x+features)]))
X_test_cnn = np.array(scaler2.transform(cnn_input.iloc[0:test_size, start_x:(end_x+features)]))
num_neighbors = graph_mat.shape[1]
X_train_2 = X_train_cnn.reshape((X_train_cnn.shape[0], 27, features))
X_test_2 = X_test_cnn.reshape((X_test_cnn.shape[0], 27, features))

max_train = 450000


sizes = [3, 5, 10]
bs = list()
auc = list()

seed_value = 0
os.environ['PYTHONHASHSEED'] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)

for s in sizes:
    print(s)
    g_model = Sequential()
    g_model.add(GraphConv(filters=s, neighbors_ix_mat = graph_mat, num_neighbors=num_neighbors, activation='elu', input_shape=(X_train_2.shape[1], features)))
    g_model.add(Dropout(0.2))
    g_model.add(GraphConv(filters=s, neighbors_ix_mat = graph_mat, num_neighbors=num_neighbors, activation='elu'))
    g_model.add(GraphConvPro())
    g_model.add(Flatten())
    g_model.add(Dense(1, activation='sigmoid'))
    g_model.compile(loss='mean_squared_error', optimizer='Adam')
    X_train_3 = X_train_2[0:max_train]
    history = g_model.fit(X_train_3.reshape(X_train_3.shape[0], X_train_3.shape[1], features),
                          y_train[0:max_train],
                          epochs=10, 
                          batch_size=512, verbose=0)
    X_train_4 = X_train_2[max_train:train_size]
    y_pred = g_model.predict(X_train_4.reshape(X_train_4.shape[0],X_train_4.shape[1],features)).flatten()

    auc.append(roc_auc_score(y_train[max_train:train_size], y_pred))
    bs.append(calc_bs(y_train[max_train:train_size], y_pred))


perf = {"size": sizes, 
        "auc": auc,
        "bs": bs}

perf_df = pd.DataFrame(data=perf)
perf_df.to_csv('../tuning_results/tune_cnn1.csv')
