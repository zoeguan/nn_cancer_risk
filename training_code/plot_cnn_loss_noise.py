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

import theano.ifelse

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
def boot(y_scores, y_true, nboot=100):
    
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

def boot_cor(y_scores, pred_brca, nboot=100):
    cor_actual = pearsonr(y_scores, pred_brca)[0]
    cor = np.zeros(nboot)

    for i in range(nboot):
        ind = np.random.choice(range(len(y_scores)), len(y_scores))
        scores_boot = y_scores[ind]
        cor[i] = pearsonr(scores_boot, pred_brca[ind])[0]

    return(cor_actual, np.percentile(cor, 2.5), np.percentile(cor, 97.5))


# https://www.kaggle.com/tomcwalker/keras-nn-with-custom-loss-function-for-gini-auc

import theano
# An analogue to AUC which takes the differences between each pair of true/false predictions
# and takes the average sigmoid of the differences to get a differentiable loss function.
# Based on code and ideas from https://github.com/Lasagne/Lasagne/issues/767
def soft_AUC_theano(y_true, y_pred):
    # Extract 1s
    pos_pred_vr = y_pred[y_true.nonzero()]
    # Extract zeroes
    neg_pred_vr = y_pred[theano.tensor.eq(y_true, 0).nonzero()]
    # Broadcast the subtraction to give a matrix of differences  between pairs of observations.
    pred_diffs_vr = pos_pred_vr.dimshuffle(0, 'x') - neg_pred_vr.dimshuffle('x', 0)
    # Get sigmoid of each pair.
    stats = theano.tensor.nnet.sigmoid(pred_diffs_vr * 2)
    # Take average 
    return theano.tensor.mean(stats) 

def oe(y_true, y_pred):
   return theano.tensor.sum(y_true)/theano.tensor.sum(y_pred)

### train CNN

cnn_input = pd.read_csv('../nn_input/cnn_input.csv', low_memory=False)
graph_mat = pd.read_csv('../nn_input/graph_mat.csv')
graph_mat = np.array(graph_mat)
num_neighbors = graph_mat.shape[1]

train_size = 800000
test_size = 87353
start_x = 2
end_x = 198-14
ind_y = 210-14
ind_brcapro = 208-14


features = 7


y_train = np.array(cnn_input.iloc[test_size:test_size+train_size, ind_y])
y_test = np.array(cnn_input.iloc[0:test_size, ind_y])

scaler2 = MinMaxScaler().fit(cnn_input.iloc[test_size:test_size+train_size, start_x:(end_x+features)])
X_train_cnn = np.array(scaler2.transform(cnn_input.iloc[test_size:test_size+train_size, start_x:(end_x+features)]))
X_test_cnn = np.array(scaler2.transform(cnn_input.iloc[0:test_size, start_x:(end_x+features)]))
num_neighbors = graph_mat.shape[1]
X_train_2 = X_train_cnn.reshape((X_train_cnn.shape[0], 27, features))
X_test_2 = X_test_cnn.reshape((X_test_cnn.shape[0], 27, features))

pred_brca = np.array(cnn_input.iloc[0:test_size, ind_brcapro])



max_train = 25000
epochs = 10
seed_value = 0

if len(sys.argv)>1:
    max_train = int(sys.argv[1])
    epochs = int(sys.argv[2])
    seed_value = int(sys.argv[3])


os.environ['PYTHONHASHSEED'] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)


g_model = Sequential()
g_model.add(GraphConv(filters=10, neighbors_ix_mat = graph_mat, num_neighbors=num_neighbors, activation='elu', input_shape=(X_train_2.shape[1], features)))
g_model.add(Dropout(0.2))
g_model.add(GraphConv(filters=5, neighbors_ix_mat = graph_mat, num_neighbors=num_neighbors, activation='elu'))
g_model.add(GraphConvPro())
g_model.add(Flatten())
g_model.add(Dense(1, activation='sigmoid'))
g_model.compile(loss='mean_squared_error', optimizer='Adam', metrics=[soft_AUC_theano,oe])

X_train_3 = X_train_2[0:max_train]
history = g_model.fit(X_train_3.reshape(X_train_3.shape[0], X_train_3.shape[1], features),
                          y_train[0:max_train],
                          epochs=epochs, 
                          batch_size=512, verbose=0, validation_split=0.1)



import matplotlib.pyplot as plt



figname = "../tuning_results/val_loss_cnn_noise" + str(max_train) + "_" + str(epochs)
fig = plt.figure()
plt.plot(history.history['val_loss'])
ymin = min(history.history['val_loss']) - 0.00005
ymax = min(history.history['val_loss']) + 0.0005
plt.ylim(ymin, ymax)
fig.savefig(figname)
plt.close(fig)

figname = "../tuning_results/val_auc_cnn_noise_" + str(max_train) + "_" + str(epochs)
fig2 = plt.figure()
plt.plot(history.history['val_soft_AUC_theano'])
fig2.savefig(figname)
plt.close(fig2)


figname = "../tuning_results/val_oe_cnn_noise_" + str(max_train) + "_" + str(epochs)
fig3 = plt.figure()
plt.plot(history.history['val_oe'])
fig3.savefig(figname)
plt.close(fig3)