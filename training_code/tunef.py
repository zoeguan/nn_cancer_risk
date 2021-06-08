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
import sys
import random
import itertools

import pandas as pd
from sklearn.metrics import roc_auc_score, precision_score, recall_score, average_precision_score
from scipy.stats.stats import pearsonr  
from scipy.stats.stats import spearmanr
import theano.ifelse

os.environ['KERAS_BACKEND'] = 'theano'

from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten, Reshape
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.layers import Dropout
from keras.optimizers import Adam
from keras import regularizers

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



def create_model(learning_rate, weight_decay, dropout, hidden_layers, activation): 
    
    model = Sequential()
    
    for i in range(len(hidden_layers)): 
        if (i==0):
            model.add(Dense(hidden_layers[i], input_shape = (input_size,), activation=activation, kernel_regularizer=regularizers.l2(weight_decay)))
            model.add(Dropout(dropout))
        else:
            model.add(Dense(hidden_layers[i], activation=activation, kernel_regularizer=regularizers.l2(weight_decay)))
      
        
    model.add(Dense(1, activation='sigmoid'))
    optimizer = Adam(lr=learning_rate)
    model.compile(optimizer=optimizer, loss='mean_squared_error')
    return model



### command line arguments

max_train = 500000
seed_value = 0
niter = 20


if len(sys.argv)>1:
    max_train = int(sys.argv[1])
    seed_value = int(sys.argv[2])
    niter = int(sys.argv[3])

### load inputs 

nn_input = pd.read_csv('../nn_input/nn_input.csv', low_memory=False)

train_size = 800000
test_size = 87353
start_x = 2
end_x = nn_input.columns.get_loc("FDR")
ind_y = nn_input.columns.get_loc("BC.10")
ind_brcapro = nn_input.columns.get_loc("brcapro10")
pred_brca = np.array(nn_input.iloc[0:test_size, ind_brcapro])

scaler = MinMaxScaler().fit(nn_input.iloc[test_size:test_size+train_size, start_x:end_x])
X_train = np.array(scaler.transform(nn_input.iloc[test_size:test_size+train_size, start_x:end_x]))
y_train = np.array(nn_input.iloc[test_size:test_size+train_size, ind_y])

X_test = np.array(scaler.transform(nn_input.iloc[0:test_size, start_x:end_x]))
y_test = np.array(nn_input.iloc[0:test_size, ind_y])

val_size = int(0.1*max_train)
X_val = X_train[(max_train-val_size):max_train]
y_val = y_train[(max_train-val_size):max_train]

input_size = X_train.shape[1]


### search space 

learning_rate_list = [1e-4, 1e-3, 1e-2]

weight_decay_list = [0, 1e-4, 1e-3, 1e-2]

activation_list = ['elu', 'relu']

dropout_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5]

nodes_list = [i for i in range(10, 110, 10)]
hidden_layers_list = list(itertools.product(nodes_list)) + list(itertools.product(nodes_list, nodes_list))  + list(itertools.product(nodes_list[0:5], nodes_list[0:5], nodes_list[0:5]))


hyperparam_list = list(itertools.product(learning_rate_list, weight_decay_list, dropout_list, hidden_layers_list, activation_list))

random.seed(seed_value)
hyp_combs = random.sample(hyperparam_list, niter)
print(hyp_combs)


### tuning


perf_val_list = [0]*len(hyp_combs)

os.environ['PYTHONHASHSEED'] = str(0)
random.seed(0)
np.random.seed(0)

for i in range(len(hyp_combs)):
    print(i)
    hyps = hyp_combs[i]
    fcnn = create_model(hyps[0], hyps[1], hyps[2], hyps[3], hyps[4])
    history = fcnn.fit(X_train[0:(max_train-val_size)], y_train[0:(max_train-val_size)], 
                   epochs=20,
                   batch_size=512, verbose=0)
    pred_fcnn_val = fcnn.predict(X_val).flatten()
    perf_val = boot(pred_fcnn_val, y_val, 1)
    perf_val_list[i] = [perf_val[x] for x in [0, 3, 6, 9]]

results ={'learning_rate': [hyp_combs[x][0] for x in range(len(hyp_combs))],
          'decay': [hyp_combs[x][1] for x in range(len(hyp_combs))],
          'dropout': [hyp_combs[x][2] for x in range(len(hyp_combs))],
          'hidden': [hyp_combs[x][3] for x in range(len(hyp_combs))],
          'act': [hyp_combs[x][4] for x in range(len(hyp_combs))],
          'auc_val': [perf_val_list[x][0] for x in range(len(hyp_combs))],
          'oe_val': [perf_val_list[x][1] for x in range(len(hyp_combs))],
          'bs_val': [perf_val_list[x][2] for x in range(len(hyp_combs))],
          'prauc_val': [perf_val_list[x][3] for x in range(len(hyp_combs))]}

results_df = pd.DataFrame(data=results)
results_df.to_csv("../tuning_results/tunef_" + str(max_train) + "_" + str(seed_value) + "_" + str(niter) + ".csv")



