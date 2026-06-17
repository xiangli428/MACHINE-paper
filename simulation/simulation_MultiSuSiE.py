#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 18 13:10:50 2024

@author: r9user9
"""

import sys
import numpy as np
import pandas as pd
import pickle
import MultiSuSiE
import numba

# numba.config.NUMBA_NUM_THREADS = 1
# numba.config.NUMBA_DEFAULT_NUM_THREADS = 1



Dir = sys.argv[1]
LD_f1 = sys.argv[2]
LD_f2 = sys.argv[3]
N = [int(sys.argv[4]), int(sys.argv[5])]

data = pd.read_csv(Dir + "data.txt", sep='\t', header = 0)
R_list = [np.loadtxt(LD_f1, delimiter = '\t'), 
          np.loadtxt(LD_f2, delimiter = '\t')]

ss_fit = MultiSuSiE.multisusie_rss(
    z_list = [data.iloc[:,1], data.iloc[:,2]],
    R_list = R_list,
    rho = np.array([[1, 0.8], [0.8, 1]]),
    population_sizes = N,
    L = 5,
    scaled_prior_variance = 0.2,
    min_abs_corr = 0.5,
    low_memory_mode = False,
    verbose = True)

with open(Dir + "ss_fit.pkl", 'wb') as p:
    pickle.dump(ss_fit, p)

data.loc[:,"PIP"] = ss_fit.pip

data.loc[:,"CS"] = 0

for l in range(0,5):
    if ss_fit.sets[3][l]:
        data.loc[ss_fit.sets[0][l],"CS"] = l + 1

data.to_csv(Dir + "data.txt", sep = '\t', 
            index = False, header = True)
