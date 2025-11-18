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
import time
import MultiSuSiE

N1 = 2e5
N2_seq = [2e4,2e5]

Dir = sys.argv[1]

data_full = pd.read_csv(Dir + "data.txt", sep='\t', header = 0)
R_list = [np.loadtxt(Dir + "LD_1.txt", delimiter = '\t'), 
          np.loadtxt(Dir + "LD_2.txt", delimiter = '\t')]

runtime = np.zeros(2)

for m in [0,1]:
    data = data_full.iloc[:,[0,1,m+2]].copy()

    tic = time.time()
    ss_fit = MultiSuSiE.multisusie_rss(
        z_list = [data.iloc[:,1], data.iloc[:,2]],
        R_list = R_list,
        rho = np.array([[1, 0.8], [0.8, 1]]),
        population_sizes = [N1,N2_seq[m]],
        L = 5,
        scaled_prior_variance = 0.2,
        min_abs_corr = 0.5,
        low_memory_mode = False)
    toc = time.time()
    runtime[m] = toc - tic
    
    with open(Dir + "ss_fit_N2-%d.pkl" %N2_seq[m], 'wb') as p:
        pickle.dump(ss_fit, p)
    
    data.loc[:,"PIP"] = ss_fit.pip
    
    data.loc[:,"CS"] = 0
    
    for l in range(0,5):
        if ss_fit.sets[3][l]:
            data.loc[ss_fit.sets[0][l],"CS"] = l + 1
    
    data.to_csv(Dir + "data_N2-%d.txt" %N2_seq[m], sep = '\t', 
                index = False, header = True)

np.savetxt(Dir + "runtime.txt", runtime, delimiter = '\t')
