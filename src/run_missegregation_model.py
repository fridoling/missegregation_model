#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pickle
import os
import sys
from datetime import datetime
import multiprocessing as mp
import itertools as it


if len(sys.argv)>1:
    res_folder = "./res/"+sys.argv[4]+'/'
else:
    res_folder = "./res/"+str(datetime.now()).split('.')[0]+'/'
if not os.path.isdir('res'):
    os.mkdir('res')
if not os.path.isdir(res_folder):
    os.mkdir(res_folder)

params = {}
params['n_gen'] = 10
params['pop_size'] = 100
params['m_vec'] = np.arange(0.01,0.1,0.01)
params['frac4_vec'] = np.arange(0.1,1.0,0.1)
params['ms_factor'] = 2
params['fertility'] = 1
params['fert_factor'] = 0.9


pool = mp.Pool()
m_range = range(len(params['m_vec']))
f_range = range(len(params['frac4_vec']))
res = pool.starmap(run_simulation_parallel, [(i, j, params) for i,j in it.product(m_range, f_range)])
pool.close()
trajs = np.stack(results)

with open(res_folder+'simulation_data.pickle', 'wb') as f:
    pickle.dump((trajs, params), f)