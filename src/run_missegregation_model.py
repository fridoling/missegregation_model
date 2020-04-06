#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pickle
import os
import sys
from datetime import datetime
import multiprocessing as mp
import itertools as it

from missegregation_model import *

if len(sys.argv)>1:
    res_folder = "./res/"+sys.argv[1]+'/'
else:
    res_folder = "./res/"+str(datetime.now()).split('.')[0]+'/'
if not os.path.isdir('res'):
    os.mkdir('res')
if not os.path.isdir(res_folder):
    os.mkdir(res_folder)

params = {}
params['n_gen'] = 10
params['pop_size'] = 100
params['m_vec'] = np.linspace(0.0,0.5,11)
params['frac4_vec'] = np.linspace(0.0,1.0,11)
params['ms_factor'] = 2
params['fertility'] = 1
params['fert_factor'] = 0.9
params['ap_loss'] = 0.1
params['ms8'] = 0.001
params['ap8_gain'] = 0


pool = mp.Pool()
m_range = range(len(params['m_vec']))
f_range = range(len(params['frac4_vec']))
sys.stdout.write('running simulation for '+str(len(params['m_vec']))+'x'+str(len(params['frac4_vec']))+' index pairs...\n')
res = pool.starmap(run_simulation, [(i, j, params) for i,j in it.product(m_range, f_range)])
pool.close()

trajs = np.reshape(res, (len(params['m_vec']), len(params['frac4_vec']), params['n_gen'], 2))
with open(res_folder+'simulation_data.pickle', 'wb') as f:
    pickle.dump((trajs, params), f)
