# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:13:37 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')

import subprocess
import scipy.constants as sc
import pandas as pd
from GPT import write_MRfile, load_gdf
from Pegasus import estimate_quads

flags = {
    'plot_traj':1,
    'estimate_quads':0,
    'estimate_R56':0,
    'beam_flag':0,
    'xband_tf':0,
    'xband_mode':0
    }

params = {
    'R':1.63e-3,
    'f0':.315e12,
    'lamu':.032,
    'ku':2*sc.pi/.032,
    'nemit':2e-6,
    'nps':100,
    'qtot':-125e-12
    }

fields = {
    'B0':.730,
    'gunloopmv':57, 
    'gunphasedeg':6*2.856,
    'linacgradient':20,
    'linacphase':0,
    'sol1':1.05,
    'quads':[3.4, -6.4, 3.4],
    'sol2':1.5,
    'chicane_field':.2
    }

laser = {
    'ztime':2e-12,
    'laser_spotsize':1e-3
    }

undulator = {
    'taperdelay':.06,
    'taper1':0,
    'taper2':-0.35
    }


mr = pd.Series({**flags, **params, **fields, **laser, **undulator})
mr['K'] = sc.elementary_charge*mr.B0/sc.m_e/sc.c/mr.ku

write_MRfile('pegasus.mr', mr)

gpt_call = 'mr -o BeamatScreen4.gdf pegasus.mr gpt Gun_to_Screen4.in'
cmd = subprocess.run(gpt_call.split(' '), capture_output=True, encoding='UTF-8')
print(cmd.stderr)

consts, params, particle = load_gdf('BeamatScreen4.gdf', [['x','Bx','y','By','z','t','G']])
G = particle.G[particle.t.notnull()].mean()
s_beam = particle[['x','Bx','y','By']][particle.t.notnull()].cov().values
quads = estimate_quads(s_beam, G, 5.16-3.191)



