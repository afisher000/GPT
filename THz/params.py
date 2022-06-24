# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 12:23:53 2022

@author: afisher
"""


import scipy.constants as sc
import pandas as pd

def read_parameters():
    flags = {
        'plot_traj':1,
        'estimate_quads':1,
        'estimate_R56':0,
        'beam_flag':0,
        'xband_tf':0,
        'xband_mode':0,
        'chicane_flag':0
        }
    
    params = {
        'R':1.63e-3,
        'f0':0.315e12,
        'lamu':.032,
        'ku':2*sc.pi/.032,
        'nemit':2e-6,
        'nps':10,
        'qtot':-125e-12
        }
    
    fields = {
        'B0':0.730,
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
    return mr