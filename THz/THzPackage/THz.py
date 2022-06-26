# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 12:23:53 2022

@author: afisher
"""

import scipy.constants as sc
import pandas as pd
from scipy.interpolate import interp2d

def estimate_linacphase(G, linacgradient):
    ''' Estimate linacphase from linacgradient and desired energy using
    data from linac scan simulations. '''
    
    file = 'C:\\Users\\afisher\\Documents\\GitHub\\GPT\\THz\\linac_scan.csv'
    linac = pd.read_csv(file)
    f = interp2d(linac.G, linac.gradient, linac.phase)
    return f(G, linacgradient)[0]
    
def parameters():
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
        'nps':100,
        'qtot':-125e-12
        }
    
    fields = {
        'B0':0.730,
        'gunloopmv':57, 
        'gunphasedeg':6*2.856,
        'linacgradient':20,
        'G':15,
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
        'nperiods':30,
        'taperdelay':.06,
        'taper1':0,
        'taper2':-0.35
        }
    
    fel = {
        'nfreq':61,
          
        }
    
    beamline = {
        'screen4_pos':3.191,
        'und_pos':5.16,
        'sol2_pos':4.764
        }
    
    mr = pd.Series({**flags, **params, **fields, **laser, **undulator, **fel, **beamline})
    mr['K'] = sc.elementary_charge*mr.B0/sc.m_e/sc.c/mr.ku
    mr['linacphase'] = estimate_linacphase(mr.G, mr.linacgradient)
    return mr