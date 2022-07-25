# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 12:23:53 2022

@author: afisher
"""

import scipy.constants as sc
import pandas as pd
from scipy.interpolate import griddata
import os
import numpy as np

# Separate quads?
# Add quad positions to parameters?


def parameters():
    flags = dict(estimate_quads=0,
                 estimate_chicane=1,
                 estimate_sol2=1,
                 beam_flag=0,
                 xband_tf=0,
                 xband_mode=0,
                 chicane_flag=2,
                 fel_flag=0,
                 GPT_warnings=1)
    
    params = dict(R=1.63e-3,
                  f0=.315e12,
                  lamu=.032,
                  ku=2*sc.pi/.032,
                  nemit=2e-6,
                  nps=1000,
                  qtot=-125e-12)
    
    fields = dict(B0=.730,
                  gunloopmv=57,
                  gunphasedeg=6*2.856,
                  linacgradient=20,
                  G=15,
                  sol1=1.05,
                  quads=[3.4, -6.4, 3.4],
                  sol2=0.8,
                  chicane_field=.2)

    chicane = dict(chicane_dx=0.005,
                   chicane_offset=.02,
                   chicane_D = .03,
                   chicane_L = .05)
    
    laser = dict(ztime=2e-12,
                 laser_spotsize=1e-3)

    undulator = dict(nperiods=30,
                     taperdelay=.06,
                     taper1=0,
                     taper2=-.35,
                     undpass=0)
    
    fel = dict(nfreq=121)
    
    beamline = dict(screen4_pos=3.191,
                    und_pos=5.16,
                    sol2_pos=4.764,
                    chicane_pos=3.984,
                    boxend_pos=4.327)
    
    
    mr = pd.Series({**flags, **params, **fields, **chicane, **laser, **undulator, **fel, **beamline})
    mr['K'] = sc.elementary_charge*mr.B0/sc.m_e/sc.c/mr.ku
    mr['kp'] = 1.8412/mr.R
    mr['kz'] = np.sqrt((2*sc.pi*mr.f0/sc.c)**2 - mr.kp**2)
    mr['end_pos'] = mr.und_pos + 1.032
    mr['vp'] = 2*sc.pi*mr.f0/mr.kz
    mr['vg'] = sc.c**2/mr.vp
    
    linac_scan = pd.read_csv(os.path.join('Parameter Scans','linac_scan.csv'))
    mr['linacphase'] = griddata((linac_scan.G, linac_scan.gradient), 
                                linac_scan.phase, 
                                (mr.G, mr.linacgradient)).tolist()

    return mr