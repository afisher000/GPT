# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 17:40:03 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\THz\\THzPackage')
from Pegasus import simulate_gun_to_screen4, simulate_screen4_to_und
from GPT import load_gdf, aggregate_gdf
from THz import parameters
import numpy as np
import pandas as pd

# Read Parameters
mr = parameters()

mr.nps = 100
mr.qtot = 0
mr.estimate_quads = 0

# Run simulations
phases = np.linspace(0,90,40)
gradients = np.linspace(10,20,5)
data = []
for phase in phases:
    for gradient in gradients:
        mr.linacphase = phase
        mr.linacgradient = gradient
        
        simulate_gun_to_screen4(mr, cmd_output=False)
        _, _, particle = load_gdf('beam_at_screen4.gdf', arrays_to_load=[['z', 'G']])
    
        data.append([mr.linacphase, mr.linacgradient, particle.G.mean()])
        
linac = pd.DataFrame(data, columns=['phase','gradient','G'])
#linac.to_csv('linac_scan.csv', index=False )
    
