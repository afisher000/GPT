# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 17:40:03 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\THz\\THzPackage')
import Pegasus as PEG
import GPT
import THz
import numpy as np
import pandas as pd
import os

# Read Parameters
mr = THz.parameters()

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
        
        PEG.simulate_gun_to_screen4(mr, cmd_output=False)
        _, _, particle = GPT.load_gdf('beam_at_screen4.gdf', arrays_to_load=[['z', 'G']])
    
        data.append([mr.linacphase, mr.linacgradient, particle.G.mean()])
        
linac = pd.DataFrame(data, columns=['phase','gradient','G'])
# linac.to_csv(os.path.join('Parameter Scans','linac_scan.csv'), index=False )
    
