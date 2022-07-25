# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 09:36:24 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\THz\\THzPackage')

import Pegasus as PEG
import GPT
import THz
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import subprocess
import os


# Read Parameters
# =============================================================================
plt.close('all')
mr = THz.parameters()


# Run simulations
# =============================================================================
mr.fel_flag = 1
mr.GPT_warnings = 0

mr.taper2 = -0.35
mr.taper1 = -.2

PEG.simulate(mr, 'boxend','und')
GPT.edit_screen_gdf('beam_at_und.gdf', target_G=13.8)
PEG.simulate(mr, 'und', 'end', undpass=0)
# PEG.simulate(mr, 'und', 'end', undpass=1)
# PEG.simulate(mr, 'und', 'end', undpass=2)


# Read in data, drop non-interesting timesteps
consts, params, part0, fel0, tfel0 = GPT.load_gdf('und_pass0.gdf', arrays_to_load=[['z','G','t','Bz'],
                                                                                  ['freq00','A00','phi00','lam00'],
                                                                                  ['t00','At00']])
_, _, part1, fel1, tfel1 = GPT.load_gdf('und_pass1.gdf', arrays_to_load=[['z','G','t','Bz'],
                                                                                  ['freq00','A00','phi00','lam00'],
                                                                                  ['t00','At00']])
part = part0[part0.t.isna()]
timesteps = part.index.get_level_values(level=0).unique()[:-1]
part0 = part0.loc[timesteps]
part1 = part1.loc[timesteps]
fel0 = fel0.loc[timesteps]
fel1 = fel1.loc[timesteps]
tfel0 = tfel0.loc[timesteps]
tfel1 = tfel1.loc[timesteps]

# Plot Final Pulses
fig, ax = plt.subplots()
PULSE0 = tfel0.loc[timesteps[-1]]
PULSE1 = tfel1.loc[timesteps[-1]]
ax.plot(PULSE0.t00, PULSE0.At00, label='Pass 0')
ax.plot(PULSE1.t00, PULSE1.At00, label='Pass 1')
ax.legend()
ax.set_xlabel('Z (m)')
ax.set_ylabel('At00')

# Final LPS
fig, ax = plt.subplots()
LPSf0 = part0.loc[timesteps[-1]]
LPSf1 = part1.loc[timesteps[-1]]
ax.scatter(LPSf0.z - LPSf0.z.mean(), LPSf0.G, label='Pass 0', s=5)
ax.scatter(LPSf1.z - LPSf1.z.mean(), LPSf1.G, label='Pass 1', s=5)
ax.legend()
ax.set_xlabel('Z (m)')
ax.set_ylabel('\gamma')

# Plot efficiencies
gdf0, _ = GPT.aggregate_gdf('und_pass0.gdf')
gdf1, _ = GPT.aggregate_gdf('und_pass1.gdf')
gdf0['efficiency'] = (gdf0.avgG - gdf0.avgG.iloc[0])/gdf0.avgG.iloc[0]*100
gdf1['efficiency'] = (gdf1.avgG - gdf1.avgG.iloc[0])/gdf1.avgG.iloc[0]*100
fig, ax = plt.subplots()
ax.plot(gdf0.avgz, gdf0.efficiency, label='Pass 0')
ax.plot(gdf1.avgz, gdf1.efficiency, label='Pass 1')
ax.legend()
ax.set_xlabel('Z (m)')
ax.set_ylabel('Efficiency (%)')

