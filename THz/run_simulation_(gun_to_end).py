# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:13:37 2022

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


# Read Parameters
# =============================================================================
plt.close('all')
mr = THz.parameters()


# Run simulations
# =============================================================================
mr.fel_flag = 1
mr.GPT_warnings = 0

# Why does chicane_estimate need fudge factor to work right?
PEG.simulate(mr, 'gun','screen4')
PEG.simulate(mr, 'screen4','boxend')
PEG.simulate(mr, 'boxend','und')
PEG.simulate(mr, 'und', 'end', undpass=0)


# Get data
# =============================================================================
print('Loading data')
locs = ['gun','screen4','boxend','und','end']
gdf_files = [f'{start}_to_{end}.gdf' for start, end in zip(locs[:-1], locs[1:])]

gdf, particle = GPT.aggregate_gdf(gdf_files)
consts, _, fel = GPT.load_gdf('und_to_end.gdf', arrays_to_load=[['freq00','A00','phi00','k00','lam00']])
chicgdf, chicpart = GPT.aggregate_gdf('screen4_to_boxend.gdf')
undgdf, undpart = GPT.aggregate_gdf('und_to_end.gdf')

# Plots
print('Making plots')
# =============================================================================
gdf.plot(x='avgz', y=['stdx', 'stdy'])
gdf.plot(x='avgz', y='counts', ylim=(0,1.1*gdf.counts.max()))
gdf.plot(x='avgz', y=['avgx','avgy'])
gdf.plot(x='avgz', y='stdz')

# At undulator entrance
get_bunching = lambda z: np.abs(np.sum(np.exp(1j*mr.kz*z)))/len(z)
undpart.loc[0].plot.scatter(x='x', y='y')
undpart.loc[0].plot.scatter(x='z', y='G')
chicgdf['bunching'] = chicpart.z.groupby(level=0).agg(get_bunching)
chicgdf.plot(x='avgz',y='bunching')

# Deceleration
Gi = undgdf.avgG.iloc[0]
undgdf['efficiency'] = 100*(undgdf.avgG - Gi)/Gi
undgdf.plot(x='avgz', y='efficiency', ylabel='Efficiency (%)')

# Initial and final LPS
LPSi = undpart.loc[0]
timesteps = undpart.index.get_level_values(level=0).unique()
LPSf = undpart.loc[timesteps[-2]]
fig, ax = plt.subplots()
ax.scatter(LPSi.z-LPSi.z.mean(), LPSi.G, label='Initial', s = 5)
ax.scatter(LPSf.z-LPSf.z.mean(), LPSf.G, label='Final', s=5)
ax.set_xlabel('Z (m)')
ax.set_ylabel('\gamma')
ax.legend()

print('\nSimulation Results:')
print(f'Bunching = {chicgdf.bunching.iloc[-1]:.2f}')
print(f'Efficiency = {undgdf.efficiency.iloc[-1]:.1f}%')

