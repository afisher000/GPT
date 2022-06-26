# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:13:37 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\THz\\THzPackage')

from Pegasus import simulate_gun_to_screen4, simulate_screen4_to_und, simulate_und_to_end
from GPT import load_gdf, aggregate_gdf
from THz import parameters
import matplotlib.pyplot as plt

# Read Parameters
plt.close('all')
mr = parameters()


# Run simulations
simulate_gun_to_screen4(mr)
simulate_screen4_to_und(mr)
simulate_und_to_end(mr, passnum=1)

# Get data
gdf_files = ['gun_to_screen4.gdf', 'screen4_to_und.gdf', 'und_to_end.gdf']
gdf, particle = aggregate_gdf(gdf_files)
consts, _, fel = load_gdf('und_to_end.gdf', arrays_to_load=[['freq00','A00','phi00','k00','lam00']])


# Plot
gdf.plot(x='avgz', y=['stdx', 'stdy'])
gdf.plot(x='avgz', y='counts', ylim=(0,1.1*gdf.counts.max()))

