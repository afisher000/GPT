# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:13:37 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')

from Pegasus import simulate_gun_to_screen4, simulate_screen4_to_und
from GPT import load_gdf, aggregate_gdf
from params import read_parameters

# Read Parameters
mr = read_parameters()


# Run simulations
simulate_gun_to_screen4(mr, cmd_output=False)
simulate_screen4_to_und(mr, cmd_output=False)


# Get data
gdf_files = ['gun_to_screen4.gdf', 'screen4_to_und.gdf']
gdf, particle = aggregate_gdf(gdf_files)


# Plot
gdf.plot(x='avgz', y=['stdx', 'stdy'])

