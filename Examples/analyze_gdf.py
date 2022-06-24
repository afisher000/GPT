# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 05:42:41 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')
import pandas as pd
import numpy as np
from GPT import load_gdf


# Read GDF file
arrays_to_load = [['z','G'],['freq00','A00']]
consts, params, particle, fel = load_gdf('fel.gdf', arrays_to_load=arrays_to_load)

# Aggregation calculations
counts = particle.z.groupby(level=0).count() # particles at each timestep
particle_grouped = particle.groupby(level=0) # Group by timestep
particle_means = particle_grouped.mean().add_prefix('avg')
particle_stds = particle_grouped.std().add_prefix('std')

# Construct gdf dataframe
gdf = pd.concat([params, counts, particle_means, particle_stds], axis=1) #concatenate columns
gdf['eff'] = (gdf.avgG - gdf.avgG[0])/gdf.avgG[0]*100


# Misc calculations
particle['z_rel'] = (particle.z.unstack() - gdf.avgz).stack() #adds Nan for missing particles



