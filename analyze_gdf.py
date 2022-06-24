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
param, part, fel = load_gdf('fel.gdf', arrays_to_load=arrays_to_load)
param = param.loc[0] # if values do not change over timesteps

# Aggregation calculations
counts = part.z.groupby(level=0).count()
part_grouped = part.groupby(level=0)
part_means = part_grouped.mean().add_prefix('avg')
part_stds = part_grouped.std().add_prefix('std')


# Construct gdf dataframe
gdf = pd.concat([counts, part_means, part_stds], axis=1)
gdf['eff'] = (gdf.avgG - gdf.avgG[0])/gdf.avgG[0]*100


# Misc calculations
part['z_rel'] = (part.z.unstack() - gdf.avgz).stack() #adds Nan for missing particles



