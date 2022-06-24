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
import pandas as pd

mr = read_parameters()
simulate_gun_to_screen4(mr)
simulate_screen4_to_und(mr)

_, params, particle = load_gdf('gun_to_screen4.gdf')
gdf = aggregate_gdf(params, particle)

_, params, particle = load_gdf('screen4_to_und.gdf')
gdf2 = aggregate_gdf(params, particle)