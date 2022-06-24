# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 07:48:27 2022

@author: afisher
"""

import subprocess
import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')
from GPT import load_gdf

cmd = subprocess.run(['gpt', '-o', 'example.gdf', 'example.in'], stdout=subprocess.PIPE)
consts, params, particle = load_gdf('example.gdf', [['z','G','t']])

# Separate screen vs tout output by checking if 't' is Nan.