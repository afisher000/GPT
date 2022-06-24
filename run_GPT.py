# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 07:48:27 2022

@author: afisher
"""

import subprocess
cmd = subprocess.run(['gpt', '-o', 'example.gdf', 'example.in'], stdout=subprocess.PIPE)