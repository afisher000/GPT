# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 17:06:32 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')
import numpy as np
import pandas as pd
import subprocess
import os
import GPT
import matplotlib.pyplot as plt
plt.close('all')


# Measures offset and R56 as a function of dx and energy
data = pd.DataFrame(columns=['dx','G0','field','offset','R56','Bx'])
counter = 0
R56_i = -.05
dx = 0.005

for G0 in np.linspace(14.1,16.1,10):
    max_R56 = 0
    for field in np.linspace(.150,.300,15): #keep max field in mind
            # Create beam
            GPT.create_chirped_beam(G0, R56=R56_i)

            # Measure offset through 1 magnet
            call = f'gpt -o one_magnet.gdf one_magnet.in dx={dx} field={field}'
            cmd = subprocess.run(call.split(' '), capture_output=True, encoding='UTF-8')
            _, _, part = GPT.load_gdf('one_magnet.gdf')
            offset = part.x.mean()
            
# =============TROUBLESHOOTING=================================================
#             gdf1, part = GPT.aggregate_gdf('one_magnet.gdf')            
#             df = GPT.load_gdf('one_magnet.gdf')[2]
#             screen = df.loc[df.index.levels[0][-1]]
#             offset = screen.x.mean()
# =============================================================================
            
            # Estimate chicane R56 from two magnets
            call = f'gpt -o two_magnets.gdf two_magnets.in dx={dx} field={field} offset={offset}'
            cmd = subprocess.run(call.split(' '), capture_output=True, encoding='UTF-8')
            _, _, part = GPT.load_gdf('two_magnets.gdf')
            print(f'Xp1 = {part.Bx.mean()*1000:.1f} mrad')
            
            
            # Make first order correction to offset
            offset = part.x.mean()/2
            call = f'gpt -o two_magnets.gdf two_magnets.in dx={dx} field={field} offset={offset}'
            cmd = subprocess.run(call.split(' '), capture_output=True, encoding='UTF-8')
            _, _, part = GPT.load_gdf('two_magnets.gdf')
            
            R56_f = 1/np.polyfit(-3e8*part.t, part.G/part.G.mean(), 1)[0]
            print(f'Xp2 = {part.Bx.mean()*1000:.1f} mrad')
            
# ============TROUBLESHOOTING=====================================================
#             #gdf, part = GPT.aggregate_gdf('two_magnets.gdf')
#             #fig, ax  =plt.subplots()
#             #ax.plot(gdf1.avgz, gdf1.avgx)
#             #ax.plot(gdf.avgz, gdf.avgx)
#             #x_exp = 2*offset
#             #xerror = (gdf.avgx.iloc[-1]-x_exp)/(x_exp)
#             #print(xerror)
# =============================================================================
            
            # If R56 does not increase with field, continue
            chicane_R56 = (R56_f-R56_i)*2
            if chicane_R56<max_R56:
                break
            
            x_exp = 2*offset
            xerror = (part.x.mean()-x_exp)/x_exp
            
            # Prepare for next iteration
            data.loc[len(data)] = [dx, G0, field, offset, chicane_R56, part.Bx.mean()]
            max_R56 = chicane_R56
            print(f'Iteration = {counter}')
            counter += 1
            
# Save to file            
file = 'chicane_parameter_scan.csv'
data.to_csv(file, index=False, mode='a', header=not os.path.exists(file))



# Data-analysis plots
cmap = 'coolwarm'


data = pd.read_csv('chicane_parameter_scan.csv')

data.plot(x='G0',y='field', kind='scatter',
          c='xerror', colormap = plt.get_cmap(cmap))

data.plot(x='dx', y='R56', kind='scatter', 
          c='field', 
          colormap=plt.get_cmap(cmap))

data.plot(x='G0', y='R56', kind='scatter', 
          c='field', 
          colormap=plt.get_cmap(cmap))

data.plot(x='dx', y='offset', kind='scatter', 
          c='field', 
          colormap=plt.get_cmap(cmap))

data.plot(x='G0', y='offset', kind='scatter', 
          c='field', 
          colormap=plt.get_cmap(cmap))

