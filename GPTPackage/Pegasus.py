# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 11:17:45 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')

import subprocess
import os
from GPT import load_gdf
import numpy as np
import scipy.optimize


def system(call, mr):
    ''' Execute subprocess call for GPT'''
    write_MRfile('pegasus.mr', mr)
    cmd = subprocess.run(call.split(' '), capture_output=True, encoding='UTF-8')
    print(cmd.stderr)
    
    return

    

def simulate_gun_to_screen4(mr):
    # Simulate gun to screen4
    gpt_call = 'mr -o gun_to_screen4.gdf pegasus.mr gpt gun_to_screen4.in'
    system(gpt_call, mr)
    
    # Create gdf for screen at screen4
    write_screens_to_gdf('gun_to_screen4.gdf', 'beam_at_screen4.gdf')

    # Update quads
    if mr.estimate_quads:
        _, _, particle = load_gdf('beam_at_screen4.gdf')
        G = particle.G.mean()
        s_beam = particle[['x','Bx','y','By']].cov().values
        mr.quads = estimate_quads(s_beam, G, 5.16-3.191)
    return 


def simulate_screen4_to_und(mr):
    # Simulate screen4 to und
    _, _, particle = load_gdf('beam_at_screen4.gdf')
    gpt_call = f'mr -o screen4_to_und.gdf pegasus.mr gpt screen4_to_und.in tstart={particle.t.mean()}'
    system(gpt_call, mr)
    
    # Create gdf for screen at und entrance
    write_screens_to_gdf('screen4_to_und.gdf', 'beam_at_und.gdf')
    return 


def simulate_und_to_end(mr, passnum=1):
    # Simulate und to end
    _, _, particle = load_gdf('beam_at_und.gdf')
    gpt_call = f'mr -o und_to_end.gdf pegasus.mr gpt und_to_end.in tstart={particle.t.mean()} passnum={passnum}'
    system(gpt_call, mr)
    
    # Create gdf for screen after und entrance
    write_screens_to_gdf('und_to_end.gdf', 'beam_after_und.gdf')
    return 


def write_screens_to_gdf(gdf_file, dest_files):
    ''' Find screen outputs in gdf file and save each to locations specified
    by dest_files. Length of dest_files must match number of screen outputs.'''
    # Ensure list type
    if not isinstance(dest_files, list):
        dest_files = [dest_files]
    
    # Loop over screen outputs
    _, _, particle = load_gdf(gdf_file)
    
    if 't' not in particle.columns:
        print(f'No screen outputs found in {gdf_file}')
        return
    
    for j, p_idx in enumerate(particle.t.dropna().unstack().index):
        particle.loc[p_idx].to_csv('temp.txt', sep=' ', index=False)
        cmd = subprocess.run(f'asci2gdf -o {dest_files[j]} temp.txt'.split(' '), capture_output=True, encoding='UTF-8')
        print(cmd.stderr)
        os.remove('temp.txt')
        
    return
    

def write_MRfile(file, mr):
    ''' Writes the numerical data in mr to the file specified.'''
    with open(file, 'w') as f:
        for name, value in mr.iteritems():
            if name=='linacphase':
                value = -value - 139
            elif name=='quads':
                f.write(f'I4 {value[0]:e}\n')
                f.write(f'I5 {value[1]:e}\n')
                f.write(f'I6 {value[2]:e}\n')
                
            if isinstance(value, int):
                f.write(f'{name} {value}\n')
            elif isinstance(value, float):
                f.write(f'{name} {value:e}\n')
    return


def get_quad_focusing_matrix(quads, gammabeta, focus):
    dz = .001
    R = np.eye(4)
    
    # Quad Parameters
    bRho = gammabeta*1.70451e-3;
    grad = 0.45
    b = 135
    leff = .0768
    k = np.array(quads)*grad/bRho
    
    # Quad centers from screen 4
    z_q4 = .104
    z_q5 = .104 + .086 
    z_q6 = .104 + .086 + .085
    
    for z in np.arange(0,focus,dz):
        grad = np.sum([k[0]/2* (np.tanh(b/2*(leff/2-(z-z_q4))) + np.tanh(b/2*(leff/2+(z-z_q4)))),
                       k[1]/2* (np.tanh(b/2*(leff/2-(z-z_q5))) + np.tanh(b/2*(leff/2+(z-z_q5)))),
                       k[2]/2* (np.tanh(b/2*(leff/2-(z-z_q6))) + np.tanh(b/2*(leff/2+(z-z_q6))))])
        
        if abs(grad)<1e-4:
            dR = np.array([[1, dz, 0, 0],
                           [0, 1, 0, 0],
                           [0, 0, 1, dz],
                           [0, 0, 0, 1]])
        else:
            sK = np.sqrt(grad*(1+0j)) # need complex argument
            sKL = sK*dz
            dR = np.array([[np.cos(sKL), np.sin(sKL)/sK, 0, 0],
                           [-sK*np.sin(sKL), np.cos(sKL), 0, 0],
                           [0, 0, np.cosh(sKL), np.sinh(sKL)/sK],
                           [0, 0, sK*np.sinh(sKL), np.cosh(sKL)]])
            
        R = dR.real.dot(R)
        
    return R

def estimate_quads(s_beam, gammabeta, focus):
    s_beam = np.array(s_beam)
    
    def cost_function(quads, s_beam, gammabeta, focus):
        R = get_quad_focusing_matrix(quads, gammabeta, focus)
        s_beam_final = R.dot(s_beam).dot(R.T)
        sigx = s_beam_final[0,0] #off by square
        sigy = s_beam_final[2,2]
        cost = sigx + sigy + abs(sigx-sigy)
        return cost
    
    opt_quads = scipy.optimize.fmin(cost_function, 
                        np.array([3,-5,3]), 
                        args=(s_beam, gammabeta, focus))
    return opt_quads