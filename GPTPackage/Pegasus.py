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



def simulate_gun_to_screen4(mr, cmd_error=True):
    # Simulate gun to screen4
    write_MRfile('pegasus.mr', mr)
    gpt_call = 'mr -o gun_to_screen4.gdf pegasus.mr gpt gun_to_screen4.in'
    cmd = subprocess.run(gpt_call.split(' '), capture_output=True, encoding='UTF-8')
    if cmd_error:
        print(cmd.stderr)
    
    # Create gdf for screen at screen4
    consts, params, particle = load_gdf('gun_to_screen4.gdf')
    particle[particle.t.notnull()].to_csv('temp.txt', sep=' ', index=False)
    cmd = subprocess.run('asci2gdf -o beam_at_screen4.gdf temp.txt'.split(' '), capture_output=True, encoding='UTF-8')
    if cmd_error:
        print(cmd.stderr)
    os.remove('temp.txt')
    
    # Update quads
    if mr.estimate_quads:
        G = particle.G[particle.t.notnull()].mean()
        s_beam = particle[['x','Bx','y','By']][particle.t.notnull()].cov().values
        mr.quads = estimate_quads(s_beam, G, 5.16-3.191)
    return 


def simulate_screen4_to_und(mr, cmd_error=True):
    # Simulate screen4 to und
    write_MRfile('pegasus.mr', mr)
    gpt_call = 'mr -o screen4_to_und.gdf pegasus.mr gpt screen4_to_und.in'
    cmd = subprocess.run(gpt_call.split(' '), capture_output=True, encoding='UTF-8')
    if cmd_error:
        print(cmd.stderr)
    
    # Create gdf for screen at und entrance
    consts, params, particle = load_gdf('gun_to_screen4.gdf')
    particle[particle.t.notnull()].to_csv('temp.txt', sep=' ', index=False)
    cmd = subprocess.run('asci2gdf -o beam_at_und.gdf temp.txt'.split(' '), capture_output=True, encoding='UTF-8')
    if cmd_error:
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
    
    for z in np.linspace(0, focus, 100):
        grad = np.sum([k[0]/2* (np.tanh(b/2*(leff/2-(z-z_q4))) + np.tanh(b/2*(leff/2)+(z-z_q4))),
                       k[1]/2* (np.tanh(b/2*(leff/2-(z-z_q5))) + np.tanh(b/2*(leff/2)+(z-z_q5))),
                       k[2]/2* (np.tanh(b/2*(leff/2-(z-z_q6))) + np.tanh(b/2*(leff/2)+(z-z_q6)))])
        
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
        s_beam_final = R.dot(s_beam).dot(R)
        sigx = s_beam_final[1,1] #off by square
        sigy = s_beam_final[3,3]
        cost = sigx + sigy + abs(sigx-sigy)
        return cost
    
    opt_quads = scipy.optimize.fmin(cost_function, 
                        np.array([3,-5,3]), 
                        args=(s_beam, gammabeta, focus))
    return opt_quads