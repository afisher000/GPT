# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 11:17:45 2022

@author: afisher
"""

import sys
sys.path.append('C:\\Users\\afisher\\Documents\\GitHub\\GPT\\GPTPackage')

import subprocess
import os
import GPT
import numpy as np
import scipy.optimize
import pandas as pd
from scipy.interpolate import griddata
import scipy.constants as sc
import shutil
import time

def sum_At00(wfm, dNm=1):
    ''' Computes At00 from wfm dictionary. wfm must contain data for 'amp', 
    'k', 'z', 'omega', 't', and 'phi' to compute amp*cos(kz*z-omega*t+phi).
    We assume frequency info on axis 0 and time/z info on axis 1.'''
    
    matrix = dNm * wfm['amp']*np.cos(wfm['kz']*wfm['z'] - wfm['omega']*wfm['t'] + wfm['phi'])
    return np.sum(matrix, axis=0)

def simulate(mr, start, end, undpass=0):
    ''' Run a GPT simulation of the Pegasus Beamline. Specify the start and end
    location (only adjacent pairs in break_locations allowed). The output of
    the input file should be a tout (where beam is destroyed just after end
    location) and a screen from which to start the next simulation. '''
    
    # Assert start and end are correct locations
    break_locations = ['gun','screen4','boxend','und','end']
    assert start in break_locations
    assert end in break_locations
    
    
    # Make estimations before simulation
    if mr.estimate_quads and start=='screen4':
        print('Estimating quad values from matrix model')
        _, _, particle = GPT.load_gdf(f'beam_at_{start}.gdf')
        G = particle.G.mean()
        s_beam = particle[['x','Bx','y','By']].cov().values
        mr.quads = estimate_quads(s_beam, G, mr.und_pos - mr.screen4_pos)
        
    if mr.estimate_chicane and start=='screen4':
        _, _, particle = GPT.load_gdf(f'beam_at_{start}.gdf')
        chic_scan = pd.read_csv(os.path.join('Parameter Scans','chicane_parameter_scan.csv'))
        G = particle.G.mean()
        R56 = -1/np.polyfit(-3e8*particle.t, particle.G/G, 1)[0]
        
        R56 = R56 * .8
        
        # From beam energy and R56, get field
        mr.chicane_field = griddata((chic_scan.G0, chic_scan.R56),
                                    chic_scan.field,
                                    (G, R56)).tolist()
        
        # From beam energy and field, get offset
        mr.chicane_offset = griddata((chic_scan.G0, chic_scan.field),
                                     chic_scan.offset,
                                     (G, mr.chicane_field)).tolist()
        print(f'Estimated chicane parameters:\nfield = {mr.chicane_field:.3f} T\noffset = {mr.chicane_offset:.3f} m\n')
        
        
    if mr.estimate_sol2 and start=='boxend':
        _, _, particle = GPT.load_gdf(f'beam_at_{start}.gdf')
        G = particle.G.mean()
        s_beam = particle[['x','Bx','y','By']].cov().values
        mr.sol2 = estimate_sol2(mr, s_beam, G)
        print(f'Estimated sol2 value: {mr.sol2:.3f}')
        
    if start=='und':
        mr.undpass = undpass
        if undpass>0:
            circulate_pulse(mr, undpass=undpass-1)
        
    # Get start time
    tstart = 0
    if start!='gun':
        _, _, particle = GPT.load_gdf(f'beam_at_{start}.gdf')
        tstart = particle.t.mean()
    
    # Run simulation
    print(f'Starting simulation from {start} to {end}')
    gpt_call = f'mr -o {start}_to_{end}.gdf pegasus.mr gpt {start}_to_{end}.in tstart={tstart}'
    GPT.write_MRfile('pegasus.mr', mr)
    start_timer = time.perf_counter()
    cmd = subprocess.run(gpt_call.split(' '), capture_output=True, encoding='UTF-8')
    stop_timer = time.perf_counter()
    if mr.GPT_warnings:
        print(cmd.stderr)    
    print(f'Finished simulation from {start} to {end} in {stop_timer-start_timer:.1f} seconds\n')
    
    # Save with undpass for undulator simulations
    if start=='und':
        shutil.copy('und_to_end.gdf', f'und_pass{undpass}.gdf')
        print(f'Copied "und_to_end.gdf" as "und_pass{undpass}.gdf"')
        
    # Write screen data to gdf files
    GPT.write_screens_to_gdf(f'{start}_to_{end}.gdf', f'beam_at_{end}.gdf')
    


    return         
        
def calculate_detuning(mr, first_idx=1):
    # Get averaged data from und_to_end.gdf
    consts, params, part = GPT.load_gdf('und_to_end.gdf')
    avg_part = part.groupby(level=0).mean()
    
    # Compute slippage up to 1 period in undulator
    n_timesteps = (avg_part.z<(mr.und_pos+mr.lamu)).sum()
    avg_part = avg_part.iloc[first_idx:n_timesteps]
    params = params[first_idx:n_timesteps]
    
    dz_electrons = avg_part.z.iloc[-1] - avg_part.z.iloc[0]
    dz_radiation = mr.vg*(params.time.iloc[-1] - params.time.iloc[0])
    phase_detuning = mr.kz * (dz_electrons - dz_radiation) - sc.pi/2
    print(f'Computed phase detuning: {phase_detuning:.2f} radians')
    return phase_detuning

def circulate_pulse(mr, undpass=0):
    # =============================================================================
    # # GPTFEL computes At00 phase as kz*(t00+ct+dz)-wt+phi.
    # We want phi_und such that kz*(t00+ct_und+dz)-wt_und+phi_und=kz*(t00+ct_end+dz)-wt_end+phi_end
    # So phi_und = phi_end + c*kz*(t_end-tund) - w*(t_end-t_und)
    # Shift to peak simply by sending z00->z00+pulse_center_dz
    # =============================================================================
    
    # Read data from end timestamp of simulation (where beam defined)
    consts, params, part, fel, tfel = GPT.load_gdf(f'und_pass{undpass}.gdf', 
                                                   arrays_to_load=[['z','G','t', 'Bz'], 
                                                                   ['freq00','A00','phi00','k00','lam00'], 
                                                                   ['t00','At00']])
    part = part[part.t.isna()].drop(columns=['t'])
    end_idx = part.index.get_level_values(level=0)[-1]
    t_end = params.time.loc[end_idx]
    
    fel = fel.unstack().loc[end_idx].unstack(level=0)
    tfel = tfel.unstack().loc[end_idx].unstack(level=0)
    # tfel.plot(x='t00', y='At00', title='Pulse after pass 0')
    
    
    # Read time from start of undulator simulation
    _, _, part = GPT.load_gdf('beam_at_und.gdf')
    t_und = part.t.mean()
    
    # Compute phi_und (phase for next pass)
    fel['omega'] = 2*sc.pi*fel.freq00
    fel['k0'] = fel.omega/sc.c
    fel['kz'] = np.sqrt( (fel.omega/sc.c)**2 - mr.kp**2 )
    
    
    # Apply dispersion here!
    
    pulse_center_dz = tfel.t00[tfel.At00.argmax()]
    pulse_detuning = calculate_detuning(mr)
    fel['phi_und'] = fel.phi00 + sc.c*fel.kz*(t_end-t_und) - fel.omega*(t_end-t_und) + fel.kz*pulse_center_dz - pulse_detuning
    
    # Sum new At00 to confirm centered correctly
    # wfm = {}
    # wfm['amp'] = fel.A00.values.reshape(-1,1)
    # wfm['kz'] = fel.kz.values.reshape(-1,1)
    # wfm['z'] = tfel.t00.values.reshape(1,-1) + sc.c*(t_und) + consts.dz
    # wfm['omega'] = fel.omega.values.reshape(-1,1)
    # wfm['t'] = t_und
    # wfm['phi'] = fel.phi_und.values.reshape(-1,1)
    # tfel['At_und'] = PEG.sum_At00(wfm, dNm=consts.dNm)
    # tfel.plot(x='t00', y='At_und', title='Pulse for pass 1')
    
    
    # Write to THzPulse.txt
    with open('THzPulse.txt', 'wt') as f:
        for idx in range(len(fel)):
            f.write('%8.8f\n%8.8f\n' % (fel.A00[idx]*consts.dNm, fel.phi_und[idx]))
    return



def optimize_quads():
    # Run GPT simulations from screen4 to boxend (large beamsizes will blow up in chicane!)
    # Read spotsize and divergence from gdf file
    # Cost function is that beam is same size (and same divergence?)
    # return optimized values
    pass

def optimize_sol2():
    # Run GPT simulations from boxend to undulator
    # Read spotsize from gdf file
    # Cost function is minimizing the round beam spotsize (minimize maximum spotsize?)
    # return optimized values
    pass


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
            dR = get_drift_matrix(dz)
        else:
            sK = np.sqrt(grad*(1+0j)) # need complex argument
            sKL = sK*dz
            dR = np.array([[np.cos(sKL), np.sin(sKL)/sK, 0, 0],
                           [-sK*np.sin(sKL), np.cos(sKL), 0, 0],
                           [0, 0, np.cosh(sKL), np.sinh(sKL)/sK],
                           [0, 0, sK*np.sinh(sKL), np.cosh(sKL)]])
            
        R = dR.real.dot(R)
    return R

def get_drift_matrix(L):
    M = np.array([[1,   L,  0,   0],
                  [0,   1,   0,   0],
                  [0,   0,   1,   L],
                  [0,   0,   0,   1]])
    return M    
    
def get_sol2_focusing_matrix(sol2, mr, gammabeta):
    # https://uspas.fnal.gov/materials/13Duke/SCL_Chap3.pdf
     
    #Hardcoded values from GPT implementation
    mu0 = 4*np.pi*1e-7
    L = 0.191784
    R = 0.0281232 
    Bz = sol2[0]*.357 - .0037 #aka sol2fac
    nI = Bz * np.sqrt(L**2 + 4*R**2) / mu0 / L
    bRho = gammabeta*1.70451e-3
    k = Bz / (2*bRho)
    kL2 = k*L*2
   
    
    # Drift from box_end to sol2 entrance
    L1 = mr.sol2_pos - L/2 - mr.boxend_pos
    drift1 = get_drift_matrix(L1)

    
    # Solenoid matrices
    entrancematrix = np.array([[1,   0,   0,   0],
                               [0,   1,   k,   0],
                               [0,   0,   1,   0],
                               [-k,  0,   0,   1]])
    
    solmatrix = np.array([[1,   np.sin(kL2)/(2*k),      0,  (1-np.cos(kL2))/(2*k)],
                          [0,   np.cos(kL2),            0,  np.sin(kL2)],
                          [0,   (np.cos(kL2)-1)/(2*k),  1,  np.sin(kL2)/(2*k)],
                          [0,   -np.sin(kL2),           0,  np.cos(kL2)]])
    
    exitmatrix = np.array([[1,   0,   0,   0],
                           [0,   1,   -k,  0],
                           [0,   0,   1,   0],
                           [k,   0,   0,   1]])
    
    # Drift from sol2 exit to focus
    L2 = mr.und_pos - (mr.sol2_pos+L/2) - .05
    drift2 = get_drift_matrix(L2)

    M = drift2.dot(exitmatrix).dot(solmatrix).dot(entrancematrix).dot(drift1)    
    return M
    
    

def estimate_sol2(mr, s_beam, gammabeta):
    s_beam = np.array(s_beam)
    def cost_function(sol2, mr, s_beam, gammabeta):
        M = get_sol2_focusing_matrix(sol2, mr, gammabeta)
        s_beam_final = M.dot(s_beam).dot(M.T)
        sigx = s_beam_final[0,0]
        sigy = s_beam_final[2,2]
        cost = sigx + sigy + abs(sigx-sigy)
        return cost
    
    opt_sol2 = scipy.optimize.fmin(cost_function, 1,
                                   args=(mr, s_beam, gammabeta))[0]
    return opt_sol2


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