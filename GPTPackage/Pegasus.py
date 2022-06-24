# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 11:17:45 2022

@author: afisher
"""
import numpy as np
import scipy.optimize


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
        
        if grad<0:
            print(z)
            print(quads)
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
            
        R = dR.dot(R)
        
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