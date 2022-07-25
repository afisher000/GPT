# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 13:43:47 2022

@author: afisher
"""

import scipy.constants as sc
import numpy as np

R = 1.63e-3
lamu = .032
B0 = .730
ku = 2*sc.pi/lamu
eta = 1.8412
K = sc.elementary_charge*B0 / (sc.m_e * sc.c * ku)
kz = eta**2/R**2/ku
omega = sc.c * np.sqrt(kz**2 + eta**2/R**2)
Bz = sc.c*kz/omega
Gres = np.sqrt( (1+K**2)/(1-Bz**2) )
f0 = omega/2/sc.pi
