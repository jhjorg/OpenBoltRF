#!/usr/bin/env python
# -*- coding: utf-8 -*-

# region
# Fatigue.py helper to calculate fatigue damage from Schmidt model, FEA
__author__ = "Jack Jorgensen"
__license__ = "GNU GPLv3"
__version__ = "Per Git commit"
__maintainer__ = "Jack Jorgensen"
__email__ = "jack.jorgensen@research.uwa.edu.au"
__status__ = "Development"
# endregion

# region Import required packages and modules

# Packages
import numpy as np
import pandas as pd

# endregion

def EN1993_fat(Stress_bins, dnom, det_cat):

    # EN1993-1-9-2005 Fatigue SN curve for Details
    if det_cat == '36*':
        sigc    =   40  # Detail category
        sigd    =   23.4 # Stress at CAFL
        sigl    =   0.549*sigd # Cut-off limit
        N_c     =   2*10**6 # Cycles at detail category
        N_d     =   1*10**7 # Cycles at sigd CAFL
    elif det_cat == '50':
        sigc    =   50  # Detail category
        sigd    =   36.84 # Stress at CAFL
        sigl    =   0.549*sigd # Cut-off limit
        N_c     =   2*10**6 # Cycles at detail category
        N_d     =   5*10**6 # Cycles at sigd
    ks          =   (30/dnom)**0.25 # Size factor
    sigc_red    =   sigc*ks # Reduced detail cat
    sigd_red    =   sigd*ks # Reduced CAFL
    sigl_red    =   sigl*ks # Reduced limit
    m1          =   3       # Slope m1
    m2          =   5       # Slope m2

    # Calculate cycles to failure for each stress block
    N                 =   np.zeros((len(Stress_bins[0]),1))
    Dam_blocks_single = np.zeros((len(Stress_bins[0]), 1)) # Damage contribution from
    for i  in range(0,len(Stress_bins[0])):
        if Stress_bins[1][i] > sigd_red:
            N[i] = (N_c)*((sigc_red)**m1)/((Stress_bins[1][i])**m1)
        elif (Stress_bins[1][i] < sigd_red)&(Stress_bins[1][i] > sigl_red):
            N[i] = (N_d)*((sigd_red)**m2)/((Stress_bins[1][i])**m2)
        else :
            N[i] = 10 ** 15
        Dam_blocks_single[i] = 1 / N[i]  # Damage contribution from single cycle in that bin

    # Calculate damage
    Dam_blocks = Stress_bins[0]/N[:]  # Damage contribution from each block bin
    Cumu_Dam   = np.sum(Dam_blocks)     # Linearly cumulative damage, with 10^8 draws from samples, so damage is damage x 10 for 10^9 cycles

    return Cumu_Dam, Dam_blocks, Dam_blocks_single


    # endregion


