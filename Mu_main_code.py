#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 15:16:17 2022

@author: Haffari and Flinf
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
from pylab import *
import xraylib
from Mu_functions import *
import math
plt.rcParams.update({'font.size': 20})


# --------------------------------------------------------------------
# create a snowpack
#---------------------------------------------------------------------

sample_depth = 0.01                # 1 cm
number_mesh_points = 1001        #
mesh_size = sample_depth/float(number_mesh_points)

### ==== X-RAY SOURCE ==========================================
print('tomography parameters often used for snow with impregnation products:')  
I_source = 131e-6 #A => 131µA
V_source = 40e3 #V => 40 kV
print("I_source =", I_source*1e6, "microA")
print("V_source =", V_source*1e-3, "kV")
P_source_th = I_source * V_source #W
P_source_eff = P_source_th * 1 / (float(100)) #QQQ (cf wikipedia https://en.wikipedia.org/wiki/X-ray#Computed_tomography)
print("P_source_th =", P_source_th, "W")
print("P_source_eff =", P_source_eff, "W")

### alternative approach (requires K!)
#http://www.imre.ucl.ac.be/rpr/RDGN3120/faisceauRX.pdf
K = 1
Z = 74 #(Tungsten) http://www.hamamatsu.com/resources/pdf/etd/L12161-07_TXPR1023E.pdf
# P_source_eff_K = K*P_source_th*Z*V_source #QQQ (cf imre)
# print("P_source_eff_K=1 =", P_source_eff_K, "W")
###

### ==== test the csv file if empty  ==== ###
chemical_products = []
list_c = []
list_symbols = []
list_molar_mass = []
list_molar_number = []
list_symbols = []
try:
    test = pd.read_csv('NIST_list.csv')
except pd.errors.EmptyDataError:
    list_names = xraylib.GetCompoundDataNISTList()
    ## Adding elements of Mendeleiev table
    for i in range(1,99):
        symbol = xraylib.AtomicNumberToSymbol(i)
        density = xraylib.ElementDensity(i)
        list_symbols.append((symbol, symbol, density))
    ## Computing chemical formulas from the existing product of the NIST list
    for name in list_names:
        chemical_formula=''
        density = xraylib.GetCompoundDataNISTByName(name)['density']
        chemical_formula = ''.join(compute_chemical_formula(name))
        chemical_products.append((name,chemical_formula,density))
    chemical_products = pd.DataFrame.from_records(list_symbols + chemical_products, columns = ['product','chemical_formula', 'density'])
    chemical_products.to_csv('NIST_list.csv', index=False)
### ============================================

condition = True
simulation_type = input('compute absorption coefficient (µ) (1) or intensity ratio (I/I0) (2): ')

if simulation_type == '1':
    simulation_type = '\"absorption coefficient (µ)\"'
    print('Computing µ (mu)')
    mu_computation(figsize = np.array([29,21]), x_y_fontsize = 22, range_interest=True)

elif simulation_type == '2':
    simulation_type = '\"intensity ratio (I/I0)\"'
    print('Computing intensity ratio (I/I0)')
    intensity_ratio(figsize = np.array([29,21]), x_y_fontsize = 22, range_interest=True)

else:
    while condition == True:
        simulation_type = input('compute absorption coefficient (µ) (1) or intensity ratio (I/I0) (2): ')
        if simulation_type == '1':
            condition = False
            mu_computation()

        elif simulation_type == '2':
            condition = False
            intensity_ratio()


print('-----------computation of ' + simulation_type + ' finished')