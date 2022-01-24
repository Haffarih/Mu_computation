#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 15:14:05 2022

@author: Haffari and Flinf
"""

#### =============== Computation of the absorption coefficient =============== ####

#### ============================================================================= ####

import sys
import matplotlib.pyplot as plt
import xraylib
import numpy as np
import math
import warnings
import pandas as pd
import matplotlib.cm as cm
warnings.filterwarnings("ignore", category=RuntimeWarning)

###-------list of all compounds available in NIST -------###
NIST_list = xraylib.GetCompoundDataNISTList()
### -----------------------------------------------------###


def compute_chemical_formula(name):
    chemical_products = []
    list_c = []
    list_symbols = []
    list_molar_mass = []
    list_molar_number = []
    list_chemical_formula = []
    for i in range(len(xraylib.GetCompoundDataNISTByName(name)['massFractions'])):
        a_mass_fraction = xraylib.GetCompoundDataNISTByName(name)['massFractions'][i]
        symbols = xraylib.AtomicNumberToSymbol(xraylib.GetCompoundDataNISTByName(name)['Elements'][i])
        molar_mass = xraylib.CompoundParser(symbols)['molarMass']
        molar_number = a_mass_fraction/molar_mass
        list_molar_number.append(molar_number)
        list_molar_mass.append(molar_mass)
        list_symbols.append(symbols)
        list_c.append(a_mass_fraction)
        list_chemical_formula.append(symbols + str(round(molar_number,4)))
    return list_chemical_formula


def find_intersection(x,y):
    try:
        s = np.abs(np.diff(np.sign(y))).astype(bool)
    except ZeroDivisionError:
        pass
    finally:
       return x[:-1][s] + np.diff(x)[s]/ (np.abs(y[1:][s]/y[:-1][s])+1)


#### =============== Computation of the coefficient of attenuation =============== ####
#### ============================================================================= ####

def mu_computation(figsize = np.array([29, 21.7]), x_y_fontsize = 18, range_interest=False):
    ######===================================================
    ###### MU: DEFINING VARUIABLES
    ######===================================================

    ###### ------- number fo products to test
    try :
        number_products = int(input('Please enter the number of materials or products to test: '))
    except ValueError:
        raise ValueError("Wrong inserted value...")

    ###### ------- several existing product
    # chemical_products = {'Chloronaphtalene': ['C10H7Cl', 1194], 'Kerosene': ['C11H21', 793.9], 'Paraffin oil': ['C15H32', 850],
    #                       'Ethyl_dodecanoate': ['C14H28O2', 862]}
    chemical_products = [('Chloronaphtalene','C10H7Cl', 1.194), ('Kerosene','C11H21', 0.7939), ('Paraffin oil','C15H32', 0.850),
                          ('Ethyl_dodecanoate','C14H28O2', 0.862)] # "name of the product", "chemical formula", "rho value"

    print('Current chemical products available in the dictionary: ', chemical_products)
    clear_products = input('clear products ? (y/n) ')

    ###### ------- test chemical products or pass through manual input
    if clear_products == 'y':
        chemical_products.clear()
    elif clear_products == 'n':
        chemical_products = chemical_products
        pass

    ###### ------- for 1st, 2nd 3rd computation in the output text
    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

    ###### ------- reading the NIST.csv file containing several products
    df = pd.read_csv('NIST_list.csv', delimiter = ',')

    ###### ------- loop over the number of products to test
    for i in range(number_products):
        name = str(input(f'name of the new {ordinal(i+1)} product: '))
        ###### --- first capital letter
        if name[0].islower():
            name= name[0].upper() + ''.join(name[1:])
        ###### --- verifying if the input name exists in NIST list
        NIST_list = xraylib.GetCompoundDataNISTList()
        list_c = []
        list_symbols = []
        list_molar_mass = []
        list_molar_number = []
        list_chemical_formula = []
        ###### --- if name exists in NIST list defined by "Xraylib": informations will be taken directly from the table
        if name in NIST_list:
            compute_chemical_formula(name)
            ## FINAL MODIF
            density = xraylib.GetCompoundDataNISTByName(name)['density']
            chemical_formula = ''.join(compute_chemical_formula(name))
            chemical_products.append((name, chemical_formula, density))
        ###### --- if name exists in the new defined NIST list: informations will be taken directly from it
        elif name in list(df['product']): 
            chemical_formula = str(list(df.loc[df['product'] == name]['chemical_formula'])[0])
            density_value = float(list(df.loc[df['product'] == name]['density'])[0])
            chemical_products.append((name, chemical_formula, density_value))
        ###### --- otherwise, passing through manual input of name, chemical formula and rho value
        else:
            single_element_or_formula = str(input('chemical formula or single element ?: (f/e) '))
            ###### --- "f" stands for the available chemical formula
            if single_element_or_formula == 'f':
                chemical_formula = input('Chemical formula: ')
                try:
                    rho_value = float(input('Please insert rho value: ')) # 850 kg/m3
                    density_value = rho_value/1000.
                except ValueError:
                    pass
                if len(name) != 0 and len(chemical_formula) !=0 and len(str(density_value)) !=0:
                    chemical_products.append((name, chemical_formula, density_value))
                else:
                    print('please recheck your insertion or insert \'none\'')
                    if input('If none required : please insert \'none\': ') =='none':
                        pass
                    elif input('If none required : please insert \'none\': ') =='':
                        print('please recheck your insertion or insert \'none\'')
                    else:
                        while len(name) == 0 and len(chemical_formula) ==0 and len(str(density_value)) !=0:
                            name= input('Please choose your impregnation product')
                            chemical_formula = input('Please insert the chemical formula')
                            chemical_products.append((name, chemical_formula, density_value))
            ###### --- "e" stands for the available chemical element
            elif single_element_or_formula =='e':
                try:
                    chemical_formula = input('element: ')
                    # rho_value = float(input('Please insert rho value: ')) # 850 kg/m3
                    # density_value = rho_value/1000
                    Element_number = xraylib.CompoundParser(chemical_formula)['Elements'][0]
                    density_value = xraylib.ElementDensity(Element_number)
                    rho_value = float(density_value*1000)
                except ValueError:
                    pass
                if len(name) != 0 and len(str(chemical_formula)) !=0 and len(str(density_value)) !=0:
                    chemical_products.append((name, chemical_formula, density_value))
                    # chemical_products_df = pd.DataFrame(chemical_products,columns = ['product','chemical_formula', 'density'])
                    # chemical_products_Nist_list = pd.read_csv('NIST_list.csv', delimiter = ',')
                    # chemical_products_new = chemical_products_Nist_list.append(chemical_products_df)
                    # chemical_products_new.to_csv('NIST_list.csv', index=False)
                else:
                    print('please recheck your insertion or insert \'none\'')
                    if input('If none required : please insert \'none\': ') =='none':
                        pass
                    elif input('If none required : please insert \'none\': ') =='':
                        print('please recheck your insertion or insert \'none\'')
                    else:
                        while len(name) == 0 and len(str(chemical_formula)) ==0 and len(str(density_value)) !=0:
                            name= input('Please choose your impregnation product')
                            chemical_formula = input('Please insert the chemical formula')



        pass
        if name in list(df['product']):
            pass
        else:
            chemical_products_df = pd.DataFrame(chemical_products,columns = ['product','chemical_formula', 'density'])
            chemical_products_Nist_list = pd.read_csv('NIST_list.csv', delimiter = ',')
            chemical_products_new = chemical_products_Nist_list.append(chemical_products_df).drop_duplicates(subset=['product']).reset_index(drop=True)
            chemical_products_new.to_csv('NIST_list.csv', index=False)

    ######===================================================
    ###### MU: BEGIN COMPUATION
    ###### UNIT = cm²/g*g/cm³ = cm^{-1} (ou cm²/g si on considère que la densité est sans unité : MAIS c'est trompeur car c'est difficile de distinguer \mu/\rho de \mu...)
    ######===================================================
    mesh_size_E = 1.0 #keV #should always be > 1 to prevent plotting problems
    E_l_min = 4 #keV => 40 kV
    E_l_max = 80 #keV => 150 kV
    number_mesh_points_E = int ((E_l_max-E_l_min)/mesh_size_E)+1
    E=[0.]*number_mesh_points_E
    E_l = E_l_min
    
    fig=plt.figure(figsize=figsize*0.39)
    ax=plt.subplot(1,1,1)
    ax.set_ylim(-0.25, 20)
    ax.set_xlim(2.5, 80)
    # color = cm.rainbow(np.linspace(0, 1, number_products))
    color = cm.rainbow(np.linspace(0, 1, len(chemical_products)))
    mu_list={}

    # for keys, values in chemical_products.items():
    for items in range(len(chemical_products)):
        E_l = E_l_min
        # mu_list[keys] = [0.]*number_mesh_points_E
        mu_list[items] = [0.]*number_mesh_points_E
        # print(mu_list)
        while(E_l <= E_l_max):
            if (E_l >= E_l_min) and (E_l <= E_l_max):
                mu_list[items][int((E_l - E_l_min)/mesh_size_E)] = xraylib.CS_Total_CP(chemical_products[items][1], E_l)#*chemical_products[items][2]
                E[int ((E_l-E_l_min)/mesh_size_E)] = E_l
                min_vertical = 13
                max_vertical = 20
                z_min = find_intersection(np.array(mu_list[items]), np.array(E)-min_vertical)
                z_max = find_intersection(np.array(mu_list[items]), np.array(E)-max_vertical)
            E_l = E_l + mesh_size_E
        # for c in color:
        plt.plot(E, mu_list[items], label = chemical_products[items][0], color = color[items])
        plt.plot(np.ones_like(z_min)*min_vertical, z_min, marker="o", ls="", ms=6, color=color[items])
        plt.plot(np.ones_like(z_max)*max_vertical, z_max, marker="o", ls="", ms=6, color=color[items])
        print('---mu of {} is equal to {} for 13 keV and {} for 20 keV (per cm)'.format(chemical_products[items][0], round(float(z_min[1]),2), round(float(z_max[1]),2)))

    ax.set_xlabel("Energy (keV)\nConversion into kV => E(kV) = E(keV) *2 or *3", fontsize=x_y_fontsize)
    ax.set_ylabel('$\mu$ (cm$^{-1}$)', fontsize=x_y_fontsize)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(0.0)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(0.0)
    ax.tick_params(axis = 'x', width=1.5)
    ax.tick_params(axis = 'y', width=1.5)

    if range_interest == True :
        ax.axvline(13, color='k', linestyle = 'dashed', linewidth= 1.)
        ax.axvline(20, color='k', linestyle = 'dashed', linewidth= 1.)
    plt.tight_layout()
    ax.legend(fontsize=10, loc='best')
    plt.savefig('Fig_mu'+'.png', dpi=300, bbox_inches='tight')


#### ================================= END of mu ================================= ####
#### ============================================================================= ####



#### ============================== Computation of I/I0 ========================== ####
#### ============================================================================= ####
def intensity_ratio(figsize = np.array([29, 21.7]), x_y_fontsize = 18, range_interest=False):
#######################################
    # density_graphite = 2.26 # https://www.entegris.com/content/dam/web/resources/brochures/brochure-properties-and-characteristics-of-graphite-7329.pdf
    # density_Al = 2.7 # https://www.thyssenkrupp-materials.co.uk/density-of-aluminium.html
    # density_Fe = 7.874 # https://www.vis-express.fr/fr/content/40-densite-des-metaux
    # density_Cu = 8.96 # https://en.wikipedia.org/wiki/Copper
    # density_Si_CH3_2_O = 1.04 #https://www.merckmillipore.com/FR/fr/product/Silicone-oil,MDA_CHEM-107742
    # density_PP = 0.85 #https://fr.wikipedia.org/wiki/Polypropyl%C3%A8ne
    # density_PMMA = 1.18 # https://fr.wikipedia.org/wiki/Polym%C3%A9thacrylate_de_m%C3%A9thyle
    # density_Dodecane = 0.750 #https://fr.wikipedia.org/wiki/Dod%C3%A9cane
    # density_kerosene = 0.800 #https://fr.wikipedia.org/wiki/K%C3%A9ros%C3%A8ne
#######################################
    # materials = {'ice': ['H2O', 917, 0.3*0.6],
    #              'PP': ['C3H6', 850, 0.1*2],
    #              'Al': [xraylib.SymbolToAtomicNumber("Al"), 2700, 0.1*2],'PMMA': ['C5O2H8', 1180, 0.1*2],
    #              'Graphite': [xraylib.SymbolToAtomicNumber("C"), 2260, 0.1*2],
    #              'Fe': [xraylib.SymbolToAtomicNumber("Fe"), 7874, 0.1*2],
    #              'PTFE': ['C2F4', 2200, 0.1*2]
    #              }
    
    materials = [('ice','H2O', 0.917, 0.3*0.6),
                 ('PP','C3H6', 0.90, 0.1*2),
                 ('PE', 'C2H4', 0.88, 0.1*2),
                 ('PMMA','C5O2H8', 1.18, 0.1*2),
                 ('Al',"Al", 2.700, 0.1*2),
                 ('Graphite',"C", 2.260, 0.1*2),
                 ('Fe',"Fe", 7.874, 0.1*2),
                 ('PTFE','C2F4', 2.200, 0.1*2)
                 ]

    try :
        number_products = int(input('Please enter the number of materials or products to test: '))
    except ValueError:
        raise ValueError("No or wrong added value...")

    print('Current chemical products available in the dictionary: ', materials)
    clear_products = input('clear products ? (y/n) ')
    if clear_products == 'y':
        materials.clear()
    elif clear_products == 'n':
        materials = materials
        pass
        

    df = pd.read_csv('NIST_list.csv', delimiter = ',')
    # print(df['product'])

    ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])

    for i in range(number_products):
        name = str(input(f'name of the new {ordinal(i+1)} product: '))
        if name[0].islower():
            name= name[0].upper() + ''.join(name[1:])
        ##" ===> verifying of name exists in NIST list
        NIST_list = xraylib.GetCompoundDataNISTList()
        list_c = []
        list_symbols = []
        list_molar_mass = []
        list_molar_number = []
        list_chemical_formula = []
        
        if name in NIST_list:
            compute_chemical_formula(name)
            density = xraylib.GetCompoundDataNISTByName(name)['density']
            chemical_formula = ''.join(compute_chemical_formula(name))
            thickness =float(input('Please insert thickness of the material in cm: '))
            materials.append((name, chemical_formula, density, thickness))
        
        elif name in list(df['product']): 
            chemical_formula = str(list(df.loc[df['product'] == name]['chemical_formula'])[0])
            density_value = float(list(df.loc[df['product'] == name]['density'])[0])
            thickness =float(input('Please insert thickness of the material in cm: '))
            materials.append((name, chemical_formula, density_value, thickness))
        ### ====> if not: passing through manual insertion
        
        else:
            single_element_or_formula = str(input('chemical formula or single element ?: (f/e) '))
    
            if single_element_or_formula == 'f':
                chemical_formula = input('Chemical formula: ')
                try:
                    rho_value = float(input('Please insert rho value: ')) # 850 kg/m3
                    density_value = rho_value/1000.
                    thickness =float(input('Please insert thickness of the material in cm: '))
                except ValueError:
                    pass
                # rho_value = float(input('rho value: '))
                if len(name) != 0 and len(chemical_formula) !=0 and len(str(density_value)) !=0 and len(str(thickness)) !=0:
                    # print('Insertion OK')
                    materials.append((name, chemical_formula, density_value, thickness))
                else:
                    print('please recheck your insertion or insert \'none\'')
                    if input('If none required : please insert \'none\': ') =='none':
                        pass
                    elif input('If none required : please insert \'none\': ') =='':
                        print('please recheck your insertion or insert \'none\'')
                    else:
                        while len(name) == 0 and len(chemical_formula) ==0 and len(str(density_value)) !=0 and len(str(thickness)) !=0:
                            name= input('Please choose your impregnation product')
                            chemical_formula = input('Please insert the chemical formula')
                            thickness =float(input('Please insert thickness of the material in cm: '))
                            materials.append((chemical_formula, density_value, thickness))
    
            elif single_element_or_formula =='e':
                try:
                    chemical_formula = input('element: ')
                    Element_number = xraylib.CompoundParser(chemical_formula)['Elements'][0]
                    density_value = xraylib.ElementDensity(Element_number)
                    rho_value = float(density_value*1000)
                    thickness =float(input('Please insert thickness of the material in cm: '))
                except ValueError:
                    pass
                if len(name) != 0 and len(str(chemical_formula)) !=0 and len(str(rho_value)) !=0 and len(str(thickness)) !=0:
                    # print('Insertion OK')
                    materials.append((name,chemical_formula, density_value, thickness))
                else:
                    print('please recheck your insertion or insert \'none\'')
                    if input('If none required : please insert \'none\': ') =='none':
                        pass
                    elif input('If none required : please insert \'none\': ') =='':
                        print('please recheck your insertion or insert \'none\'')
                    else:
                        while len(name) == 0 and len(str(chemical_formula)) ==0 and len(str(rho_value)) !=0 and len(str(thickness)) !=0:
                            name= input('Please choose your impregnation product')
                            chemical_formula = input('Please insert the chemical formula')
                            thickness =float(input('Please insert thickness of the material in cm: '))
                            materials.append((name,chemical_formula, density_value, thickness))

        if name in list(df['product']):
            pass
        else:
            chemical_products_df = pd.DataFrame(materials,columns = ['product','chemical_formula', 'density', 'thickness'])
            chemical_products_Nist_list = pd.read_csv('NIST_list.csv', delimiter = ',')
            chemical_products_new = chemical_products_Nist_list.append(chemical_products_df).drop_duplicates(subset=['product']).reset_index(drop=True)
            chemical_products_new.to_csv('NIST_list.csv', index=False)

    #===================================================
    # I/I0 COMPUTATION: BEGIN
    #===================================================
    mesh_size_E = 1.0 #keV #should always be > 1 to prevent plotting problems
    E_l_min = 4 #keV => 40 kV
    E_l_max = 80 #keV => 150 kV
    number_mesh_points_E = int ((E_l_max-E_l_min)/mesh_size_E)+1
    E=[0.]*number_mesh_points_E
    E_l = E_l_min
    fig=plt.figure(figsize=figsize*0.39)
    ax=plt.subplot(1,1,1)
    ax.set_ylim(-0.25, 4.5)
    ax.set_xlim(2.5, 42)

    color = cm.rainbow(np.linspace(0, 1, len(materials)))

    mu_list={}
    for items in range(len(materials)):
        E_l = E_l_min
        mu_list[items] = [0.]*number_mesh_points_E
        while(E_l <= E_l_max):
            if (E_l >= E_l_min) and (E_l <= E_l_max):
                mu_list[items][int((E_l - E_l_min)/mesh_size_E)] =math.exp(-xraylib.CS_Total_CP(materials[items][1], E_l)*(materials[items][2])*materials[items][3]) # [0]: chemical formula or element number; [1]: rho value; [2]: thickness
                E[int ((E_l-E_l_min)/mesh_size_E)] = E_l
                min_vertical = 13
                max_vertical = 20
                z_min = find_intersection(np.array(mu_list[items]), np.array(E)-min_vertical)
                z_max = find_intersection(np.array(mu_list[items]), np.array(E)-max_vertical)
            E_l = E_l + mesh_size_E

        # ax.plot(E, mu_list[keys], label = keys)
        plt.plot(E, mu_list[items], label = materials[items][0], color = color[items])
        plt.plot(np.ones_like(z_min)*min_vertical, z_min, marker="o", ls="", ms=6, color=color[items])
        plt.plot(np.ones_like(z_max)*max_vertical, z_max, marker="o", ls="", ms=6, color=color[items])
        # print('---Intensity ratio (I/I0) of {} is equal to {} for 13 keV and {} for 20 keV (per cm)'.format(keys, round(float(z_min[1]),2), round(float(z_max[1]),2)))
        print('---Intensity ratio (I/I0) of {} is equal to {} for 13 keV and {} for 20 keV (per cm)'.format(materials[items][0], round(float(z_min[1]),2), round(float(z_max[1]),2)))


        ax.set_xlabel("Energy (keV)\nConversion into kV => E(kV) = E(keV) *2 or *3", fontsize=x_y_fontsize)
        ax.set_ylabel('I/I$_0$ (-)', fontsize=x_y_fontsize)
        plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize = 10)
        plt.tight_layout()
        ax.set_xlim(E_l_min, E_l_max)
        ax.set_ylim(0,1)



        if range_interest == True :
            ax.axvline(13, color='k', linestyle = 'dashed', linewidth= 1.)
            ax.axvline(20, color='k', linestyle = 'dashed', linewidth= 1.)
            
        plt.savefig('Fig_I_I0'+'.png', dpi=300, bbox_inches='tight')
        
#### ================================= END of I/I0 =============================== ####
#### ============================================================================= ####




# #### ============================== Computation of I/I0 ========================== ####
# #### ============================================================================= ####
# def intensity_ratio(figsize = np.array([29, 21.7]), x_y_fontsize = 18, range_interest=False):
# #######################################
#     # density_graphite = 2.26 # https://www.entegris.com/content/dam/web/resources/brochures/brochure-properties-and-characteristics-of-graphite-7329.pdf
#     # density_Al = 2.7 # https://www.thyssenkrupp-materials.co.uk/density-of-aluminium.html
#     # density_Fe = 7.874 # https://www.vis-express.fr/fr/content/40-densite-des-metaux
#     # density_Cu = 8.96 # https://en.wikipedia.org/wiki/Copper
#     # density_Si_CH3_2_O = 1.04 #https://www.merckmillipore.com/FR/fr/product/Silicone-oil,MDA_CHEM-107742
#     # density_PP = 0.85 #https://fr.wikipedia.org/wiki/Polypropyl%C3%A8ne
#     # density_PMMA = 1.18 # https://fr.wikipedia.org/wiki/Polym%C3%A9thacrylate_de_m%C3%A9thyle
#     # density_Dodecane = 0.750 #https://fr.wikipedia.org/wiki/Dod%C3%A9cane
#     # density_kerosene = 0.800 #https://fr.wikipedia.org/wiki/K%C3%A9ros%C3%A8ne
# #######################################
#     materials = {'ice': ['H2O', 917, 0.3*0.6],
#                  'PP': ['C3H6', 850, 0.1*2],
#                  'Al': [xraylib.SymbolToAtomicNumber("Al"), 2700, 0.1*2],'PMMA': ['C5O2H8', 1180, 0.1*2],
#                  'Graphite': [xraylib.SymbolToAtomicNumber("C"), 2260, 0.1*2],
#                  'Fe': [xraylib.SymbolToAtomicNumber("Fe"), 7874, 0.1*2],
#                  'PTFE': ['C2F4', 2200, 0.1*2]
#                  }

#     try :
#         number_products = int(input('Please enter the number of materials or products to test: '))
#     except ValueError:
#         raise ValueError("No or wrong added value...")

#     print('Current chemical products available in the dictionary: ', materials)
#     clear_products = input('clear products ? (y/n) ')
#     if clear_products == 'y':
#         materials.clear()

#     ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
#     for i in range(number_products):
#         name = str(input(f'name of the new {ordinal(i+1)} product: '))
#         if name[0].islower():
#             name= name[0].upper() + ''.join(name[1:])
#         ##" ===> verifying of name exists in NIST list
#         NIST_list = xraylib.GetCompoundDataNISTList()
#         list_c = []
#         list_symbols = []
#         list_molar_mass = []
#         list_molar_number = []
#         list_chemical_formula = []
#         if name in NIST_list:
#             for i in range(len(xraylib.GetCompoundDataNISTByName(name)['massFractions'])):
#                 a_mass_fraction = xraylib.GetCompoundDataNISTByName(name)['massFractions'][i]
#                 symbols = xraylib.AtomicNumberToSymbol(xraylib.GetCompoundDataNISTByName(name)['Elements'][i])
#                 molar_mass = xraylib.CompoundParser(symbols)['molarMass']
#                 molar_number = a_mass_fraction/molar_mass
#                 list_molar_number.append(molar_number)
#                 list_molar_mass.append(molar_mass)
#                 list_symbols.append(symbols)
#                 list_c.append(a_mass_fraction)
#                 list_chemical_formula.append(symbols + str(round(molar_number,4)))
                
#             density = xraylib.GetCompoundDataNISTByName(name)['density']*1000
#             print('mass fraction: ', list_c)
#             print('Symbols : ', list_symbols)
#             print('Molar mass: ', list_molar_mass)
#             print('Molar number: ', list_molar_number)
#             print('Density: ', density)
#             print('Chemical formula: ', list_chemical_formula)
#             chemical_formula = ''.join(list_chemical_formula)
#             thickness =float(input('Please insert thickness of the material in cm: '))
#             materials[name] = [chemical_formula, density, thickness]
#         ### ====> if not: passing through manual insertion
        
#         else:
#             single_element_or_formula = str(input('chemical formula or single element ?: (f/e) '))
    
#             if single_element_or_formula == 'f':
#                 chemical_formula = input('Chemical formula: ')
#                 try:
#                     rho_value = float(input('Please insert rho value: ')) # 850 kg/m3
#                     density_value = rho_value/1000.
#                     thickness =float(input('Please insert thickness of the material in cm: '))
#                 except ValueError:
#                     pass
#                 # rho_value = float(input('rho value: '))
#                 if len(name) != 0 and len(chemical_formula) !=0 and len(str(rho_value)) !=0 and len(str(thickness)) !=0:
#                     # print('Insertion OK')
#                     materials[name] = [chemical_formula, rho_value, thickness]
#                 else:
#                     print('please recheck your insertion or insert \'none\'')
#                     if input('If none required : please insert \'none\': ') =='none':
#                         pass
#                     elif input('If none required : please insert \'none\': ') =='':
#                         print('please recheck your insertion or insert \'none\'')
#                     else:
#                         while len(name) == 0 and len(chemical_formula) ==0 and len(str(rho_value)) !=0 and len(str(thickness)) !=0:
#                             name= input('Please choose your impregnation product')
#                             chemical_formula = input('Please insert the chemical formula')
#                             thickness =float(input('Please insert thickness of the material in cm: '))
#                             materials[name] = [chemical_formula, rho_value, thickness]
    
#             elif single_element_or_formula =='e':
#                 chemical_formula = input('element: ')
#                 try:
#                     rho_value = float(input('Please insert rho value: ')) # 850 kg/m3
#                     density_value = rho_value/1000.
#                     thickness =float(input('Please insert thickness of the material in cm: '))
#                 except ValueError:
#                     pass
#                 # rho_value = float(input('rho value: '))
#                 if len(name) != 0 and len(str(chemical_formula)) !=0 and len(str(rho_value)) !=0 and len(str(thickness)) !=0:
#                     # print('Insertion OK')
#                     materials[name] = [chemical_formula, rho_value, thickness]
#                 else:
#                     print('please recheck your insertion or insert \'none\'')
#                     if input('If none required : please insert \'none\': ') =='none':
#                         pass
#                     elif input('If none required : please insert \'none\': ') =='':
#                         print('please recheck your insertion or insert \'none\'')
#                     else:
#                         while len(name) == 0 and len(str(chemical_formula)) ==0 and len(str(rho_value)) !=0 and len(str(thickness)) !=0:
#                             name= input('Please choose your impregnation product')
#                             chemical_formula = input('Please insert the chemical formula')
#                             thickness =float(input('Please insert thickness of the material in cm: '))
#                             materials[name] = [chemical_formula, rho_value, thickness]
#     #===================================================
#     # I/I0 COMPUTATION: BEGIN
#     #===================================================
#     mesh_size_E = 1.0 #keV #should always be > 1 to prevent plotting problems
#     E_l_min = 4 #keV => 40 kV
#     E_l_max = 80 #keV => 150 kV
#     number_mesh_points_E = int ((E_l_max-E_l_min)/mesh_size_E)+1
#     E=[0.]*number_mesh_points_E
#     E_l = E_l_min
#     #(in cm)
#     # 0.393701 corresponds to inch to cm conversion
#     # h = 29.7 * 0.393701
#     # w = 21 * 0.393701
#     fig=plt.figure(figsize=figsize*0.39)
#     ax=plt.subplot(1,1,1)
#     ax.set_ylim(-0.25, 4.5)
#     ax.set_xlim(2.5, 42)

#     mu_list={}
#     for keys, values in materials.items():
#         E_l = E_l_min
#         mu_list[keys] = [0.]*number_mesh_points_E
#         while(E_l <= E_l_max):
#             if (E_l >= E_l_min) and (E_l <= E_l_max):
#                 mu_list[keys][int((E_l - E_l_min)/mesh_size_E)] =math.exp(-xraylib.CS_Total_CP(values[0], E_l)*(values[1]/1000)*values[2]) # [0]: chemical formula or element number; [1]: rho value; [2]: thickness
#                 E[int ((E_l-E_l_min)/mesh_size_E)] = E_l
#                 min_vertical = 13
#                 max_vertical = 20
#                 z_min = find_intersection(np.array(mu_list[keys]), np.array(E)-min_vertical)
#                 z_max = find_intersection(np.array(mu_list[keys]), np.array(E)-max_vertical)
#             E_l = E_l + mesh_size_E

#         ax.plot(E, mu_list[keys], label = keys)
#         ax.set_xlabel("Energy (keV)\nConversion into kV => E(kV) = E(keV) *2 or *3", fontsize=x_y_fontsize)
#         ax.set_ylabel('I/I$_0$ (-)', fontsize=x_y_fontsize)
#         plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize = 10)
#         plt.tight_layout()
#         ax.set_xlim(E_l_min, E_l_max)
#         ax.set_ylim(0,1)
#         plt.plot(np.ones_like(z_min)*min_vertical, z_min, marker="o", ls="", ms=6, color="limegreen")
#         plt.plot(np.ones_like(z_max)*max_vertical, z_max, marker="o", ls="", ms=6, color="limegreen")
#         print('---Intensity ratio (I/I0) of {} is equal to {} for 13 keV and {} for 20 keV (per cm)'.format(keys, round(float(z_min[1]),2), round(float(z_max[1]),2)))

#         if range_interest == True :
#             ax.axvline(13, color='k', linestyle = 'dashed', linewidth= 1.)
#             ax.axvline(20, color='k', linestyle = 'dashed', linewidth= 1.)
            
#         plt.savefig('Fig_I_I0'+'.png', dpi=300, bbox_inches='tight')
        
# #### ================================= END of I/I0 =============================== ####
# #### ============================================================================= ####