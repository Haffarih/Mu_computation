#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Copyright (c) 2009, 2010, 2011 Tom Schoonjans
#All rights reserved.

#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#    * The names of the contributors may not be used to endorse or promote products derived from this software without specific prior written permission.

#THIS SOFTWARE IS PROVIDED BY Tom Schoonjans ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Tom Schoonjans BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
Created on Fri Jan 21 12:58:38 2022

@author: haffari
"""
"""Example of using various xraylib functionality in python."""

import xraylib
import math
import numpy as np

## converting atomic numbr to symbol
symbol = xraylib.AtomicNumberToSymbol(13)
print(symbol)

## converting symbol of an element to its atmic number
atomic_number = xraylib.SymbolToAtomicNumber('Al')
print(atomic_number)

## calculating the density of an element: requires atomic number
density_element = xraylib.ElementDensity(13)
print(density_element)

## compound parser: (1) calculating characteristic of an element 
compound_element = xraylib.CompoundParser('Al')
print(compound_element)
## compound parser: (2) calculating characteristic of a chemical formula
compound_formula = xraylib.CompoundParser('C10H7Cl')
print(compound_formula)

## example of the properties of a chemical product
compound = "H2O"
compound_properties = xraylib.CompoundParser("H2O")
print("{} contains {} atoms, {} elements and has a molar mass of {} g/mol".format(compound, compound_properties['nAtomsAll'], compound_properties['nElements'], compound_properties['molarMass']))
for i in range(compound_properties['nElements']):
    print("Element {} ({}): {} % and {} atoms".format(compound_properties['Elements'][i],xraylib.AtomicNumberToSymbol(compound_properties['Elements'][i]), 'mass Fraction = ' + str(round(compound_properties['massFractions'][i]*100.0,2)), compound_properties['nAtoms'][i]))

## example of all products available in the NIST list
nist_compounds = xraylib.GetCompoundDataNISTList()
print("List of available NIST compounds:")
for i in range(len(nist_compounds)):
    print("  Compound {}: {}".format(i,nist_compounds[i]))


## calling a product of the NIST list
compound = xraylib.GetCompoundDataNISTByName("Polypropylene")
print("  Name: {}".format(compound['name']))
print("  Density: {}".format(compound['density']))
for i in range(compound['nElements']):
    print("  Element {} ({}): {} %".format(compound['Elements'][i],xraylib.AtomicNumberToSymbol(compound['Elements'][i]) , round(compound['massFractions'][i]*100.0),2))


##------------------ END OF XRLEXAMPLE5 -------------------------")




