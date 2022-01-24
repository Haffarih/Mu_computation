#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 12:58:38 2022

@author: haffari
"""

#Copyright (c) 2009, 2010, 2011 Tom Schoonjans
#All rights reserved.

#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#    * The names of the contributors may not be used to endorse or promote products derived from this software without specific prior written permission.

#THIS SOFTWARE IS PROVIDED BY Tom Schoonjans ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Tom Schoonjans BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""Example of using various xraylib functionality in python."""

import xraylib
import math
import numpy as np


test = xraylib.AtomicNumberToSymbol(13)
print(test)

test_x = xraylib.ElementDensity(13)
print(test_x)

test4 = xraylib.SymbolToAtomicNumber('Al')
print(test4)

test2 = xraylib.CompoundParser("Al")['Elements'][0]
print(test2)

test3 = xraylib.ElementDensity(xraylib.CompoundParser("Al")['Elements'][0])
print(test3)

list_symbols = []
for i in range(1,99):
    symbol = xraylib.AtomicNumberToSymbol(i)
    density = xraylib.ElementDensity(int(i))
    list_symbols.append((symbol,density))
print(list_symbols)


xraylib.XRayInit()

print("Example of python program using xraylib")
print("xraylib version: {}".format(xraylib.__version__))
print("Density of pure Al : {} g/cm3".format(xraylib.ElementDensity(13)))




cdtest = xraylib.CompoundParser("Ca(HCO3)2")
print("Ca(HCO3)2 contains {} atoms, {} elements and has a molar mass of {} g/mol".format(cdtest['nAtomsAll'], cdtest['nElements'], cdtest['molarMass']))


for i in range(cdtest['nElements']):
    print("Element {}: {} % and {} atoms".format(cdtest['Elements'][i], cdtest['massFractions'][i]*100.0, cdtest['nAtoms'][i]))

cdtest = xraylib.CompoundParser("SiO2")
print("SiO2 contains {} atoms, {} elements and has a molar mass of {} g/mol".format(cdtest['nAtomsAll'], cdtest['nElements'], cdtest['molarMass']))
for i in range(cdtest['nElements']):
    print("Element {}: {} % and {} atoms".format(cdtest['Elements'][i], 'mass Fraction = ' + str(cdtest['massFractions'][i]*100.0), cdtest['nAtoms'][i]))



symbol = xraylib.AtomicNumberToSymbol(26)
print("Symbol of element 26 is: {}".format(symbol))
print("Number of element Fe is: {}".format(xraylib.SymbolToAtomicNumber("Fe")))


cdn = xraylib.GetCompoundDataNISTByName("Uranium Monocarbide")
print("Uranium Monocarbide")
print("  Name: {}".format(cdn['name']))
print("  Density: {}".format(cdn['density']))
for i in range(cdn['nElements']):
    print("  Element {}: {} %%".format(cdn['Elements'][i], cdn['massFractions'][i]*100.0))

cdn = xraylib.GetCompoundDataNISTByIndex(xraylib.NIST_COMPOUND_BRAIN_ICRP)
print("NIST_COMPOUND_BRAIN_ICRP")
print("  Name: {}".format(cdn['name']))
print("  Density: {}".format(cdn['density']))
for i in range(cdn['nElements']):
    print("  Element {}: {} %%".format(cdn['Elements'][i], cdn['massFractions'][i]*100.0))

nistCompounds = xraylib.GetCompoundDataNISTList()
print("List of available NIST compounds:")
for i in range(len(nistCompounds)):
    print("  Compound {}: {}".format(i,nistCompounds[i]))

print("")

# radioNuclideData tests
rnd = xraylib.GetRadioNuclideDataByName("109Cd")
print("109Cd")
print("  Name: {}".format(rnd['name']))
print("  Z: {}".format(rnd['Z']))
print("  A: {}".format(rnd['A']))
print("  N: {}".format(rnd['N']))
print("  Z_xray: {}".format(rnd['Z_xray']))
print("  X-rays:")

------------------- END OF XRLEXAMPLE5 -------------------------")
# print("")



