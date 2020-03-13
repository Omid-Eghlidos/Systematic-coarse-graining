#!/usr/bin/env pyhton3

import numpy

# Determine which carbn atoms are bonded to each other.
def carbon_bonds(C):
    ''' This function is to find the carbon bonding between carbon-carbon '''
    C_bonds = numpy.empty([len(C),4])
    C_bonds.fill(-11)
    for i in range(len((C))):
        counter = 0
        for j in range(len(C)):
            d = numpy.linalg.norm(C[i,:]-C[j,:])
            if i != j and d < 1.6:
                C_bonds[i,counter] = j
                counter = counter + 1
    return C_bonds

# Determine which hydrogen is connected to what carbon
def hydrogen_bonds(C,H):
    ''' This function is to find the hydrogen bonding between carbon hydrogen atoms '''
    H_bonds = numpy.empty([len(C),4])
    H_bonds.fill(-11)
    for i in range(len(C)):
        counter = 0
        for j in range(len(H)):
            dist = numpy.linalg.norm(C[i,:]-H[j,:])
            if abs(abs(dist) - 1.1) < 0.1:
                H_bonds[i,counter] = j + len(C) + 1
                counter = counter + 1
        if i % 3 == 1:
            H_bonds[i,1] = -11
    return H_bonds
        

