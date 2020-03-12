#!/usr/bin/env python3

import numpy
import math
from matplotlib import pyplot
import output as out
import plotipp as plt
import applyH
import bonds
import angles as ang
import dihedrals as dihed

''' Input the size of the system in fractional coordinates '''
a = 1
b = 1
c = 1 

# Turn on(key = 1) and off (key = 0) plotting the figures
key = 0

''' output filename '''
crystallinity = 'Alpha'
filename = crystallinity+'_iPP_'+str(a)+'a'+str(b)+'b'+str(c)+'c'
filename = filename + '.data'

# Creating the backbone carbons in a unit cell
''' Atomic fractional coordiantes x/a, y/b, and z/c
from (1988 Macromolecules Immiczi and Iannelil) 
Comments indicate bonding. '''

C1 = numpy.array([[-0.0727, 0.2291, 0.2004],  #C3 0: 1
                  [-0.0765, 0.1592, 0.2788],  #C1 1: 0 2
                  [-0.1021, 0.1602, 0.5098],  #C2 2: 1 4
                  [-0.3087, 0.0589, 0.4941],  #C3 3: 4
                  [-0.1146, 0.0928, 0.6057],  #C1 4: 2 3 5
                  [-0.1044, 0.0854, 0.8428],  #C2 5: 4 7
                  [ 0.2775, 0.0797, 0.9260],  #C3 6: 7
                  [ 0.0872, 0.1156, 0.9730],  #C1 7: 5 6 8
                  [ 0.1026, 0.1221, 1.2109]]) #C2 8: 7

def apply_space_group(C1):
    ''' Adding chains in a unitcell using symmetry operations '''
    C2, C3, C4 = C1.copy(), C1.copy(), C1.copy()

    C2[:,0] =  C2[:,0]
    C2[:,1] = -C2[:,1]
    C2[:,2] =  C2[:,2] - 0.5

    C3[:,0] =  C3[:,0] - 0.5
    C3[:,1] =  C3[:,1] - 0.5
    C3[:,2] =  C3[:,2] 

    C4[:,0] =  C3[:,0] 
    C4[:,1] = -C3[:,1] 
    C4[:,2] =  C2[:,2]    
    
    return C1, C2, C3, C4

''' Angle between c vector and x-axis. '''
beta = 99.5 * numpy.pi / 180.0

''' Columns of unit_cell are the a, b, and c cell vectors. '''
unit_cell = numpy.array([[6.63,  0.00, 6.50*numpy.cos(beta)],
                         [0.00, 20.78, 0.0],
                         [0.00,  0.00, 6.50*numpy.sin(beta)]])

C = numpy.vstack([numpy.dot(C, unit_cell) for C in apply_space_group(C1)])

''' Creating a dictionary for existing atoms '''
atom_types_dict = {'c1': -4, 'c2': -1, 'c3': -2, 'h': -3} 

def sign(x):
    ''' function to change the sign of coordinates '''
    if x < 0:
        a = -1
    elif x > 0:
        a = 1
    return a

def unit_vector(v):
    ''' creating a unit normal vector of vector v'''
    return v / numpy.linalg.norm(v)

H = numpy.zeros([72,3])
for i in range(4):
    H[18*i:18*(i+1),:] = applyH.apply_hydrogens(C[9*i:9*(i+1),:],i)

if key == 1:
    plt.plot_unit_cell(C,H)

def crystal(a,b,c):
    ''' Creating the initial crystalline system - a & b & c are the dimensions of the system '''    
    C = numpy.vstack([C for C in apply_space_group(C1)])    
    Carbons = numpy.zeros([a*b*c*36,3])    
    Hydrogens = numpy.zeros([a*b*c*72,3])
                
    for i in range(a):
        Carbons[(36*i):(36*(i+1)),0] = C[:,0] + float(i)
        Carbons[(36*i):(36*(i+1)),1] = C[:,1]
        Carbons[(36*i):(36*(i+1)),2] = C[:,2]   
    
    for j in range(1,b):
        Carbons[36*(a*j):36*(a*(j+1)),0] = Carbons[0:36*a,0]
        Carbons[36*(a*j):36*(a*(j+1)),1] = Carbons[0:36*a,1] + float(j)         
        Carbons[36*(a*j):36*(a*(j+1)),2] = Carbons[0:36*a,2] 
        
    for k in range(1,c):
        Carbons[36*(a*b*k):36*(a*b*(k+1)),0] = Carbons[0:36*a*b,0]
        Carbons[36*(a*b*k):36*(a*b*(k+1)),1] = Carbons[0:36*a*b,1]          
        Carbons[36*(a*b*k):36*(a*b*(k+1)),2] = Carbons[0:36*a*b,2] + float(k)                    
    
    Carbons = numpy.dot(Carbons, unit_cell) 
    
    # number of chains = a*b*c*4
    for i in range(a*b*c*4):        
        Hydrogens[18*i:18*(i+1),:] = applyH.apply_hydrogens(Carbons[9*i:9*(i+1),:],i)            
        
    return Carbons, Hydrogens

C, H = crystal(a,b,c)

if key == 1:
    plt.plot_crystal(C, H)

# Determine bonds
CC_bonds = bonds.carbon_bonds(C)
CH_bonds = bonds.hydrogen_bonds(C,H)
bonds = numpy.append(CC_bonds,CH_bonds, axis = 0)

# Determine the angle between atoms
total_angles = ang.angles(CC_bonds, CH_bonds)


# Determine the dihedral angles in the molecule
# Proper dihedrdal angle
proper_dihedrals = dihed.proper_dihedral(CC_bonds, CH_bonds)

# Improper dihedral angle
improper_dihedrals = dihed.improper_dihedral(CC_bonds, CH_bonds)


out.outputfile(filename,unit_cell,C,H,CC_bonds,CH_bonds,atom_types_dict,total_angles,proper_dihedrals,improper_dihedrals,a,b,c)

