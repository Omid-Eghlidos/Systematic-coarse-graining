#!usr/bin/env python3

import numpy
import math
from matplotlib import pyplot

# Plotting one unit cell
def plot_unit_cell(C,H):    
    ''' this function only plots the carbons and hydrogens of one unitcell 
        in xy, xz, and yz planes'''   
    # plot in x-y plane 
    figure = pyplot.figure(figsize=(9,5))
    aC = pyplot.subplot(221)
    aH = pyplot.subplot(221)   
    # plotting Carbons
    aC.plot(C[:9,1],    -C[:9,0],    '.', ms=30, alpha=0.5, label='C1')
    aC.plot(C[9:18,1],  -C[9:18,0],  '.', ms=30, alpha=0.5, label='C2')
    aC.plot(C[18:27,1], -C[18:27,0], '.', ms=30, alpha=0.5, label='C3')
    aC.plot(C[27:36,1], -C[27:36,0], '.', ms=30, alpha=0.5, label='C4')
    aC.set_xlabel('y')
    aC.set_ylabel('-x')
    aC.set_aspect('equal')
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()                
    # plotting Hydrogens
    aH.plot(H[:18,1],    -H[:18,0],    '.', ms=20, alpha=0.5, label='H1')
    aH.plot(H[18:36,1],  -H[18:36,0],  '.', ms=20, alpha=0.5, label='H2')
    aH.plot(H[36:54,1], -H[36:54,0], '.', ms=20, alpha=0.5, label='H3')
    aH.plot(H[54:72,1], -H[54:72,0], '.', ms=20, alpha=0.5, label='H4')
    aH.set_xlabel('y')
    aH.set_ylabel('-x')
    aH.set_aspect('equal')    
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()        
       
    # plot in x-z plane 
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    aC = pyplot.subplot(222)
    aH = pyplot.subplot(222)
    # plotting Carbons
    aC.plot(C[:9,2],    -C[:9,0],    '.', ms=30, alpha=0.5, label='C1')
    aC.plot(C[9:18,2],  -C[9:18,0],  '.', ms=30, alpha=0.5, label='C2')
    aC.plot(C[18:27,2], -C[18:27,0], '.', ms=30, alpha=0.5, label='C3')
    aC.plot(C[27:36,2], -C[27:36,0], '.', ms=30, alpha=0.5, label='C4')
    aC.set_xlabel('z')
    aC.set_ylabel('-x')
    aC.set_aspect('equal')
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()            
    # plotting Hydrogens
    aH.plot(H[:18,2],    -H[:18,0],    '.', ms=20, alpha=0.5, label='H1')
    aH.plot(H[18:36,2],  -H[18:36,0],  '.', ms=20, alpha=0.5, label='H2')
    aH.plot(H[36:54,2], -H[36:54,0], '.', ms=20, alpha=0.5, label='H3')
    aH.plot(H[54:72,2], -H[54:72,0], '.', ms=20, alpha=0.5, label='H4')
    aH.set_xlabel('y')
    aH.set_ylabel('-x')
    aH.set_aspect('equal')   
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()        
       
    # plot in y-z plane
    aC = pyplot.subplot(223)
    aH = pyplot.subplot(223)    
    # plotting Carbons
    aC.plot(C[:9,1],    C[:9,2],    '.', ms=30, alpha=0.5, label='C1')
    aC.plot(C[9:18,1],  C[9:18,2],  '.', ms=30, alpha=0.5, label='C2')
    aC.plot(C[18:27,1], C[18:27,2], '.', ms=30, alpha=0.5, label='C3')
    aC.plot(C[27:36,1], C[27:36,2], '.', ms=30, alpha=0.5, label='C4')
    aC.set_xlabel('y')
    aC.set_ylabel('z')
    aC.set_aspect('equal')
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()        
    # plotting Hydrogens
    aH.plot(H[:18,1],    H[:18,2],    '.', ms=20, alpha=0.5, label='H1')
    aH.plot(H[18:36,1],  H[18:36,2],  '.', ms=20, alpha=0.5, label='H2')
    aH.plot(H[36:54,1], H[36:54,2], '.', ms=20, alpha=0.5, label='H3')
    aH.plot(H[54:72,1], H[54:72,2], '.', ms=20, alpha=0.5, label='H4')
    aH.set_xlabel('y')
    aH.set_ylabel('-x')
    aH.set_aspect('equal')    
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()        
   
# Plotting the crystalline system
def plot_crystal(C, H):
    ''' This function plots the whole created crystaline system in 
        xy, xz, and yz planes'''
    # plot in x-y plane
    fig1 = pyplot.figure(figsize = (10,10))
    aC1 = pyplot.plot(C[:,1], -C[:,0],'.', ms=20, alpha=0.5, label='Carbons')
    aH1 = pyplot.plot(H[:,1], -H[:,0],'.', ms=15, alpha=0.5, label='Hydrogens')
    pyplot.xlabel('y')
    pyplot.ylabel('x')
    pyplot.axis('equal')
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()
    # plot in y-z plane
    fig2 = pyplot.figure(figsize = (10,10))
    aC2 = pyplot.plot(C[:,1], -C[:,2],'.', ms=20, alpha=0.5, label='Carbons')
    aH2 = pyplot.plot(H[:,1], -H[:,2],'.', ms=15, alpha=0.5, label='Hydrogens')
    pyplot.xlabel('y')
    pyplot.ylabel('z')
    pyplot.axis('equal')
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.show()
    pyplot.tight_layout()
