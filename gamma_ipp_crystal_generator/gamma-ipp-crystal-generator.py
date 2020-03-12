#!/usr/bin/env python3

import numpy
import math
from matplotlib import pyplot

''' Energy Calculations for Isotactic Polypropylene: 
A Comparison between Models of the alpha and Gamma 
Crystalline Structures from D.R. Ferro 1992'''


C1 = numpy.array([[-0.1391, -0.0921, 0.0452],
                  [-0.1098, -0.0289, 0.0779],
                  [-0.2625, 0.0248 , 0.0925],
                  [0.0121 , 0.0860 , 0.0763],
                  [0.1769 , 0.0487 , 0.0641],
                  [0.2501 , -0.0646, 0.0840],
                  [0.2886 , 0.1711 , 0.0638],
                  [0.2397 , 0.2927 , 0.0434],
                  [0.2156 , 0.2524 , 0.0087]])


''' Angle between c vector and x-axis. '''
beta = 90.0 * numpy.pi / 180.0

''' Columns of unit_cell are the a, b, and c cell vectors. '''
unit_cell = numpy.array([[8.54,  0.00, 0.00],
                         [0.00,  9.93, 0.00],
                         [0.00,  0.00, 42.41]])


def apply_space_group(C1):
    ''' Edited symmetry operations to make figure match - not checked '''
    C2, C3, C4 = C1.copy(), C1.copy(), C1.copy()

    C2[:,0] = -C2[:,0]
    C2[:,1] =  C2[:,1]
    C2[:,2] =  C2[:,2] 

    C3[:,0] =  C3[:,0] - 0.5 
    C3[:,1] =  C3[:,1] 
    C3[:,2] =  C3[:,2] 

    C4[:,0] = -C3[:,0] 
    C4[:,1] =  C3[:,1] 
    C4[:,2] =  C3[:,2]    
    
    return C1, C2, C3, C4

def apply_hydrogens(C,i):
    ''' Construction of hydrogens from Theodorou (MD of aPP Melts) with 
    l_H=0.11 nm and theta_H=1.28 rad and c=+1 in C1 & C4 and c=-1 in C2 & C3'''   
    
    temp = numpy.zeros(3)
    H = numpy.zeros([9,3])
    l_H = 1.1
    theta_H = 1.28
    c = 1
    
    
    # Non-Methyl Hydrogens (Chiral and Gemini)
    # Chiral(Pendant) hydrogen
    # r_Hi = r_Ci+c*l_H*[((r_Ri)-(r_Ci-1))cross((r_Ri)-(r_Ci+1))/|((r_Ri)-(r_Ci-1))cross((r_Ri)-(r_Ci+1))|]
    
    def chiral_H_direction(Ci,Cj,Ck):      
        return unit_vector(numpy.cross(Ci-Cj, Ci-Ck))
    
    temp[0] = C[8,0]+1.07  # deducing along the unit cell c axis and transform it into orthonormal coordinates
    temp[1] = C[8,1]   
    temp[2] = C[8,2]-6.41
    H[0,:] = C[1,:] + c*l_H*chiral_H_direction(C[0,:],temp[:],C[2,:])
    H[3,:] = C[4,:] + c*l_H*chiral_H_direction(C[3,:],C[2,:],C[5,:])
    H[6,:] = C[7,:] + c*l_H*chiral_H_direction(C[6,:],C[5,:],C[8,:])
    
     
    # Gemini hydrogens 
    # b(i) = rc(i) - rc(i-1)/|rc(i)-rc(i-1)|
    b = [None]*6
        
    b[0] = unit_vector(C[2,:] - C[1,:])
    b[1] = unit_vector(C[4,:] - C[2,:])
    
    b[2] = unit_vector(C[5,:] - C[4,:])
    b[3] = unit_vector(C[7,:] - C[5,:])
    
    b[4] = unit_vector(C[8,:] - C[7,:])
    temp[0] = C[1,0]+1.07
    temp[1] = C[1,1]
    temp[2] = C[1,2]-6.41
    b[5] = unit_vector(temp[:] - C[8,:])


    # u(i) = (b(i) - b(i+1))/sqrt(2*(1-b(i).b(i+1))    
    u = [None]*3

    def bisector_vector(b1,b2):
        return (b1 - b2) / math.sqrt(2.0*(1.0 - numpy.dot(b1,b2)))
      
    # v(i) = (b(i) X b(i+1))/|b(i) X b(i+1)|
    v = [None]*3

    def plane_normal(b1,b2):
        return numpy.cross(b1,b2)/numpy.linalg.norm(numpy.cross(b1,b2))

    # rH(i) = rc(i) + lH*(sin(thetaH/2)*u(i) +- cos(thetaH/2)*v(i))
    def gemini_hydrogen(u,v,n):
        if n == 1:
            return (numpy.sin(theta_H/2)*u+numpy.cos(theta_H/2)*v)
        else:
            return (numpy.sin(theta_H/2)*u-numpy.cos(theta_H/2)*v)
    
    v[0] = plane_normal(b[0],b[1])
    u[0] = bisector_vector(b[0],b[1])
    H[1,:] = C[2,:]+l_H*gemini_hydrogen(u[0],v[0],1)
    H[2,:] = C[2,:]+l_H*gemini_hydrogen(u[0],v[0],2)

    v[1] = plane_normal(b[2],b[3])
    u[1] = bisector_vector(b[2],b[3])
    H[4,:] = C[5,:]+l_H*gemini_hydrogen(u[1],v[1],1)
    H[5,:] = C[5,:]+l_H*gemini_hydrogen(u[1],v[1],2)

    v[2] = plane_normal(b[4],b[5])
    u[2] = bisector_vector(b[4],b[5])
    H[7,:] = C[8,:]+l_H*gemini_hydrogen(u[2],v[2],1)
    H[8,:] = C[8,:]+l_H*gemini_hydrogen(u[2],v[2],2)  
    
    '''
    # Methyl Group Hydrogens 
    for i in range(0,9,3):
        H[9+i,0] = C[1+i,0] + 1.037
        H[9+i,1] = C[1+i,1] 
        H[9+i,2] = C[1+i,2] - 0.37
        H[10+i,0] = C[1+i,0] - 0.5185        
        H[10+i,1] = C[1+i,1] + 0.898
        H[10+i,2] = C[1+i,2] - 0.37
        H[11+i,0] = C[1+i,0] - 0.5185
        H[11+i,1] = C[1+i,1] - 0.898       
        H[11+i,2] = C[1+i,2] - 0.37
    '''
    
    return H

def plot_unit_cell(CC):
    
    C = numpy.vstack([numpy.dot(C, unit_cell) for C in apply_space_group(CC)])    
    H = numpy.zeros([36,3])
    
    for i in range(4):
        H[9*i:9*(i+1),:] = apply_hydrogens(C[9*i:9*(i+1),:],i) 
    
    lx, ly = unit_cell[0,0], unit_cell[1,1]
    
    figure = pyplot.figure(figsize=(5,5))
    aC = pyplot.subplot(221)
    aH = pyplot.subplot(221)
    
    aC.plot(C[:9,1],    -C[:9,0],    '.', ms=30, alpha=0.5, label='C1')
    '''aC.plot(C[9:18,1],  -C[9:18,0],  '.', ms=30, alpha=0.5, label='C2')
    aC.plot(C[18:27,1], -C[18:27,0], '.', ms=30, alpha=0.5, label='C3')
    aC.plot(C[27:36,1], -C[27:36,0], '.', ms=30, alpha=0.5, label='C4')
    aC.set_xlabel('y')
    aC.set_ylabel('-x')
    aC.set_aspect('equal')
    
    aH.plot(H[:9,1],    -H[:9,0],    '.', ms=20, alpha=0.5, label='H1')
    aH.plot(H[9:18,1],  -H[9:18,0],  '.', ms=20, alpha=0.5, label='H2')
    aH.plot(H[18:27,1], -H[18:27,0], '.', ms=20, alpha=0.5, label='H3')
    aH.plot(H[27:36,1], -H[27:36,0], '.', ms=20, alpha=0.5, label='H4')
    aH.set_xlabel('y')
    aH.set_ylabel('-x')
    aH.set_aspect('equal')
    '''
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    aC = pyplot.subplot(222)
    aH = pyplot.subplot(222)
    
    aC.plot(C[:9,2],    -C[:9,0],    '.', ms=30, alpha=0.5, label='C1')
    aC.plot(C[9:18,2],  -C[9:18,0],  '.', ms=30, alpha=0.5, label='C2')
    aC.plot(C[18:27,2], -C[18:27,0], '.', ms=30, alpha=0.5, label='C3')
    aC.plot(C[27:36,2], -C[27:36,0], '.', ms=30, alpha=0.5, label='C4')
    aC.set_xlabel('z')
    aC.set_ylabel('-x')
    aC.set_aspect('equal')
    
    aH.plot(H[:9,2],    -H[:9,0],    '.', ms=20, alpha=0.5, label='H1')
    aH.plot(H[9:18,2],  -H[9:18,0],  '.', ms=20, alpha=0.5, label='H2')
    aH.plot(H[18:27,2], -H[18:27,0], '.', ms=20, alpha=0.5, label='H3')
    aH.plot(H[27:36,2], -H[27:36,0], '.', ms=20, alpha=0.5, label='H4')
    aH.set_xlabel('y')
    aH.set_ylabel('-x')
    aH.set_aspect('equal')

    aC = pyplot.subplot(223)
    aH = pyplot.subplot(223)

    aC.plot(C[:9,1],    C[:9,2],    '.', ms=30, alpha=0.5, label='C1')
    aC.plot(C[9:18,1],  C[9:18,2],  '.', ms=30, alpha=0.5, label='C2')
    aC.plot(C[18:27,1], C[18:27,2], '.', ms=30, alpha=0.5, label='C3')
    aC.plot(C[27:36,1], C[27:36,2], '.', ms=30, alpha=0.5, label='C4')
    aC.set_xlabel('y')
    aC.set_ylabel('z')
    aC.set_aspect('equal')
    
    aH.plot(H[:9,1],    H[:9,2],    '.', ms=20, alpha=0.5, label='H1')
    aH.plot(H[9:18,1],  H[9:18,2],  '.', ms=20, alpha=0.5, label='H2')
    aH.plot(H[18:27,1], H[18:27,2], '.', ms=20, alpha=0.5, label='H3')
    aH.plot(H[27:36,1], H[27:36,2], '.', ms=20, alpha=0.5, label='H4')
    aH.set_xlabel('y')
    aH.set_ylabel('-x')
    aH.set_aspect('equal')
    
    pyplot.tight_layout()
    
    return CC
    

def unit_vector(v):
    return v / numpy.linalg.norm(v)

        
def crystal(a,b,c):
    ''' Creating the initial crystalline system - a & b & c are the dimensions of the system '''   
           
    C = numpy.vstack([C for C in apply_space_group(C1)])    
    Carbons = numpy.zeros([a*b*c*36,3])    
    Hydrogens = numpy.zeros([a*b*c*36,3])
                
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
        Hydrogens[9*i:9*(i+1),:] = apply_hydrogens(Carbons[9*i:9*(i+1),:],i)            
    
    
    
    return Carbons, Hydrogens

def plot_crystal(C, H):
    
    fig1 = pyplot.figure(figsize = (10,5))      
    aC1 = pyplot.plot(C[:,1], -C[:,0],'.', ms=30, alpha=0.5, label='Carbons')
    aH1 = pyplot.plot(H[:,1], -H[:,0],'.', ms=20, alpha=0.5, label='Hydrogens')
    pyplot.xlabel('y')
    pyplot.ylabel('x')   
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    
    fig2 = pyplot.figure(figsize = (10,5))    
    aC2 = pyplot.plot(C[:,1], -C[:,2],'.', ms=30, alpha=0.5, label='Carbons')
    aH2 = pyplot.plot(H[:,1], -H[:,2],'.', ms=20, alpha=0.5, label='Hydrogens')
    pyplot.xlabel('y')
    pyplot.ylabel('z')    
    pyplot.legend(ncol=4, bbox_to_anchor=(0,1), loc='lower left')
    pyplot.tight_layout()

CC = plot_unit_cell(C1)
C, H = crystal(1,1,1)

plot_crystal(C, H)

