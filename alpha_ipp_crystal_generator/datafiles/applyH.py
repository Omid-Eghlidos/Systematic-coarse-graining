#!/usr/bin/env python3

import numpy
import math

def sign(x):
    ''' a function to change the sign '''
    if x < 0:
        a = -1
    elif x > 0:
        a = 1
    return a

def unit_vector(v):
    ''' creating a unit normal vector '''
    return v / numpy.linalg.norm(v)

def apply_hydrogens(C,i):
    ''' Construction of hydrogens from Theodorou (MD of aPP Melts) with
    l_H=0.11 nm and theta_H=1.28 rad and c=+1 for iPP - all the following
    formulas are taken from this paper '''

    temp = numpy.zeros(3)
    H = numpy.zeros([18,3])
    l_H = 1.1
    theta_H = 1.28
    c = 1

    # Non-Methyl Hydrogens (Chiral and Gemini)
    # Chiral(Pendant) hydrogen
    def chiral_H_direction(Ci,Cj,Ck):
        return unit_vector(numpy.cross(Ci-Cj, Ci-Ck))

    def __chiral__():
        ''' r_Hi = r_Ci+c*l_H*[((r_Ri)-(r_Ci-1))cross((r_Ri)-(r_Ci+1))/
                   |((r_Ri)-(r_Ci-1))cross((r_Ri)-(r_Ci+1))|]'''

        # deducting along the unit cell c axis and transform it 
        # into orthonormal coordinates
        temp[0] = C[8,0]+1.07
        temp[1] = C[8,1]
        temp[2] = C[8,2]-6.41
        H[0,:] = C[1,:] + c*l_H*chiral_H_direction(C[0,:],temp[:],C[2,:])
        H[3,:] = C[4,:] + c*l_H*chiral_H_direction(C[3,:],C[2,:],C[5,:])
        H[6,:] = C[7,:] + c*l_H*chiral_H_direction(C[6,:],C[5,:],C[8,:])

    # Gemini hydrogens
    def __gemini__():
        ''' adding the gemini hydrogens on C2 '''
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

    # Methyl Group Hydrogens
    def __methyl__():
        ''' Adding methyl hydrogens to the C3 - these hydrone coordinates 
            are derived based on geometric calculation with respect to the 
            carbon ther are connectde to in 3D space'''
        for i in range(0,9,3):
            H[9+i,0] = C[i,0] + sign(C[1+i,0])*1.037
            H[9+i,1] = C[i,1] + 0.0
            H[9+i,2] = C[i,2] - sign(C[1+i,2])*0.37
            H[10+i,0] = C[i,0] - sign(C[1+i,0])*0.5185
            H[10+i,1] = C[i,1] + sign(C[1+i,1])*0.898
            H[10+i,2] = C[i,2] - sign(C[1+i,2])*0.37
            H[11+i,0] = C[i,0] - sign(C[1+i,0])*0.5185
            H[11+i,1] = C[i,1] - sign(C[1+i,1])*0.898
            H[11+i,2] = C[i,2] - sign(C[1+i,2])*0.37

    __chiral__()
    __gemini__()
    __methyl__()

    return H




