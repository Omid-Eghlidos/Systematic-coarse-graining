#!/usr/bin/env python3

import numpy

# Function to determine the proper dihedral angle in the molecule
def proper_dihedral(CC_bonds, CH_bonds):
    proper_dihedral_angle = []
    for i,b in enumerate(CC_bonds):
        atomsC1 = [int(a) for a in b if a!= -11]
        if len(atomsC1) < 2:
            continue
        elif len(atomsC1) == 2:
            for j,d in enumerate(CC_bonds):
                atomsC2 = [int(c) for c in d if c != -11]
                if j != i:
                    if j == atomsC1[0] and j < i:
                        if len(atomsC2) < 2:
                            proper_dihedral_angle.append(numpy.array([atomsC2[0], atomsC1[0], i, atomsC1[1], 1]))
                        elif len(atomsC2) == 2:
                            proper_dihedral_angle.append(numpy.array([atomsC2[0], atomsC1[0], i, atomsC1[1], 1]))
                            proper_dihedral_angle.append(numpy.array([atomsC2[1], atomsC1[0], i, atomsC1[1], 1]))
                        else:
                            proper_dihedral_angle.append(numpy.array([atomsC2[0], atomsC1[0], i, atomsC1[1], 1]))
                            proper_dihedral_angle.append(numpy.array([atomsC2[1], atomsC1[0], i, atomsC1[1], 1]))
                            proper_dihedral_angle.append(numpy.array([atomsC2[2], atomsC1[0], i, atomsC1[1], 1]))
                    elif j == atomsC1[1] and j > i:
                        if len(atomsC2) < 2:
                            proper_dihedral_angle.append(numpy.array([atomsC1[0], i, atomsC1[1], atomsC2[0], 1]))
                        elif len(atomsC2) == 2:
                            proper_dihedral_angle.append(numpy.array([atomsC1[0], i, atomsC1[1], atomsC2[0], 1]))
                            proper_dihedral_angle.append(numpy.array([atomsC1[0], i, atomsC1[1], atomsC2[1], 1]))
                        else:
                            proper_dihedral_angle.append(numpy.array([atomsC1[0], i, atomsC1[1], atomsC2[0], 1]))
                            proper_dihedral_angle.append(numpy.array([atomsC1[0], i, atomsC1[1], atomsC2[1], 1]))
                            proper_dihedral_angle.append(numpy.array([atomsC1[0], i, atomsC1[1], atomsC2[2], 1]))

    print(proper_dihedral_angle)

    return proper_dihedral_angle



'''                
    for i,b in enumerate(CC_bonds):
        atomsC = [int(a) for a in b if a != -11]
        for j,d in enumerate(CH_bonds):
            atomsH1 = [int(c) for c in d if c != -11]           
            if i == j:
                for l,h in enumerate (CH_bonds):
                    atomsH2 = [int(g) for g in h if g != -11]
                    if len(atomsC) < 2:
                        if l == atomsC1[0]:
                            for m in range(len(atomsH1)):
                                for n in range(len(atomsH2)):
                                    # atom id 1, atom id 2, atom id 3, atom id 4, dihedral_type 
                                    peoper_dihedral_angle.append(numpy.array([atomsH1[m], i, atomsC[0], atomsH2[n], 3]))
                    elif len(atomsC) == 2:
                        if l == atomsC[0]:
                            for m in range(len(atomsH1)):
                                for n in range(len(atomsH2)):
                                    # atom id 1, atom id 2, atom id 3, atom id 4, dihedral_type 
                                    peoper_dihedral_angle.append(numpy.array([atomsH1[m], i, atomsC[0], atomsH2[n], 3]))



                            
                            




        if i == j:
            for k,f in enumerate(CC_bonds):
                atomsC2 = [int(e) for e in f if e != -11]
                    if k == l:


                
                
                


    return proper_dihedral_angle
'''



# Function to determine the improper dihedral angle in the molecule
def improper_dihedral(CC_bonds, CH_bonds):
    improper_dihedral_angle = []


    return improper_dihedral_angle



