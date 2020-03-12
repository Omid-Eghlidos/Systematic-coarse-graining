#!/usr/bin/env python3

import numpy

# Dteremine the angle between atoms
def angles(CC_bonds, CH_bonds):
    #print(CC_bonds)
    #print(CH_bonds)
    total_angles = []
    for i,b in enumerate(CC_bonds):
        atomsC = [int(a) for a in b if a != -11]
        if len(atomsC) < 2:
            continue
        elif len(atomsC) == 2:
            total_angles.append(numpy.array([atomsC[0], i, atomsC[1], 1]))
        else:
            total_angles.append(numpy.array([atomsC[0], i, atomsC[1], 1]))
            total_angles.append(numpy.array([atomsC[0], i, atomsC[2], 1]))
            total_angles.append(numpy.array([atomsC[1], i, atomsC[2], 1]))


    for k,l in enumerate(CC_bonds):
        atomsC = [int(e) for e in l if e != -11]
        for m,n in enumerate(CH_bonds):
            atomsH = [int(f) for f in n if f != -11]
            if len(atomsC) < 2:
                if k == m:
                    if len(atomsH) < 2:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                    elif len(atomsH) == 2:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[1]-1, 2]))
                    else:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[1]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[2]-1, 2]))
            if len(atomsC) == 2:
                if k == m:                     
                    if len(atomsH) < 2:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[0]-1, 2]))
                    elif len(atomsH) == 2:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[1]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[1]-1, 2]))
                    else:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[1]-1, 2]))             
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[1]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[2]-1, 2]))     
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[2]-1, 2]))
            if len(atomsC) > 2:
                if k == m:                    
                    if len(atomsH) < 2:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[2], k, atomsH[0]-1, 2]))
                    elif len(atomsH) == 2:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[2], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[1]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[1]-1, 2]))
                        total_angles.append(numpy.array([atomsC[2], k, atomsH[0]-1, 2]))
                    else:
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[2], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[1]-1, 2]))             
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[1]-1, 2]))
                        total_angles.append(numpy.array([atomsC[2], k, atomsH[0]-1, 2]))
                        total_angles.append(numpy.array([atomsC[0], k, atomsH[2]-1, 2]))     
                        total_angles.append(numpy.array([atomsC[1], k, atomsH[2]-1, 2]))
                        total_angles.append(numpy.array([atomsC[2], k, atomsH[0]-1, 2]))
    
    for j,d in enumerate(CH_bonds):
        atoms = [int(c) for c in d if c != -11]
        if len(atoms) < 2:
            continue
        if len(atoms) == 2:
            total_angles.append(numpy.array([atoms[0]-1, j, atoms[1]-1, 3]))
        else:
            total_angles.append(numpy.array([atoms[0]-1, j, atoms[1]-1, 3]))
            total_angles.append(numpy.array([atoms[0]-1, j, atoms[2]-1, 3]))    
            total_angles.append(numpy.array([atoms[1]-1, j, atoms[2]-1, 3]))    

    return total_angles

