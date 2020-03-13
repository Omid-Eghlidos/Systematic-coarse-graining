#!/usr/bin/env python3

import numpy

# Function to determine the proper dihedral angle in the molecule
def proper_dihedral(CC_bonds, CH_bonds):
    ''' this function find the dihedral bonding in our iPP system '''
    dihedral_angle = []
    # Adding type 1 dihedrals which only conatins carbons
    for i,b in enumerate(CC_bonds):
        atomsC1 = [int(a) for a in b if a!= -11]
	# First centering carbon bonded to one other carbon 
        if len(atomsC1) == 1:
            continue
	# First centering carbon bonded to two other carbons  
        elif len(atomsC1) == 2:
            for j,d in enumerate(CC_bonds):
                atomsC2 = [int(c) for c in d if c != -11]
                # if the carbon index is equal to any of the 
                # bondings of the first carbon it means its mutual 
                # bond and is already counted 
                if j == atomsC1[1]:
		    # The second centering carbon bonded to one other carbon  
                    if len(atomsC2) == 1:
                        continue                        
	            # The second centering carbon bonded to two other carbons 
                    elif len(atomsC2) == 2:
                        if atomsC2[0] == i:
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[1], atomsC2[1], 1]))
                        else:
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[1], atomsC2[0], 1]))
		    # The third centering carbon bonded to three other carbons
                    elif len(atomsC2) == 3:
                        # if the first central carbon is equal to any of the 
                        # bondings of the second carbon it means its mutual 
                        # bond and is already counted 
                        if atomsC2[0] == i:
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[1], atomsC2[1], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[1], atomsC2[2], 1]))
                        elif atomsC2[1] == i:
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[1], atomsC2[0], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[1], atomsC2[2], 1]))
	# First centering carbon bonded to three other carbons 
        elif len(atomsC1) == 3:
            for j,d in enumerate(CC_bonds):
                atomsC2 = [int(c) for c in d if c != -11]
                if j == atomsC1[2]:
                    # The second centering atom bonded to one other carbon
                    if len(atomsC2) == 1:
                        continue
                    # The second centering atom bonded to two other carbons
                    elif len(atomsC2) == 2:
                        if atomsC2[0] == i:
                            # atom id 1, atom id 2,
                            # atom id 3, atom id 4, dihedral_type 
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
                                                  atomsC1[2], atomsC2[1], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[1], i,
						  atomsC1[2], atomsC2[1], 1]))
                        else:
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[2], atomsC2[0], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[1], i,
						  atomsC1[2], atomsC2[1], 1]))                   
                    # The second centering atom bonded to three other carbons
                    elif len(atomsC2) == 3:
                        if atomsC2[0] == i:
                            # atom id 1, atom id 2,
                            # atom id 3, atom id 4, dihedral_type 
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[2], atomsC2[1], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[2], atomsC2[2], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[1], i,
						  atomsC1[2], atomsC2[1], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[1], i,
						  atomsC1[2], atomsC2[2], 1]))
                        elif atomsC2[1] == i:
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[2], atomsC2[0], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
						  atomsC1[2], atomsC2[2], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[1], i,
						  atomsC1[2], atomsC2[0], 1]))
                            dihedral_angle.append(numpy.array([atomsC1[1], i,
						  atomsC1[2], atomsC2[2], 1]))
    
    # Adding type 2 dihedral angles with two centering carbons and
    # one carbon bonded to one and one hydrogen bonded to the other one    
    for i,b in enumerate(CC_bonds):
        atomsC1 = [int(a) for a in b if a != -11]
        # First centering carbon bonded to one other carbon
        if len(atomsC1) == 1:
            for j,d in enumerate(CH_bonds):
                atomsH1 = [int(c) for c in d if c != -11]
                for l,h in enumerate(CC_bonds):
                    atomsC2 = [int(g) for g in h if g != -11]
                    if i == j and l == atomsC1[0]:
                        for m in range(len(atomsH1)):
                            for n in range(len(atomsC2)):
                                if atomsC2[n] != i and atomsC1[0] > i:
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i,
                                                           atomsC1[0],
                                                           atomsC2[n],
                                                           2]))           
        # First centering carbon bonded to two other carbons
        elif len(atomsC1) == 2:             
            for j,d in enumerate(CH_bonds):
                atomsH1 = [int(c) for c in d if c != -11]
                for l,h in enumerate(CC_bonds):
                    atomsC2 = [int(g) for g in h if g != -11]
                    if i == j and l == atomsC1[1]:
                        for m in range(len(atomsH1)):
                            for n in range(len(atomsC2)):
                                if atomsC2[n] != i and atomsC1[1] > i:
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i,
                                                           atomsC1[1],
                                                           atomsC2[n],
                                                           2]))
                    elif j == atomsC1[1]  and l == atomsC1[0]:
                        for m in range(len(atomsH1)):
                            dihedral_angle.append(numpy.array([atomsC1[0], i,
                                                               atomsC1[1], 
                                                               atomsH1[m]-1,
                                                               2]))
        elif len(atomsC1) == 3:
            for j,d in enumerate(CH_bonds):
                atomsH1 = [int(c) for c in d if c != -11]
                for l,h in enumerate(CC_bonds):
                    atomsC2 = [int(g) for g in h if g != -11]
                    if i == j and l == atomsC1[2]:
                        for m in range(len(atomsH1)):
                            for n in range(len(atomsC2)):
                                if atomsC2[n] != i and atomsC1[2] > i:
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i,
                                                           atomsC1[2],
                                                           atomsC2[n],
                                                           2]))
                    elif l == i and j == atomsC1[2]:
                        for m in range(len(atomsC2)):
                            for n in range(len(atomsH1)):
                                if atomsC2[n] < i and atomsC1[2] != atomsC2[m]:
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsC2[m], i,
                                                           atomsC1[2],
                                                           atomsH1[n]-1,
                                                           2]))

    # Adding type 3 dihedral angles with two centering carbons and
    # one hydrogen bonded to one and one hydeogen bonded to the other one  
    for i,b in enumerate(CC_bonds):
        atomsC = [int(a) for a in b if a != -11]
        for j,d in enumerate(CH_bonds):
            atomsH1 = [int(c) for c in d if c != -11]
            # First centering carbon bonded to one other carbon
            if len(atomsC) == 1:
                for l,h in enumerate(CH_bonds):
                    atomsH2 = [int(g) for g in h if g != -11]
                    if i == j and atomsC[0] == l and l > i:
                        for m in range(len(atomsH1)):
                            for n in range(len(atomsH2)):
                                # atom id 1, atom id 2,
                                # atom id 3, atom id 4, dihedral_type 
                                dihedral_angle.append(numpy.array([atomsH1[m]-1,
    							           i, atomsC[0], 
							           atomsH2[n]-1,
							           3]))
            # First centering carbon bonded to two other carbons 
            elif len(atomsC) == 2:
                for l,h in enumerate(CH_bonds):
                    atomsH2 = [int(g) for g in h if g != -11]
                    if i == j:
                        if atomsC[0] == l and l > j:
                            for m in range(len(atomsH1)):
                                for n in range(len(atomsH2)):
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i, 
                                                           atomsC[0], 
                                                           atomsH2[n]-1,
                                                           3]))                        
                        elif atomsC[1] == l and l > j:
                            for m in range(len(atomsH1)):
                                for n in range(len(atomsH2)):
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i,
                                                           atomsC[1],
                                                           atomsH2[n]-1,
			    				   3]))                        
            # First centering carbon bonded to three other carbons 
            elif len(atomsC) == 3:
                for l,h in enumerate(CH_bonds):
                    atomsH2 = [int(g) for g in h if g != -11]
                    if i == j:                        
                        if atomsC[0] == l and l > j:
                            for m in range(len(atomsH1)):
                                for n in range(len(atomsH2)):
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i, 
                                                           atomsC[0], 
                                                           atomsH2[n]-1,
                                                           3]))                                 
                        elif atomsC[1] == l and l > j:
                            for m in range(len(atomsH1)):
                                for n in range(len(atomsH2)):
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i,
                                                           atomsC[1],
                                                           atomsH2[n]-1,
			    				   3]))
                        elif atomsC[2] == l and l > j:
                            for m in range(len(atomsH1)):
                                for n in range(len(atomsH2)):
                                    # atom id 1, atom id 2,
                                    # atom id 3, atom id 4, dihedral_type 
                                    dihedral_angle.append(numpy.array(
                                                          [atomsH1[m]-1, i,
                                                           atomsC[2],
                                                           atomsH2[n]-1,
			    				   3]))

    return dihedral_angle 


# Function to determine the improper dihedral angle in the molecule
def improper_dihedral(CC_bonds, CH_bonds):
    improper_angle = []


    return improper_angle



