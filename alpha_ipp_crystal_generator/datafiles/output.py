#!/usr/bin/env python3

import numpy
import math
import angles as ang
import dihedrals as dihed

''' define the potential energy '''
import pcff as potential

def sign(x):
    ''' a function to change the sign '''
    if x < 0:
        a = -1
    elif x > 0:
        a = 1
    return a

def outputfile(filename,unit_cell,C,H,CC_bonds,CH_bonds,atom_types_dict,
               total_angles,proper_dihedrals,improper_dihedrals,a,b,c, qkey):
    ''' Create a LAMMPS data file as an output '''    
    f = open(filename,'w+')

    # initialization - the size of the system (atoms, bonds, angles, etc.)
    atoms = numpy.append(C,H, axis=0)
    bonds = numpy.append(CC_bonds,CH_bonds, axis=0)        
    atom_types = len(potential.paircoef)
    bond_types = len(potential.bondcoef)
    angles = len(total_angles)
    angle_types = len(potential.anglecoef)
    dihedrals = len(proper_dihedrals)
    dihedral_types = len(potential.dihedralcoef)

    '''
    impropers = len(improper_dihedrals)
    improper_types = len(potential.impropercoef)
    '''

    def __initial_part__():
        ''' this function write the initial part of LAMMPS data file including
            number of atoms, bonds, angles, dihedrals, impropers, size of the 
            simulation box, etc. '''
            
        f.write('LAMMPS data file via write_data, version 7 Aug 2019\n\n')
        # Number of atoms and number of their types in the system
        f.write(str(len(atoms))+' atoms\n')        
        f.write(str(atom_types)+' atom types\n')       
        # Number of bonds in the system and their types 
        bonds_count = 0
        for i in range(len(bonds)):
            if i < len(CC_bonds):
                for j in range(4):
                    if bonds[i,j] >= 0 and bonds[i,j]>i:
                        bonds_count = bonds_count + 1
            elif i >= len(CC_bonds) and i < (len(bonds)):                
                for j in range(4):
                    if bonds[i,j] >=0:
                        bonds_count = bonds_count + 1
        f.write(str(bonds_count)+' bonds\n')        
        f.write(str(bond_types)+'  bond types\n')
        # Number of angles and their types 
        f.write(str(angles)+'  angles\n')
        f.write(str(angle_types)+'  angle types\n')
        # Number of dihedrals and their types 
        f.write(str(dihedrals)+'  dihedrals\n')        
        f.write(str(dihedral_types)+'  dihedral types\n\n')

        '''
        # Number of impropers and their types
        f.write(str(impropers)+'  impropers\n')        
        f.write(str(improper_types)+'  improper types\n\n')       
        '''
    # The size of the simulation box
    def __box__(unit_cell):
        ''' This function find the size of the simulation box based on maximum 
            and minimum size of the x, y, and z. Also the orientation of the 
            box, if it is triclinic '''
        xlo = sign(numpy.amin(atoms[:,0])) + numpy.amin(atoms[:,0])
        xhi = sign(numpy.amax(atoms[:,0])) + numpy.amax(atoms[:,0])
        ylo = sign(numpy.amin(atoms[:,1])) + numpy.amin(atoms[:,1])
        yhi = sign(numpy.amax(atoms[:,1])) + numpy.amax(atoms[:,1])
        zlo = sign(numpy.amin(atoms[:,2])) + numpy.amin(atoms[:,2])
        zhi = sign(numpy.amax(atoms[:,2])) + numpy.amax(atoms[:,2])
        xy = 0.0
        xz = unit_cell[2,2]*numpy.cos(99.5)
        yz = unit_cell[2,2]*numpy.sin(99.5)
        
        f.write(str(xlo)+' '+str(xhi)+'  xlo xhi\n'+str(ylo)+' '+str(yhi)+
                '  ylo yhi\n'+str(zlo)+' '+str(zhi)+'  zlo zhi\n'+
                str(xy)+' '+str(xz)+' '+str(yz)+' xy xz yz\n\n')
    
    def __masses__():
        ''' this function write the mass of each type of atoms that we have '''
        Masses = numpy.array([12.0112, 1.00797])   
        f.write('Masses\n\n1 '+str(Masses[0])+'\n2 '+str(Masses[1])+'\n\n')        
               # +'\n3 '+str(Masses[2])+'\n4 '+str(Masses[3])+'\n\n')


    def __coeffs__():
        ''' this function will write the coefficient of the potential from the 
            respective file containing them '''

        # Bond related coefficient
        f.write('Pair Coeffs # lj/class2/coul/long\n\n')
        potential.pair_coeffs(f)
        f.write('\nBond Coeffs #  class2\n\n')
        potential.bond_coeffs(f)

        # Angle related coefficients 
        f.write('\nAngle Coeffs # class2\n\n')
        potential.angle_coeffs(f)       
        f.write('\nBondBond Coeffs\n\n')
        potential.bond_bond_coeffs(f)
        f.write('\nBondAngle Coeffs\n\n')
        potential.bond_angle_coeffs(f)
        
        # Dihedral related coefficients 
        f.write('\nDihedral Coeffs # class2\n\n')
        potential.dihedral_coeffs(f)
        f.write('\nAngleAngleTorsion Coeffs\n\n')
        potential.angle_angle_torsion_coeffs(f)
        f.write('\nEndBondTorsion Coeffs\n\n')
        potential.end_bond_torsion_coeffs(f)
        f.write('\nMiddleBondTorsion Coeffs\n\n')
        potential.middle_bond_torsion_coeffs(f)
        f.write('\nBondBond13 Coeffs\n\n')
        potential.bond_bond13_coeffs(f)
        f.write('\nAngleTorsion Coeffs\n\n')
        potential.angle_torsion_coeffs(f)
        
        '''
        # Improper related coefficients 
        f.write('\nImproper Coeffs # class2\n\n')
        potential.improper_coeffs(f)
        f.write('\nAngleAngle Coeffs\n\n')
        potential.angle_angle_coeffs(f)
        '''

    def __atoms__():
        ''' Write atomic properties: id molecular_tag charge x y z '''
        
        f.write('\nAtoms # full\n\n')
        def __charge__(i, qkey):
            ''' this function add the charge of each atom to the data file '''
            # if we use the real charges and the system might not be neutral 
            if qkey == 1:
                if i < len(C):
                    if i % 3 == 0:
                        q = -1.590000
                    elif i % 3 == 1:
                        q = -5.299999
                    elif i % 3 == 2:
                        q = -1.060000
                else:
                    q = 5.299999
            # if we do not need the charges and we want a neutral system
            elif qkey == 0:
                if i % 2 == 0:
                    q = 1
                else:
                    q = -1
            return q

        C_counter = 0
        H_counter = 0
        mol_tag = numpy.zeros(len(atoms))

        for i in range(len(atoms)):
            if i < len(C):
                atom_type = 1
                if i >= C_counter*9 and i < (C_counter+1)*9:
                    mol_tag[i] = C_counter+1
                    f.write(str(int(i+1))+' '+str(int(mol_tag[i]))+' '+str(atom_type)+' '+
                            str(__charge__(i))+' '+str(atoms[i,0])+' '+str(atoms[i,1])+' '+
                            str(atoms[i,2])+'\n')
                else:                    
                    C_counter = C_counter+1
                    mol_tag[i] = C_counter+1
                    f.write(str(int(i+1))+' '+str(int(mol_tag[i]))+' '+str(atom_type)+' '+
                            str(__charge__(i))+' '+str(atoms[i,0])+' '+str(atoms[i,1])+' '+
                            str(atoms[i,2])+'\n')
            elif i >= len(C):
                atom_type = 2
                if i >= (H_counter*18+len(C)) and i < ((H_counter+1)*18+len(C)):
                    mol_tag[i] = H_counter+1
                    f.write(str(int(i+1))+' '+str(int(mol_tag[i]))+' '+str(atom_type)+' '+
                            str(__charge__(i))+' '+str(atoms[i,0])+' '+str(atoms[i,1])+
                            ' '+str(atoms[i,2])+'\n')
                else:
                    H_counter = H_counter+1
                    mol_tag[i] = H_counter+1
                    f.write(str(int(i+1))+' '+str(int(mol_tag[i]))+' '+str(atom_type)+' '+
                            str(__charge__(i))+' '+str(atoms[i,0])+' '+str(atoms[i,1])+
                            ' '+str(atoms[i,2])+'\n')
                    
    def __velocities__():
        ''' Writing velocities: id vx vy vz '''
        f.write('\nVelocities\n\n')

    def __bonds__():
        ''' Write bonding details: id bond_type atom-1 atom-2'''
        f.write('\nBonds\n\n')
        counter = 0
        for i in range(len(bonds)):
            # carbon - carbon bonds
            if i < len(CC_bonds):
                bond_type = 1
                for j in range(4):
                    if bonds[i,j] >= 0 and bonds[i,j]>i:
                        counter = counter + 1
                        f.write(str(int(counter))+' '+str(int(bond_type))+' '+str(int(i+1))+' '+
                                str(int(bonds[i,j])+1)+'\n')
            elif i >= len(CC_bonds):
                # carbon - hydrogen bonds
                bond_type = 2
                for j in range(4):
                    if bonds[i,j] >=0:
                        counter = counter + 1
                        f.write(str(int(counter))+' '+str(int(bond_type))+' '+
                                str(int(i+1-len(CC_bonds)))+' '+str(int(bonds[i,j]))+'\n')

    def __angles__(total_angles):
        ''' Writing angle details: id angle-type atom-1 atom-2 atom-3 
            (atom-2 is the center atom)'''
        f.write('\nAngles\n\n')
        for i in range(len(total_angles)):
            angle_type = total_angles[i][3]
            f.write(str(i+1)+' '+str(angle_type)+' '+str(total_angles[i][0]+1)+' '+
                    str(total_angles[i][1]+1)+' '+str(total_angles[i][2]+1)+'\n')

    def __dihedrals__():
        ''' Writing dihedral details: id dihedral_type atom-1 atom-2 atom-3 atom-4 
            (atoms 2-3 form central bond)'''
        f.write('\nDihedrals\n\n')
        for i in range(len(proper_dihedrals)):
            dihedral_type = proper_dihedrals[i][4]
            f.write(str(i+1)+' '+str(dihedral_type)+' '+str(proper_dihedrals[i][0]+1)+' '+
                    str(proper_dihedrals[i][1]+1)+' '+str(proper_dihedrals[i][2]+1)+' '+
                    str(proper_dihedrals[i][3]+1)+'\n')

    def __impropers__():
        ''' Writing impropers details: id improper_type atom-1 atom-2 atom-3 atom-4 
            (atom-2 is the central atom)'''
        f.write('\n\nImpropers\n\n')



    __initial_part__()
    __box__(unit_cell)
    __masses__()
    __coeffs__()
    __atoms__()
    #__velocities__()
    __bonds__()
    __angles__(total_angles)
    __dihedrals__()
    #__impropers__()

