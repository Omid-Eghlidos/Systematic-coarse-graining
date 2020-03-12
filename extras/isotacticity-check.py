#!/usr/bin/env python3

import numpy
import io
import math
from matplotlib import pyplot


def read_lammps_timestep(dump_file, stop_step=1000,start_step=0, specific_step=None, only_first_step=False):
    ''' Given path to dump file, read the next timestep. '''
    if isinstance(dump_file, file):
        f = dump_file
    else:
        f = open(dump_file)
    dump = []
    step = 0
    while f and step <= stop_step:
        line = f.readline()
        if step < start_step:
            if line.startswith('ITEM: TIMESTEP'):
                print 'Skipping step ',step
                step += 1
                if step == start_step:
                    dump.append(lammps_timestep())
                    dump[-1].timestep = int(f.readline())
            continue

        if not line:
            break
        elif line.startswith('ITEM: TIMESTEP'):
            step += 1
            dump.append(lammps_timestep())
            dump[-1].timestep = int(f.readline())
            print 'Reading timestep {}'.format(dump[-1].timestep)
        elif line.startswith('ITEM: BOX BOUNDS'):
            for i in range(3):
                dump[-1].box[i] = [float(s) for s in f.readline().split()]
        elif line.startswith('ITEM: NUMBER OF ATOMS'):
            dump[-1].set_num_atoms(int(f.readline()))
        # ITEM: ATOMS id mol type x y z
        elif line.startswith('ITEM: ATOMS'):
            if specific_step:
                if dump[-1].timestep != specific_step:
                    for _ in dump[-1].type:
                        row = f.readline().split()
                    continue
            var = {v:i for i,v in enumerate(line.split()[2:])}
            for _ in dump[-1].type:
                row = f.readline().split()
                i = int(row[var['id']]) - 1
                dump[-1].coord[i,0] = float(row[var['x']])
                dump[-1].coord[i,1] = float(row[var['y']])
                dump[-1].coord[i,2] = float(row[var['z']])
                dump[-1].molecule[i] = int(row[var['mol']])
                dump[-1].type[i] = int(row[var['type']])
            if specific_step:
                break
            if only_first_step:
                return dump
    del dump[-1]
    return dump


dump = read_lammps_timestep('sample-01.lammpstrj', stop_step=1000,start_step=0, specific_step=None, only_first_step=False)

print(dump[0,:])



def isotacticity():

    ''' This function check to see after any simulation if the system remains in its initial tacticity form'''
    return data







