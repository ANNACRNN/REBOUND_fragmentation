#!/usr/bin/env python

import numpy as np
import time
import sys

start_time = time.time()
def track_composition():
    f = open("composition_input.txt", 'r')
    init_compositions = [line.split() for line in f.readlines()]
    for line in init_compositions:
        try:
            sum([float(x) for x in line[2:]]) == 1.0
        except:
            print ('ERROR: Realative abundances do not add up to 1.0')
            sys.exit(1)
    try:
        init_hashes = [int(x[0]) for x in init_compositions]
    except:
        init_hashes = [x[0].value for x in init_compositions]
    no_species = len(init_compositions[0])-2
    file = open("collision_report.txt", 'r')
    blocks = file.read().split("\n")
    blocks = [block for block in blocks if len(block) > 0]
    compositions = init_compositions
    for i in range(len(blocks)):
        block = blocks[i].split()
        time = float(block[0])
        collision_type = int(block[1])
        if collision_type == 0:
            continue
        target_hash = int(block[2])
        target_mass = float(block[3])
        projectile_hash = int(block[4])
        targ_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==target_hash][0]
        proj_idx = [i for i in range(len(compositions)) if int(compositions[i][0])==projectile_hash][0]    
        last_target_mass = float(compositions[targ_idx][1])
        last_proj_mass = float(compositions[proj_idx][1])
        last_target_abundances = compositions[targ_idx][-no_species:]   
        last_projectile_abundances = compositions[proj_idx][-no_species:] 
        no_frags = int((len(block)-5)/2)
        frag_hashes = [int(block[i*2+3]) for i in range(1,no_frags+1)]
        frag_masses = [float(block[i*2+4]) for i in range(1,no_frags+1)]
        if collision_type == 1: #perfect merger
            for i in range(no_species):
                compositions[targ_idx][i+2]=(float(last_target_abundances[i])*last_target_mass+float(last_projectile_abundances[i])*last_proj_mass)/target_mass
        if collision_type == 2: #partial accretion
            mass_accreted = target_mass-last_target_mass
            for i in range(no_species):
                compositions[targ_idx][i+2]=(float(last_target_abundances[i])*last_target_mass+float(last_projectile_abundances[i])*mass_accreted)/target_mass
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j]]+last_projectile_abundances
                compositions.append(frag_data)
                try:
                     any(n < 0 for n in frag_data) == False
                except:
                    print ('ERROR: Negative value encountered in frag data at', time)
                    sys.exit(1)
        if collision_type == 3 or collision_type == 4: #partial errosion, target abundances stay the same
            mass_lost = last_target_mass-target_mass
            frag_abundances = [(float(last_target_abundances[i])*mass_lost+float(last_projectile_abundances[i])*last_proj_mass)/np.sum(frag_masses) for i in range(no_species)]
            for j in range(no_frags):
                frag_data = [frag_hashes[j], frag_masses[j]]+frag_abundances
                try:
                     any(n < 0 for n in frag_data) == False
                except:
                    print ('ERROR: Negative value encountered in frag data at', time)
                    sys.exit(1)
                compositions.append(frag_data)
        compositions[targ_idx][1]=target_mass
        try: 
            any(n < 0 for n in compositions[targ_idx]) == False
        except:
            print ('ERROR: Negative value encountered at', time)
            sys.exit(1)

    
    return compositions

def write_output(compositions):
    f = open("composition_output.txt", "w")
    for body in compositions:
        proper = '[%s]' % ' '.join(map(str, body))
        f.write(proper[1:-1]+"\n")
    f.close()

write_output(track_composition())


print("--- %s seconds ---" % (time.time() - start_time))
