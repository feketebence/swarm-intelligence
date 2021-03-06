import math
import random
import time
from pprint import pprint

import numpy as np

g = {}

def init_simulation():

    g['start_time_epoch'] = time.time()

    g['dt'] = 0.001
    g['particle_particle_screening_length'] = 4.0
    g['particle_driving_force'] = 0.5
    print(f"Timestep (dt) = {g['dt']}")
    print(f"Screening length = {g['particle_particle_screening_length']}")
    print(f"Driving force on particles = {g['particle_driving_force']}")

    g['Verlet_cutoff_distance'] = 1.5 * g['particle_particle_screening_length']
    g['Verlet_cutoff_distance_squared'] = g['Verlet_cutoff_distance'] ** 2
    g['Verlet_intershell_squared'] = g['Verlet_cutoff_distance'] - g['particle_particle_screening_length']
    g['Verlet_intershell_squared'] = g['Verlet_intershell_squared'] / 2.0
    g['Verlet_intershell_squared'] *= g['Verlet_intershell_squared']

    print(f"Verlet cutoff distance = {g['Verlet_cutoff_distance']}")
    print(f"Half of Verlet intershell distance = {math.sqrt(g['Verlet_intershell_squared'])}")
    print(f"Half of Verlet intershell distance squared = {g['Verlet_intershell_squared']}")

    # zero everything so rebuild Verlet can find this the first time
    g['Verletlisti'] = None
    g['Verletlistj'] = None
    g['N_verlet'] = 0
    g['N_verlet_max'] = 0

    print()

def init_simulation_box():
    print("Initializing the simulation box")
    
    g['SX'] = 12.0
    g['SY'] = 12.0
    g['halfSX'] = g['SX'] / 2.0
    g['halfSY'] = g['SY'] / 2.0

    print(f"SX = {g['SX']}, SY = {g['SY']}")
    print(f"halfSX = {g['halfSX']}, halfSY = {g['halfSY']}")

    print()

def distance_folded_PBC(x0: float, y0: float, x1: float, y1: float):
    """
    Calculates the shortest distance between 2 points in a PBC configuration
    only used in the random deposition check which happens once at initialization
    so this is not so time crucial left the square root inside
    """

    dx: float = x1 - x0
    dy: float = y1 - y0

    # PBC foldback: if the distance is larger than half the box the copy in the neighboring box is closer
    if dx > g['halfSX']: dx -= g['SX']
    if dx <= -g['halfSX']: dx += g['SX']
    
    if dy > g['halfSY']: dx -= g['SY']
    if dy <= -g['halfSY']: dx += g['SY']

    r = math.sqrt(dx * dx + dy * dy)

    return r

def init_particles_randomly():
    r_min = 0.25

    for i in range(g['N_particles']):
        x_try = 0.0
        y_try = 0.0

        # check overlap with previous particles, assume there is overlap to get into the cycle
        overlap = True
        N_trials = 0

        while overlap == True and N_trials < g['N_particles']:
            # attempt to place the particle
            x_try = g['SX'] * random.uniform(0, 1)
            y_try = g['SY'] * random.uniform(0, 1)

            # assume this was good
            overlap = False

            for j in range(i):
                # calculate distance
                dr = distance_folded_PBC(x_try, y_try, g['particle_x'][j], g['particle_y'][j])

                if dr < r_min:
                    overlap = True
                    N_trials += 1
                    break
        
        if N_trials == g['N_particles'] ** 2:
            print('Cannot place particles randomly, system too dense')
            print('Quitting...')
            exit(1)

        g['particle_x'][i] = x_try
        g['particle_y'][i] = y_try
        g['particle_fx'][i] = 0.0
        g['particle_fy'][i] = 0.0
        g['particle_dx_so_far'][i] = 0.0
        g['particle_dy_so_far'][i] = 0.0
        
        if (random.uniform(0, 1) < 0.5):
            g['particle_direction'][i] = -1.0
            g['particle_color'][i] = 2
        else:
            g['particle_direction'][i] = +1.0
            g['particle_color'][i] = 3

    print("Random arrangement of particles initialized")
    print(f"N_particles = {g['N_particles']} placed")
    print()

def init_particles():
    g['N_particles'] = 30

    n_particles = g['N_particles']
    g['particle_x'] = np.ones(n_particles)
    g['particle_y'] = np.ones(n_particles)
    g['particle_fx'] = np.ones(n_particles)
    g['particle_fy'] = np.ones(n_particles)
    g['particle_color'] = np.zeros(n_particles)
    g['particle_dx_so_far'] = np.ones(n_particles)
    g['particle_dy_so_far'] = np.ones(n_particles)

    g['particle_direction'] = np.ones(n_particles)

    init_particles_randomly()

    print("Particles initialized\n")
    print(f"N_particles = {g['N_particles']}")

def init_files():

    g['moviefile_name'] = "particles.mvi"
    g['moviefile'] = open(g['moviefile_name'], "wb") # use wt for text mode, b is binary mode
    if g['moviefile'] == None:
        print(f"Could not create/open movie file: {g['moviefile_name']}")
        exit(2)

    g['statisticsfile_name'] = "statistics.txt"
    g['statisticsfile'] = open(g['statisticsfile_name'], "wt")
    if g['statisticsfile'] == None:
        print(f"Could not create/open movie file: {g['statisticsfile_name']}")
        exit(2)

    g['stat_all_movement_red_x'] = 0.0
    g['stat_all_movement_blue_x'] = 0.0
    g['stat_all_movement'] = 0.0

def print_globals():
    print("Globals:")
    pprint(g)

def main():
    init_simulation()
    init_simulation_box()
    init_particles()
    init_files()
    print_globals()

if __name__ == "__main__":
    main()
