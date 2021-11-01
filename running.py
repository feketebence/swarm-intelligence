import math

from initializer import g
import initializer

def distance_squared_folded_PBC(x0: float, y0: float, x1: float, y1: float):

    dx = x1 - x0
    dy = y1 - y0

    # PBC fold back
    # if any distance is larger than half the box
    # the copy in the neighboring box is closer
    if dx > g['halfSX']: dx -= g['SX']
    if dx <= -g['halfSX']: dx += g['SX']
    if dy > g['halfSY']: dy -= g['SY']
    if dy <= -g['halfSY']: dy += g['SY']

    dr2 = dx*dx + dy*dy

    return [dr2, dx, dy]

def fold_particle_back_PBC(i: int):
    # fold back the particle into the PBC simulation box
    # assumes it did not jump more than a box length
    # if it did the simulation is already broken anyhow

    if (g['particle_x'][i] < 0): g['particle_x'][i] += g['SX']
    if (g['particle_y'][i] < 0 ): g['particle_y'][i] += g['SY']
    if (g['particle_x'][i] >= g['SX']): g['particle_x'][i] -= g['SX']
    if (g['particle_y'][i] >= g['SY']): g['particle_y'][i] -= g['SY']

def move_particles():
    for i in range(g['N_particles']):
        dx = g['particle_fx'] * g['dt']
        dy = g['particle_fy'] * g['dt']

        g['particle_x'][i] += dx
        g['particle_y'][i] += dy

        g['particle_dx_so_far'][i] += dx
        g['particle_dy_so_far'][i] += dy

        fold_particle_back_PBC(i)

        g['particle_fx'][i] = 0.0
        g['particle_fy'][i] = 0.0

def calculate_external_forces_on_particles():
    for i in range(g['N_particles']):
        g['particle_fx'][i] += g['particle_direction'][i] * g['particle_driving_force']

def calculate_pairwise_forces():
    for idx in range(g['N_verlet']):
        i = g['Verletlisti'][idx]
        j = g['Verletlistj'][idx]

        # calculate the distance between the particles
        r2, dx, dy = distance_squared_folded_PBC(
            g['particle_x'][i], 
            g['particle_y'][i], 
            g['particle_x'][j], 
            g['particle_y'][j]
        )

        # the particles are too close to interact
        if r2 < 16.0:
            r = math.sqrt(r2)

            if r < 0.2:
                f = 100.0
                print(f"[ WARNING ] Particles too close. Lower cut-off force used, f = {f}")
            else:
                f = 1/r2 * math.exp(-r / g['particle_particle_screening_length'])

            # projection to the x and y axis
            f = f/r

            g['particle_fx'][i] -= f * dx
            g['particle_fy'][i] -= f * dy

            g['particle_fx'][j] += f * dx
            g['particle_fy'][j] += f * dy




def run_simulation():
    g['total_time'] = 1_000
    g['echo_time'] = 100        # ennyi lepesenkent ir konzolra
    g['movie_time'] = 100       # ennyi lepesenkent ir .mvi fileba
    g['statistics_time'] = 100  # ennyi lepesenkent ir statisztikas .txt fileba

    g['temperature'] = 0.0

def print_globals_keys():
    print("===========-- global g's keys --===========")
    for elem in g.keys():
        print('\t', elem)
    print("===========-- -- -- -- -- -- --===========")

def main():
    initializer.init_simulation()
    initializer.init_simulation_box()
    initializer.init_particles()
    initializer.init_files()
    print_globals_keys()

if __name__ == "__main__":
    main()
