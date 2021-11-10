import math
import random
import datetime
import time

import numpy as np

from initializer import g
import initializer


def distance_squared_folded_PBC(x0: float, y0: float, x1: float, y1: float):

    dx = x1 - x0
    dy = y1 - y0

    # PBC fold back
    # if any distance is larger than half the box
    # the copy in the neighboring box is closer
    if dx > g["halfSX"]:
        dx -= g["SX"]
    if dx <= -g["halfSX"]:
        dx += g["SX"]
    if dy > g["halfSY"]:
        dy -= g["SY"]
    if dy <= -g["halfSY"]:
        dy += g["SY"]

    dr2 = (dx * dx) + (dy * dy)

    return dr2, dx, dy


def fold_particle_back_PBC(i: int):
    # fold back the particle into the PBC simulation box
    # assumes it did not jump more than a box length
    # if it did the simulation is already broken anyhow

    if g["particle_x"][i] < 0:
        g["particle_x"][i] += g["SX"]
    if g["particle_y"][i] < 0:
        g["particle_y"][i] += g["SY"]
    if g["particle_x"][i] >= g["SX"]:
        g["particle_x"][i] -= g["SX"]
    if g["particle_y"][i] >= g["SY"]:
        g["particle_y"][i] -= g["SY"]


def move_particles():
    for i in range(g["N_particles"]):
        dx = g["particle_fx"][i] * g["dt"]
        dy = g["particle_fy"][i] * g["dt"]

        g["particle_x"][i] += dx
        g["particle_y"][i] += dy

        g["particle_dx_so_far"][i] += dx
        g["particle_dy_so_far"][i] += dy

        fold_particle_back_PBC(i)

        g["particle_fx"][i] = 0.0
        g["particle_fy"][i] = 0.0


def calculate_external_forces_on_particles():
    for i in range(g["N_particles"]):
        g["particle_fx"][i] += g["particle_direction"][i] * g["particle_driving_force"]


def calculate_pairwise_forces():
    for idx in range(g["N_verlet"]):
        i = g["Verletlisti"][idx]
        j = g["Verletlistj"][idx]

        # calculate the distance between the particles
        r2, dx, dy = distance_squared_folded_PBC(
            g["particle_x"][i],
            g["particle_y"][i],
            g["particle_x"][j],
            g["particle_y"][j],
        )

        # the particles are too close to interact
        if r2 < 16.0:
            r = math.sqrt(r2)
            # print("[ DEBUG , line 85] r =", r)
            # print("[ DEBUG , line 85] r2 =", r2)

            if r < 0.2:
                f = 100.0
                print(f"[ WARNING ] Particles too close. Lower cut-off force used, f = {f}")
            else:
                f = 1 / r2 * math.exp(-r / g["particle_particle_screening_length"])

            # projection to the x and y axis
            if r == 0:
                print("[ WARNING ] line 96, r = 0, using bigger r = 0.00001")
                r = 0.00001
            f /= r

            g["particle_fx"][i] -= f * dx
            g["particle_fy"][i] -= f * dy

            g["particle_fx"][j] += f * dx
            g["particle_fy"][j] += f * dy


def write_cmovie_frame():
    # open the file
    f_movie = open(g["moviefile_name"], "a")

    # write into the file
    f_movie.write(f"{str(g['N_particles'])}\n")
    f_movie.write(f"{str(g['time'])}\n")

    for i in range(g["N_particles"]):
        s = f"{g['particle_color'][i]}, {i}, {g['particle_x'][i]}, {g['particle_y'][i]}, {1.0}\n"
        f_movie.write(s)

    f_movie.write("\n")

def write_cmovie_frame_binary():

    out_file = g['moviefile_name']
    with open(out_file, 'wb') as out_f:
        out_f.write(bytes(g['N_particles']))
        out_f.write(bytes(g['time']))
        
        for i in range(g["N_particles"]):
            out_f.write(bytes(g['particle_color'][i]))
            out_f.write(bytes(i))
            out_f.write(bytes(g['particle_x'][i]))
            out_f.write(bytes(g['particle_y'][i]))
            out_f.write(bytes(1.0))






def check_Verlet_rebuild_condition_and_set_flag():
    g["flag_to_rebuild_verlet"] = False

    for i in range(g["N_particles"]):
        dr2 = g["particle_dx_so_far"][i] ** 2 + g["particle_dy_so_far"][i] ** 2

        if dr2 >= g["Verlet_intershell_squared"]:
            g["flag_to_rebuild_verlet"] = True
            break


def rebuild_Verlet_list():

    # building the Verlet list for the first time in the simulation
    if g["N_verlet_max"] == 0:
        print("Verlet list will be built for the first time.")

        estimation = g["N_particles"] / g["SX"] / g["SY"]
        print(f"System density is {estimation}")

        estimation *= math.pi * g["Verlet_cutoff_distance"] * g["Verlet_cutoff_distance"]
        print(f"Particles in a R = {g['Verlet_cutoff_distance']}, Shell = {estimation}")

        g["N_Verlet_max"] = int(int(estimation) * g["N_particles"] / 2)
        print(f"Estimated N_Verlet_max = {g['N_Verlet_max']}")

        g["Verletlisti"] = np.zeros(g["N_Verlet_max"], dtype=np.int32)
        g["Verletlistj"] = np.zeros(g["N_Verlet_max"], dtype=np.int32)

    # Build the verlet list
    g["N_verlet"] = 0
    n_particles = g["N_particles"]
    
    i = 0
    while i < n_particles:
        
        j = i + 1
        while j < n_particles:
            dr2, dx, dy = distance_squared_folded_PBC(
                g["particle_x"][i],
                g["particle_y"][i],
                g["particle_x"][j],
                g["particle_y"][j],
            )
            # print("dr2:", dr2)

            if (dr2 < 36.0):
                g["Verletlisti"][g["N_verlet"]] = i
                g["Verletlistj"][g["N_verlet"]] = j

                g["N_verlet"] += 1

                if g["N_verlet"] >= g["N_verlet_max"]:
                    print(f"Verlet list reallocated from {g['N_verlet_max']}")
                    g["N_verlet_max"] = int(1.1 * g["N_verlet"])

                    print(g["Verletlisti"].shape)
                    # g["Verletlisti"].resize(g["N_Verlet_max"])
                    # g["Verletlistj"].resize(g["N_Verlet_max"])
                    # g["Verletlisti"] = np.zeros(g["N_Verlet_max"], dtype=np.int32)
                    # g["Verletlistj"] = np.zeros(g["N_Verlet_max"], dtype=np.int32)
                    g["Verletlisti"] = np.pad(g["N_Verlet_max"], [0])
                    g["Verletlistj"] = np.pad(g["N_Verlet_max"], [0])


                    print(f"New Verlet list max size = {g['N_Verlet_max']}")
        j += 1

    i += 1


    print(f"Counted Verlet list length = {g['N_verlet']}")

    g["flag_to_rebuild_verlet"] = False

    for i in range(g["N_particles"]):
        g["particle_dx_so_far"][i] = 0.0
        g["particle_dy_so_far"][i] = 0.0


def calculate_statistics():
    all_movement_x_blue = 0.0
    all_movement_x_red = 0.0
    all_movement = 0.0

    N_red = 0
    N_blue = 0

    for i in range(g["N_particles"]):
        if g["particle_color"][i] == 2:  # RED goes left
            all_movement_x_red += g["particle_fx"][i]
            all_movement -= g["particle_fx"][i]
            N_red += 1
        elif g["particle_color"][i] == 3:  # BLUE goes right
            all_movement_x_blue += g["particle_fx"][i]
            all_movement += g["particle_fx"][i]
            N_blue += 1

    all_movement_x_red = all_movement_x_red / N_red
    all_movement_x_blue = all_movement_x_blue / N_blue
    all_movement = all_movement / (N_red + N_blue)

    g["stat_all_movement_red_x"] += all_movement_x_red
    g["stat_all_movement_blue_x"] += all_movement_x_blue
    g["stat_all_movement"] += all_movement


def write_statistics():
    g["stat_all_movement_red_x"] = g["stat_all_movement_red_x"] / g["statistics_time"]
    g["stat_all_movement_blue_x"] = g["stat_all_movement_blue_x"] / g["statistics_time"]
    g["stat_all_movement"] = g["stat_all_movement"] / g["statistics_time"]

    # write statistics to file
    f_stat = open(g["statisticsfile_name"], "a")

    # f_stat.write(f"{g['time']} {g['stat_all_movement_red_x']} {g['stat_all_movement_blue_x']} {g['stat_all_movement']}")
    f_stat.write(f"{g['time']} {g['stat_all_movement']}\n")

    g["stat_all_movement_red_x"] = 0.0
    g["stat_all_movement_blue_x"] = 0.0
    g["stat_all_movement"] = 0.0


def calculate_thermal_forces():
    for i in range(g["N_particles"]):
        g["particle_fx"][i] += g["temperature"] * (2.0 * random.uniform(0, 1) - 1.0)
        g["particle_fx"][i] += g["temperature"] * (2.0 * random.uniform(0, 1) - 1.0)


def run_simulation():
    g["total_time"] = 1_000
    g["echo_time"] = 100  # ennyi lepesenkent ir konzolra
    g["movie_time"] = 100  # ennyi lepesenkent ir .mvi fileba
    g["statistics_time"] = 100  # ennyi lepesenkent ir statisztikas .txt fileba

    g["temperature"] = 0.0

    g["time"] = 0

    rebuild_Verlet_list()  # build verlet list for the first time

    start_timestamp = time.time()
    start_datetime = datetime.datetime.fromtimestamp(start_timestamp)
    print("Run started on", start_datetime)

    while g["time"] < g["total_time"]:

        calculate_external_forces_on_particles()
        calculate_pairwise_forces()
        calculate_thermal_forces()

        calculate_statistics()

        move_particles()

        check_Verlet_rebuild_condition_and_set_flag()

        # if one particle moved enough to rebuild the verlet list
        if g["flag_to_rebuild_verlet"] == True:
            rebuild_Verlet_list()

        # print time
        if g["time"] % g["echo_time"] == 0:
            print(f"Timestep: {g['time']} / {g['total_time']}")

        # write to movie file
        if g["time"] % g["movie_time"] == 0:
            write_cmovie_frame()

        # write to statistics file
        if g["time"] % g["statistics_time"] == 0:
            write_statistics()

        g["time"] += 1

    end_timestamp = time.time()
    end_datetime = datetime.datetime.fromtimestamp(end_timestamp)
    print("Run ended on", end_datetime)

    run_time = end_timestamp - start_timestamp
    print("Total runtime in seconds: %s" % (run_time))


def print_globals_keys():
    print("===========-- global g's keys --===========")
    for elem in g.keys():
        print("\t", elem)
    print("===========-- -- -- -- -- -- --===========")


def main():
    initializer.init_simulation()
    initializer.init_simulation_box()
    initializer.init_particles()
    initializer.init_files()
    print_globals_keys()
    run_simulation()


if __name__ == "__main__":
    main()
