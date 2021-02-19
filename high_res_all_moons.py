from all_funcs import save_moon_trace
from multiprocessing import Pool
from timeit import default_timer as timer

# print("Start Unparallel")
# start = timer()
# all_moons = ['Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon', 'Triton']
# num_fieldlines = 0
# num_orbit_points = 2

# for moon in all_moons:
#     if moon == 'Triton':
#         num_orbits = 2
#     else:
#         num_orbits = 1
#     save_moon_trace(moon, num_orbit_points, num_orbits, num_fieldlines)

# print("End Unparallel")
# end = timer()
# print(end-start)

def run(moon):
    num_fieldlines = 0
    num_orbit_points = 500
    if moon == 'Triton':
        num_orbits = 2
    else:
        num_orbits = 1
    save_moon_trace(moon, num_orbit_points, num_orbits, num_fieldlines)

def main():
    all_moons = ['Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon', 'Triton']
    print("Start Parallel")
    start = timer()
    with Pool() as pool:
        pool.map(run, all_moons)
    end = timer()
    print("End Parallel")
    print(end-start)

# if __name__ == '__main__':
#     main()