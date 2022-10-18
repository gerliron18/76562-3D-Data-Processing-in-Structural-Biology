import numpy as np
from brownian_dynamics import Particle
from brownian_dynamics import BrownianDynamics
from brownian_dynamics.BDConfiguration import BDConfiguration

TCR_DENSITY = 0.02 #300  # number of TCRs in 1 micron^2
CD45_DENSITY = 0.06 #1000  # number of CD45s in 1 micron^2
PIXEL_TO_MICRON = 0.01  # number of micron^2 in 1 pixel


def convert_pixels_to_micron(max_x, max_y):
    """

    :param max_x: maximum x coordinate value
    :param max_y: maximum x coordinate value
    :return: the number of micron^2 in the given grid size
    """
    num_pixels = max_x * max_y
    return num_pixels * PIXEL_TO_MICRON


def random_BD(max_x, max_y):
    """
    Initializes a Brownian Dynamics configuration, where the particles are
    uniformly distributed in the cell
    :param max_x: maximum size of grid for x axis
    :param max_y: maximum size of grid for y axis
    :return: BDConfiguration objects with the randomly generated particles
    """

    num_micron = convert_pixels_to_micron(max_x, max_y)
    num_TCR = int(TCR_DENSITY * num_micron)
    num_CD45 = int(CD45_DENSITY * num_micron)

    loc_range = np.array([[0, max_x], [0, max_y]])

    tcr_loc = np.random.uniform(loc_range[:, 0], loc_range[:, 1],
                                size=(num_TCR, loc_range.shape[0]))

    cd45_loc = np.random.uniform(loc_range[:, 0], loc_range[:, 1],
                                 size=(num_CD45, loc_range.shape[0]))

    tcr_particles = [Particle.TCR(tcr_loc[i][0], tcr_loc[i][1]) for i in
                     range(num_TCR)]
    cd45_particles = [Particle.CD45(cd45_loc[i][0], cd45_loc[i][1]) for i in
                      range(num_CD45)]

    # If adding additional particle types, need to implement the following
    # lines:
    # num_NEW_PARTICLE = int(NEW_PARTICAL_DENSITY * num_micron)
    # NEW_PARTICLE_loc = np.random.uniform(loc_range[:, 0], loc_range[:, 1],
    #                              size=(num_NEW_PARTICLE, loc_range.shape[0]))
    # NEW_PARTICLE_particles = [Particle.CD45(cd45_loc[i][0], cd45_loc[i][1])
    # for i in range(num_NEW_PARTICLE)]
    # # add the NEW_PARTICLE_particles to the all_particles list


    all_particles = tcr_particles + cd45_particles


    return BDConfiguration(all_particles, (max_x, max_y))



if __name__ == '__main__':
    print(random_BD(2, 6).get_particles())
