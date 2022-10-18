import numpy as np
import pandas as pd
from scipy.spatial import distance


def location_interaction(particle_coord_array, particle_type_array, max_grid_x, max_grid_y, mc_cooficients):
    """
    Calculates the total forces acting on each particle in a given list of particles.
    :param particle_coord_array: numeric 2 * num_of_coords numpy array, the first column should contain the x coordinates
    and the second the y coordinates
    :param particle_type_array: string 1 * num_of_coords numpy array, holding the type of each particle
    :param max_grid_x: maximum x coordinate of the grid
    :param max_grid_y: maximum y coordinate of the grid
    :param mc_cooficients: the mcmc coefficients needed for this simulation, for this function one coefficient for the
    interaction force and two more for the location function.
    :return: numeric 2 * num_of_coords numpy array, the first column contains the forces along the x axis for each
    particle, and the second the forces along the y axis
    """
    particle_types = np.unique(particle_type_array)
    particle_inds = {}
    for pt in particle_types:
        particle_inds[pt] = np.where(particle_type_array == pt)
    interaction_mat = _interaction_force(particle_coord_array, particle_type_array, particle_inds, max_grid_x)

    radius_list = _location_force(particle_coord_array, max_grid_x, max_grid_y, particle_inds,
                                  mc_cooficients["location_b"])
    total_force_1d = mc_cooficients["interaction"] * interaction_mat + mc_cooficients["location_a"] * radius_list
    angles = _get_angles(particle_coord_array, max_grid_x, max_grid_y)
    split_force = np.empty((total_force_1d.shape[0], 2))
    split_force[:, 0] = total_force_1d * np.cos(angles)
    split_force[:, 1] = total_force_1d * np.sin(angles)
    return split_force


def _interaction_force(particle_coord_array, particle_type_array, particle_inds, max_grid_x):
    """
    For every particle in the particle list calculates all the forces each particle applies on it.
    :param particle_coord_array: Array of the particles coordinates
    :param particle_type_array: Array of the particles types
    :param particle_inds: A dictionary, keys are the particle types and the value a list of their indices in the
    coordinates and types array
    :param max_grid_x: The maximum x coordinate
    :return: 1D array containing the interaction forces
    """
    value_mat = _get_interaction_mat(particle_inds.keys(), max_grid_x)
    type_arr = np.empty((particle_type_array.shape[0], particle_type_array.shape[0]))
    for pt_x in particle_inds:
        for pt_y in particle_inds:
            val = value_mat.loc[pt_x, pt_y]
            type_arr[np.ix_(particle_inds[pt_x][0], particle_inds[pt_y][0])] = val

    np.fill_diagonal(type_arr, 0)
    dist_arr = distance.cdist(particle_coord_array[:, 0:2], particle_coord_array[:, 0:2], 'euclidean')
    np.fill_diagonal(dist_arr, 1)  # so we dont divide by zero
    dist_arr[dist_arr < 0.001] = 0.001
    interaction_mat = type_arr / dist_arr
    return np.sum(interaction_mat, axis=0)


def _get_interaction_mat(particle_types, max_grid_x):
    """
    Calculates the force factor for every pair of particles, if a particle pair isn't implemented, gives a zero
    :param particle_types: List of the unique particle types
    :param max_grid_x: The maximum x coordinate
    :return: A pandas dataframe with the forces for every pair of particles
    """
    saved_vals = {}
    value_mat = pd.DataFrame(columns=particle_types, index=particle_types)
    for pt_x in particle_types:
        for pt_y in particle_types:
            if pt_x == pt_y == "TCR":
                value_mat.loc[pt_x, pt_y] = np.random.normal(max_grid_x / 10, 1)
            elif pt_x == pt_y == "CD45":
                value_mat.loc[pt_x, pt_y] = np.random.normal(max_grid_x / 30, 1)
            elif (pt_x == "CD45" and pt_y == "TCR") or (pt_x == "TCR" and pt_y == "CD45"):
                if "cd_tcr" not in saved_vals:
                    saved_vals["cd_tcr"] = np.random.normal(- max_grid_x / 50, 1)
                value_mat.loc[pt_x, pt_y] = saved_vals["cd_tcr"]
            else:
                value_mat.loc[pt_x, pt_y] = 0
    return value_mat


def _location_force(particle_coord_array, max_grid_x, max_grid_y, particle_inds, b):
    """
    Calculates the force on each particle in the array depending on it's location, if the particle isn't implemented
    gives zero.
    :param particle_coord_array: Array of the particles coordinates
    :param max_grid_x: The maximum x coordinate
    :param max_grid_y: The maximum y coordinate
    :param particle_inds: A dictionary, keys are the particle types and the value a list of their indices in the
    coordinates array
    :param b: mcmc factor for the CD45s max location force
    :return: 1D array containing the location forces
    """
    radius_list = np.sqrt(np.power(particle_coord_array[:, 0] - (max_grid_x / 2), 2) +
                          np.power(particle_coord_array[:, 1] - (max_grid_y / 2), 2))
    for pt in particle_inds:
        if pt == "CD45":
            max_energy_cd = max(max_grid_x / b, max_grid_y / b)
            radius_list[particle_inds[pt]] = (radius_list[particle_inds[pt]] * 2) - max_energy_cd
        elif pt == "TCR":
            radius_list[particle_inds[pt]] = radius_list[particle_inds[pt]] * 2
        else:
            radius_list[particle_inds[pt]] = 0
    return radius_list


def _get_angles(particle_coord_array, max_grid_x, max_grid_y):
    """
    Finds the angle of each point in the array to the cells centre
    :param particle_coord_array: Array of the particles coordinates
    :param max_grid_x: The maximum x coordinate
    :param max_grid_y: The maximum y coordinate
    :return: Each coordinates angle towards the cells centre
    """
    y_coor = particle_coord_array[:, 1] - (max_grid_y / 2)
    x_coor = particle_coord_array[:, 0] - (max_grid_x / 2)
    x_coor[x_coor == 0] = 1
    return np.arctan(y_coor / x_coor)
