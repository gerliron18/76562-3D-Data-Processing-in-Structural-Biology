import copy

import numpy as np
import tqdm


class BrownianDynamics(object):
    def __init__(self, membrane_borders,
                 parameters_vec=None, force_function=None, kbt=1e-11):
        """
        A class to represent a Brownian Dynamics algorithm
        :param parameters_vec: vector of param to be used in different parts useful for MC
        :param membrane_borders: The membrane borders, can be whatever the user decided
        :param force_function: the force function
        """
        self.kbt = kbt
        self.max_x, self.max_y = membrane_borders
        self.parameters_vec = parameters_vec
        self.force_function = force_function
        self.snapshots = []

    def get_final_snapshot(self):
        """

        :return: a list of all the particles' location at the end of the
        simulation
        """
        return self.snapshots[-1]

    def get_snapshots(self):
        """

        :return: a list of the particles' location vectors for each
        iteration of the simulation
        """
        return self.snapshots

    def run_brownian_dynamics(self, particle_list, num_iter, dt):
        """
        :param particle_list: list with all the particles
        :param num_iter: number of iterations of the simulation
        :param dt: the deltaT of each step
        :return: None
        """
        if num_iter <= 1:
            raise ValueError("Number of iterations needs to be bigger than 1")

        self.snapshots = []
        location_vec = self.get_location_vector(particle_list)  # vec with

        # two columns - first for x coordinate and second for y coordinate
        type_vector = self.get_type_vector(particle_list)
        diffusion_vector = self.get_diffusion_vector(particle_list)
        equation_first_part = ((dt / self.kbt) * diffusion_vector).reshape(-1,1)
        equation_second_part = (np.sqrt(6 * diffusion_vector * dt)).reshape(-1,1)

        self.snapshots.append((self.get_location_vector(particle_list), type_vector))

        for _ in tqdm.trange(int(num_iter)):
            f = self.force_function(particle_coord_array=location_vec, particle_type_array=type_vector,
                                    max_grid_x=self.max_x, max_grid_y=self.max_y,
                                    mc_cooficients=self.parameters_vec)

            location_vec = location_vec + equation_first_part * f + equation_second_part * np.random.normal(
                size=location_vec.shape)

            location_vec = self.handle_membrane(location_vec)

            self.snapshots.append((copy.deepcopy(location_vec), type_vector))

    def handle_membrane(self, location_vector):
        """
        Function that handles what happens if particles are moved out of
        the frame of the membrane. In this case the particle gets "stuck" in
        the membrane (i.e., their x or y location becomes the maximum value
        of the x or y coordinate
        :param location_vector: ndarray where rows indicated particles and
        column 0 is the x coordinate and column 1 is the y coordinate of the particle
        :return: updated location vector with legal locations
        """
        ind_x = np.where(location_vector[:, 0] > self.max_x)
        location_vector[ind_x, 0] = self.max_x
        ind_x = np.where(location_vector[:, 0] < 0)
        location_vector[ind_x, 0] = 0

        ind_y = np.where(location_vector[:, 1] > self.max_y)
        location_vector[ind_y, 1] = self.max_y
        ind_y = np.where(location_vector[:, 1] < 0)
        location_vector[ind_y, 1] = 0

        return location_vector

    @staticmethod
    def get_location_vector(particle_list):
        """
        Convert a list of particles to a numpy array of locations (x,y)
        :param particle_list: A list of particles
        :return: A numpy array of (x,y)
        """
        return np.array(
            [particle.get_location() for particle in particle_list])

    @staticmethod
    def get_type_vector(particle_list):
        """
        Convert a list of particles to a numpy array of the Particle types

        :param particle_list: list of particles
        :return: numpy array with each of the particles type
        """
        return np.array(
            [particle.get_particle_type() for particle in particle_list])

    @staticmethod
    def get_diffusion_vector(particle_list):
        """
        Convert a list of particles to a numpy array of diffusion coefficients

        :param particle_list: list of Particles objects
        :return: numpy array with each of the particles diffusion coefficient
        """
        return np.array(
            [particle.get_diffusion_coef() for particle in particle_list])
