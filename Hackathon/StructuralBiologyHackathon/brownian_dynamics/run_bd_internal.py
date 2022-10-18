import time

import numpy as np
from matplotlib import pyplot as plt

from brownian_dynamics.Particle import TCR
from brownian_dynamics.BrownianDynamics import BrownianDynamics


def b_d():
    np.random.seed(5)

    T = 1
    N = 501
    dt = T / (N - 1)
    t = np.linspace(dt, T, N)

    dX = np.sqrt(dt) * np.random.randn(1, N)
    # X = np.cumsum(dX)
    X = np.cumsum(dX, axis=1)

    dY = np.sqrt(dt) * np.random.randn(1, N)
    Y = np.cumsum(dY, axis=1)

    fig, ax = plt.subplots()
    ax.plot(X[0, :], Y[0, :])
    ax.plot(X[0, 0], Y[0, 0], 'ro')
    ax.plot(X[0, -1], Y[0, -1], 'yo')
    ax.set_xlabel('X(t)')
    ax.set_ylabel('Y(t)')
    ax.set_title('2D Discretized Brownian Path')
    plt.tight_layout()
    plt.show()


def main():
    b_d()
    t0 = time.time()

    force_matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    membrane_borders = [-100, 100, -100, 100]
    particle_list = []

    # generate random particles
    for i in range(5):
        x, y = np.random.randint(0, 100, 2)
        particle_list.append(TCR(x, y))

    def force_func(x, y): return x + y

    parameters = [0, 1, 2, 3]

    bd = BrownianDynamics(particle_list, force_matrix, membrane_borders,
                          parameters, force_func)

    bd.run_brownian_dynamics(10)

    print((time.time() - t0) / 60)


max_x = 2
max_y = 5


def fix_loc_vec(location_vector):
    ind_x = np.where(location_vector[:, 0] > max_x)
    location_vector[ind_x, 0] = max_x
    ind_x = np.where(location_vector[:, 0] < 0)
    location_vector[ind_x, 0] = 0

    ind_y = np.where(location_vector[:, 1] > max_y)
    location_vector[ind_y, 1] = max_y
    ind_y = np.where(location_vector[:, 1] < 0)
    location_vector[ind_y, 1] = 0


    return location_vector


def test():
    a = np.array([[1, 2, 3, -4, 5], [-2, 3, 6, 7, 4]]).T
    print(a)
    fix_loc_vec(a)
    print(a)


if __name__ == '__main__':
    test()

# if __name__ == '__main__':
#     import Particle
#
#     force_matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
#     membrane_borders = [-100, 100, -100, 100]
#     particle_list = []
#
#     mhc = Particle.TCR(3, 4)
#     print(mhc.diffusion_coef)
#     # generate random particles
#     for i in range(5):
#         x, y = np.random.randint(0, 100, 2)
#         particle_list.append(Particle.TCR(x, y))
#
#
#     def force_func(x, y):
#         return x + y
#
#
#     parameters = [0, 1, 2, 3]
#
#     bd = BrownianDynamics(membrane_borders, parameters=parameters,
#                           force_function=force_func)
#
#     bd.run_brownian_dynamics(particle_list, 10, 0.1)
