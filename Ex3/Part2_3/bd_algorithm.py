import argparse
import math
import sys

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sbn
from scipy.special import logsumexp

MIN_POS = -0.5
MAX_POS = 4.5
FORCE_VECTOR = np.array([-1, -0.5])


def get_cmdline_parser():
    """
    Get the parmaters from the cmdline
    :return:
    """
    parser = argparse.ArgumentParser(
        description='Run BD on a discrete n x n configuration space.')
    parser.add_argument('n', type=int, default=5, nargs='?',
                        help='number of rows/columns of grid')
    parser.add_argument('m', type=int, default=1000000, nargs='?',
                        help='number of iterations of MCMC optimization')
    parser.add_argument('kT', type=float, default=1, nargs='?',
                        help='kT - the denominator for the metropolis criterion'
                             ' (Boltzmann constant times temperature)')
    parser.add_argument('dT', type=float, default=0.1, nargs='?',
                        help='dT - the time step')
    parser.add_argument('D', type=float, default=1, nargs='?',
                        help='the diffusion coefficient')
    return parser


def get_random_float(low_limit, high_limit):
    """
    Get random float between low_limit and high_limit
    :param low_limit: The minimum number
    :param high_limit: The maximum number
    :return: A random number x s.t low_limit <=x<=high_limit
    """
    r = np.random.random()  # from 0 to 1
    r = r * (high_limit - low_limit)  # 0 to 5
    r = r + low_limit  # -.5 to 4.5
    return r


def get_random_2D_vector():
    """
    Create a 2D random vector
    :return:
    """
    return np.array([np.random.normal(0, 1), np.random.normal(0, 1)])


def run_brownian_dynamics(num_of_steps, delta_t, kbt, diffusion):
    """Run a BD algorithm"""
    current = np.array([get_random_float(MIN_POS, MAX_POS), get_random_float(MIN_POS, MAX_POS)])
    movement_list = [current]

    for i in range(num_of_steps):
        current = current + delta_t / kbt * diffusion * FORCE_VECTOR + math.sqrt(
            6 * diffusion * delta_t) * get_random_2D_vector()

        current[0] = max(MIN_POS, current[0])
        current[0] = min(MAX_POS, current[0])
        current[1] = max(MIN_POS, current[1])
        current[1] = min(MAX_POS, current[1])

        movement_list.append(current)

    return movement_list


def get_grid(movement_list, table_size):
    """
    Construct a grid with specific table size based on the movement list of configuration
    :param movement_list: A list of configurations
    :param table_size: The grid size
    :return: The grid, each cell is the number of times we visited the cell
    """
    grid = np.zeros(shape=(table_size, table_size))
    for c in movement_list:
        rounded_c = np.floor(c + 0.5).astype(np.int)
        rounded_c[0] = min(4, rounded_c[0])
        rounded_c[1] = min(4, rounded_c[1])

        grid[rounded_c[0], rounded_c[1]] += 1

    return grid


# Internal to create the heatmap
def get_heatmap(grid, kt_value):
    plt.figure()
    m = np.sum(grid)  # stats is an nxn matrix
    sbn.heatmap(grid / m, annot=True, fmt=".2f", linewidths=.5)
    plt.xlabel("Column index")
    plt.ylabel("Row index")
    plt.title("Heatmap for BD for kt=%s" % (kt_value))
    plt.savefig("BD_heatmap_kt_%s.png" % (kt_value))


# Internal to create the heatmap for the log1p
def get_log_heatmap(grid, kt_value):
    grid = kt_value * np.max(np.log1p((grid))) - np.log1p(grid)

    plt.figure()
    # m = logsumexp(grid)  # stats is an nxn matrix
    sbn.heatmap(grid, annot=True, fmt=".2f", linewidths=.5)
    plt.xlabel("Column index")
    plt.ylabel("Row index")
    plt.title("Heatmap(log1p) for BD KT=%s" % kt_value)
    plt.savefig("normalized_BD_heatmap_kt_%s.png" % (kt_value))


def main():
    """
    Validate the user input and run the mcmc algorithm
    :return:
    """
    args = get_cmdline_parser().parse_args()
    table_size = args.n
    num_iteration = args.m
    kt_value = args.kT
    diffusion_coeff = args.D
    dT = args.dT

    if table_size <= 0 or num_iteration <= 0:
        sys.exit("Some of the input values are invalid")

    movement_list = run_brownian_dynamics(num_iteration, dT, kt_value, diffusion_coeff)
    grid = get_grid(movement_list,table_size)
    # get_heatmap(grid, kt_value)
    get_log_heatmap(grid, kt_value)


if __name__ == '__main__':
    main()
