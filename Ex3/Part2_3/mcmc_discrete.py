import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sbn


from scipy.special import logsumexp


def get_cmdline_parser():
    """
    Get the parmaters from the cmdline
    :return:
    """
    parser = argparse.ArgumentParser(
        description='Run MCMC on a discrete n x n configuration space.')
    parser.add_argument('n', type=int, default=5, nargs='?',
                        help='number of rows/columns of grid')
    parser.add_argument('m', type=int, default=1000, nargs='?',
                        help='number of iterations of MCMC optimization')
    parser.add_argument('kT', type=float, default=1.0, nargs='?',
                        help='kT - the denominator for the metropolis criterion'
                             ' (Boltzmann constant times temperature)')
    return parser


def is_valid(c, n):
    """
    Return True if c is a valid 2-D coordinate on an n x n grid
    with 0-based indices
    """
    if len(c) != 2:
        return False
    return (c[0] >= 0 and c[1] >= 0 and c[0] < n and c[1] < n)


def get_p_accept_metropolis(dE, kT, p_forward, p_backward):
    """
    return the probability to accept the metropolis criteria
    for the specified conditions
    dE - change in energy from current to proposed configuration
    kT - the factor of Boltzmann constant (kB) and the temperature (T)
    p_forward - probability to propose a move from current to proposed configuration
    p_backward - probability to propose a move from proposed to current configuration
    """
    p = np.exp(-dE / kT) * p_backward / p_forward
    return min(p, 1.0)


def E(c):
    """
    Calculate the energy of a position in the grid
    :param c: The position
    :return: The energy of this position
    """
    assert (len(c) == 2)
    return 1.0 * c[0] + 0.5 * c[1]


def get_neighbours(c, n):
    ''' get up/down/left/right neighbours on an n x n grid with 0-based indices'''
    assert (is_valid(c, n))
    ret_value = []
    if c[0] > 0:
        ret_value.append((c[0] - 1, c[1]))
    if c[0] < n - 1:
        ret_value.append((c[0] + 1, c[1]))
    if c[1] > 0:
        ret_value.append((c[0], c[1] - 1))
    if c[1] < n - 1:
        ret_value.append((c[0], c[1] + 1))
    return ret_value


def mcmc(table_size, num_iteration, kt_value):
    """
    Run an mcmc algorithm given the parmaters
    :param table_size: the grid size
    :param num_iteration: the number of iteration to run the algorithm
    :param kt_value: the KbT value
    :return: the grid of the mcmc run, each cell is a counter of the number of iteration spend in this cell
    """
    grid = np.zeros(shape=(table_size, table_size))
    start_configuration = (np.random.randint(0, table_size), np.random.randint(0, table_size))
    current = start_configuration

    for i in range(num_iteration):
        grid[current[0], current[1]] += 1

        neighbours = get_neighbours(current, table_size)
        next_neighbour = neighbours[np.random.choice(len(neighbours), 1)[0]]
        p_forward = 1 / len(neighbours)
        p_backward = 1 / len(get_neighbours(next_neighbour, table_size))

        dE = E(next_neighbour) - E(current)
        p_accept = get_p_accept_metropolis(dE, kt_value, p_forward, p_backward)

        if p_accept > np.random.random():
            current = next_neighbour

    return grid


def print_grid(grid):
    """
    Print the grid as requested in the pdf
    :param grid: the grid
    :return: None
    """
    rounded_grid = grid.astype(np.int)
    str_grid = rounded_grid.astype(np.str)
    for i in range(str_grid.shape[0]):
        print(",".join(list(str_grid[i, :].astype(np.int).astype(np.str))))


def main():
    """
    Validate the user input and run the mcmc algorithm
    :return:
    """
    args = get_cmdline_parser().parse_args()
    table_size = args.n
    num_iteration = args.m
    kt_value = args.kT

    if table_size <= 0 or num_iteration <= 0:
        sys.exit("Some of the input values are invalid")

    grid = mcmc(table_size, num_iteration, kt_value)
    assert np.sum(grid) == num_iteration
    print_grid(grid)

# Interal to create the heatmap
def get_heatmap(table_size, num_iteration, kt_value, question_name):
    grid = mcmc(table_size, num_iteration, kt_value)
    plt.figure()
    m = np.sum(grid)  # stats is an nxn matrix
    sbn.heatmap(grid / m, annot=True, fmt=".2f", linewidths=.5)
    plt.xlabel("Column index")
    plt.ylabel("Row index")
    plt.title("Heatmap for MCMC for %s iteration kt=%s" % (num_iteration, kt_value))
    plt.savefig("%s_%s_iteration_%s_kt.png" % (question_name, num_iteration, kt_value))


# Interal to create the heatmap for the log1p
def get_log_heatmap(table_size, num_iteration, kt_value, question_name):
    grid = mcmc(table_size, num_iteration, kt_value)
    grid = kt_value * np.max(np.log1p((grid))) - np.log1p(grid)

    plt.figure()
    m = logsumexp(grid)  # stats is an nxn matrix
    sbn.heatmap(grid, annot=True, fmt=".2f", linewidths=.5)
    plt.xlabel("Column index")
    plt.ylabel("Row index")
    plt.title("Heatmap(log1p) for MCMC KT=%s" % kt_value)
    plt.savefig("%s_%s_iteration_%s_kt.png" % (question_name, num_iteration, kt_value))


if __name__ == '__main__':
    main()
    # get_heatmap(table_size=5, num_iteration=10, kt_value=1.0, question_name="q4")
    # get_heatmap(table_size=5, num_iteration=100, kt_value=1.0, question_name="q4")
    # get_heatmap(table_size=5, num_iteration=500000, kt_value=1.0, question_name="q4")
    # get_log_heatmap(table_size=5, num_iteration=500000, kt_value=1.0, question_name="q6")
    # get_heatmap(table_size=5, num_iteration=500000, kt_value=0.1, question_name="q9")
    # get_heatmap(table_size=5, num_iteration=500000, kt_value=1, question_name="q9")
    # get_heatmap(table_size=5, num_iteration=500000, kt_value=2, question_name="q9")
    # get_heatmap(table_size=5, num_iteration=500000, kt_value=50, question_name="q9")
    # get_log_heatmap(table_size=5, num_iteration=500000, kt_value=0.1, question_name="q10")
    # get_log_heatmap(table_size=5, num_iteration=500000, kt_value=2, question_name="q10")
    # get_log_heatmap(table_size=5, num_iteration=500000, kt_value=50, question_name="q10")

