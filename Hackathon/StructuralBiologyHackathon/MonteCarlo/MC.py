import numpy as np

from brownian_dynamics import BrownianDynamics
from brownian_dynamics.BDConfiguration import BDConfiguration
from gui import gen_png_mp4

window = 30
SIZE = 930


def compare_bd_result(wanted_result, bd_result, n=SIZE, m=SIZE):
    """
    loss function for the bd result. compare the result to the wanted_result
    from the given picture
    :param wanted_result - grid SIZE*SIZE represent the elements on the given picture
    :param bd_result - grid SIZE*SIZE represent the elements on the result of bd run
    :return: loss value
    """
    loss = 0
    for i in range(0, n, window):  # go over the grid in window*window size
        for j in range(0, m, window):
            num_tcr_first, num_tcr_second = 0, 0  # number of TCRs in window for each grid
            num_cd_first, num_cd_second = 0, 0  # number of cds in window for each grid
            num_sum_f, num_sum_s = 0, 0  # number of elements in window for each grid
            # in each window we calculate the difference between the % of tcr and cd45
            for k in range(i, min(i + window, n)):
                for l in range(j, min(j + window, m)):
                    if wanted_result[k, l]:
                        num_sum_f += 1
                        if wanted_result[k, l] == 1:
                            num_tcr_first += 1
                        else:
                            num_cd_first += 1
                    if bd_result[k, l]:
                        num_sum_s += 1
                        if bd_result[k, l] == 1:
                            num_tcr_second += 1
                        else:
                            num_cd_second += 1
            # calculate %s
            p_tcr_first, p_cd_first, p_tcr_second, p_cd_second = 0, 0, 0, 0
            if num_sum_f:
                p_tcr_first = num_tcr_first / num_sum_f
                p_cd_first = num_cd_first / num_sum_f
            if num_sum_s:
                p_tcr_second = num_tcr_second / num_sum_s
                p_cd_second = num_cd_second / num_sum_s
            # update loss
            loss += abs(p_tcr_first - p_tcr_second) + abs(p_cd_first - p_cd_second)
    return loss


def crate_grid_elements(coordinates, types, size=SIZE):
    """ Create grid SIZE*SIZE with elements on it
    :param coordinates - 2*n (n number of elements) with [x,y] coordinates
    :param types - n sized array with the type of the elements
    """
    grid_i = np.zeros((size, size))
    for i in range(len(coordinates)):
        x = int(round(coordinates[i][0]))
        y = int(round(coordinates[i][1]))
        if 0<=x<size and 0<=y<size:
            type = types[i]
            if type == "CD45":
                grid_i[(x, y)] = 2
            if type == "TCR":
                grid_i[(x, y)] = 1
    return grid_i


def is_valid(c, list_dim):
    """ Return True if c is a valid coordinate on an (n) grid """
    if len(c) != len(list_dim):
        return False
    for i in range(len(c)):
        if c[i] < 0 or c[i] >= list_dim[i]:
            return False
    return True


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


def calculate_vector_of_parameters(c, list_dim, list_of_ranges):
    """ indexes vector c convert to vector of values for BD
    :param c: indexes vector of grid
    :param list_dim: the shape of the grid
    :param list_of_ranges: ranges for each parameter
    :return: vector of values for BD """
    v = []  # parameter's vector for BD for c configuration
    for i in range(len(list_dim)):
        w = c[i] / (list_dim[i] - 1)
        v.append((list_of_ranges[i][1] - list_of_ranges[i][0]) * w + list_of_ranges[i][0])
    return v


def E(grid_configuration, list_dim, start_conf, list_of_ranges, num_bd_iterations, bd_deltaT,
      wanted_result, parameters_vector, force_function,parameters_names, size=SIZE):
    """
    Calculate E-loss for c set of parameters
    :param grid_configuration: configuration on the grid
    :param list_dim: the shape of the grid
    :param start_conf: start_conf for BD simulation
    :param list_of_ranges: ranges for each parameter
    :param num_bd_iterations: number of iterations for bd
    :param bd_deltaT: the time delta for bd
    :param wanted_result - grid SIZE*SIZE represent the elements on the given picture
    :return: loss value
    """
    assert (len(grid_configuration) == len(list_dim))
    parameters_values = calculate_vector_of_parameters(grid_configuration, list_dim, list_of_ranges)
    parameters_vector = dict(zip(parameters_names,parameters_values))
    bd = BrownianDynamics.BrownianDynamics(force_function=force_function, parameters_vec=parameters_vector,
                                           membrane_borders=start_conf.get_borders())
    bd.run_brownian_dynamics(num_iter=num_bd_iterations, dt=bd_deltaT,
                             particle_list=start_conf.get_particles())
    final_snapshot = bd.get_final_snapshot()
    return loss_function(final_snapshot, wanted_result, size)


def loss_function(snapshot, wanted_result, size):
    coordinates, types = snapshot
    bd_result = crate_grid_elements(coordinates, types)
    return compare_bd_result(wanted_result, bd_result, size, size)


def get_neighbours(c, n):
    """ get up/down/left/right neighbours on an vector (n) grid"""
    assert (is_valid(c, n))
    ret_value = []
    for i in range(len(c)):
        if c[i] > 0:
            copy = c.copy()
            copy[i] -= 1
            ret_value.append(copy)
        if c[i] < n[i] - 1:
            copy = c.copy()
            copy[i] += 1
            ret_value.append(copy)
    return ret_value


def calculate_dE(current_grid_conf, next_grid_conf, grid_loss_values, list_dim, start_conf, list_of_ranges,
                 bd_iterations, bd_deltaT, wanted_result, parameters_vector, force_function,
                 parameters_names, size=SIZE):
    """
    calculates dE(dLOSS)
    :param current_grid_conf: current configuration on the grid
    :param next_grid_conf: configuration on the grid we want to go over
    :param grid_loss_values: grid of the saves E(loss) values
    :param list_dim: the shape of the grid
    :param start_conf: start_conf for BD simulation
    :param list_of_ranges: ranges for each parameter
    :param bd_iterations: number of iterations for bd
    :param bd_deltaT: time delta for the bd
    :param wanted_result - grid SIZE*SIZE represent the elements on the given picture
    :return:
    """
    # take saved values if zeros calculate
    e1 = grid_loss_values.item(tuple(current_grid_conf))
    e2 = grid_loss_values.item(tuple(next_grid_conf))
    # if no saved values calculate them and save in grid_loss_values
    if e1 == 0.0:
        new_v = E(grid_configuration=current_grid_conf, list_dim=list_dim, start_conf=start_conf,
                  list_of_ranges=list_of_ranges, bd_deltaT=bd_deltaT, num_bd_iterations=bd_iterations,
                  wanted_result=wanted_result, parameters_vector=parameters_vector,
                  force_function=force_function, parameters_names=parameters_names, size=size)
        grid_loss_values[tuple(current_grid_conf)] = new_v
        e1 = new_v
    if e2 == 0.0:
        new_v = E(grid_configuration=next_grid_conf, list_dim=list_dim, start_conf=start_conf,
                  list_of_ranges=list_of_ranges, bd_deltaT=bd_deltaT, num_bd_iterations=bd_iterations,
                  wanted_result=wanted_result, parameters_vector=parameters_vector,
                  force_function=force_function,parameters_names=parameters_names, size=size)
        grid_loss_values[tuple(next_grid_conf)] = new_v
        e2 = new_v
    # return dE
    return e2 - e1


def markov_chain_mc(mc_iter, n_iter, grid, grid_loss_values, kT, list_dim, start_conf, list_of_ranges,
                    bd_deltaT, wanted_result, force_function,parameters_names, size=SIZE):
    """
    run a markov chain mc
    :param mc_iter: Number of iteration for the mc
    :param n_iter: number of iterations for the bd
    :param grid: the configuration space grid
    :param grid_loss_values: grid of the saves E(loss) values
    :param kT: the kT value
    :param list_dim: the shape of the grid
    :param start_conf: start configuration for the BD simulation
    :param list_of_ranges: ranges for each parameter
    :param bd_deltaT: a time delta for the bd
    :param wanted_result - grid SIZE*SIZE represent the elements on the given picture
    :return:
    """
    # random start configuration
    c = []
    for j in range(len(list_dim)):
        c.extend(np.random.randint(0, list_dim[j], 1))

    # Go over the grid mc_iter times
    for _ in range(mc_iter):
        cur_neighbours = get_neighbours(c, list_dim)

        next_c = cur_neighbours[np.random.randint(0, len(cur_neighbours))]  # chose neighbour
        p_forward = 1 / len(cur_neighbours)

        next_neighbours = get_neighbours(next_c, list_dim)  # neighbours of next
        p_backward = 1 / len(next_neighbours)

        # calculate the dE to the chosen neighbour
        dE = calculate_dE(
            current_grid_conf=c, next_grid_conf=next_c, grid_loss_values=grid_loss_values,
            list_dim=list_dim, list_of_ranges=list_of_ranges,
            bd_iterations=n_iter, bd_deltaT=bd_deltaT, start_conf=start_conf,
            wanted_result=wanted_result, force_function=force_function, parameters_vector=c,
            parameters_names=parameters_names, size=size)

        p = get_p_accept_metropolis(dE, kT, p_forward, p_backward)
        if np.random.choice(np.array([0, 1]), p=[1 - p, p]):  # chose if move or not
            # according to probability p
            c = next_c
        grid[tuple(c)] += 1  # update the grid

    return grid


def run_mc(mc_iter, num_of_bd_iter, list_of_ranges, list_dim, start_conf, output, force_function,
           bd_deltaT, result_pic_xy, result_pic_types, parameters_names, kT=1, num_random_starts=1):
    """
    :param mc_iter: number of iterations for mc
    :param num_of_bd_iter: number of iteration for the bd
    :param list_of_ranges: list of ranges for each parameter [[s1,e1], [s2,s2], ..., [sn,sn]]
    :param list_dim: tuple of dims for each parameter (d1, d2,..., dn) --> this is the step size
    :param start_conf: the start configuration, (borders, list with all the particles)
    :type start_conf: BDConfiguration
    :param output: the path for the output file
    :param bd_deltaT: deltaT for each bd run
    :param kT: the kd for the MC
    :param result_pic_xy - 2*n (n- num of elements in the result pic) with [x,y] coordinates
    :param result_pic_types - array with the type of the elements in the result pic
    :param num_random_starts - number of random start in the simulation
    :return:
    """
    # create grids
    grid = np.zeros(list_dim)  # for visiting
    grid_loss_values = np.zeros(list_dim)  # for saving E values
    x = start_conf.get_borders()[0]
    y = start_conf.get_borders()[1]
    size = max(x, y)
    wanted_result = crate_grid_elements(result_pic_xy, result_pic_types, SIZE)

    for _ in range(num_random_starts):  # multiple starts
        markov_chain_mc(mc_iter=mc_iter, n_iter=num_of_bd_iter, grid=grid,
                        grid_loss_values=grid_loss_values,
                        kT=kT, list_dim=list_dim, start_conf=start_conf,
                        list_of_ranges=list_of_ranges, bd_deltaT=bd_deltaT,
                        wanted_result=wanted_result, force_function=force_function,
                        parameters_names=parameters_names, size=SIZE)

    # return best set of parameters
    index_best = np.where(grid == grid.max())
    c_best = []
    for i in index_best:
        c_best.append(list(i)[0])
    v_best = calculate_vector_of_parameters(c_best, list_dim, list_of_ranges)
    return v_best
    

