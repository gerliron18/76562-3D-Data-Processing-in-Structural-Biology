import argparse
import re
import sys

import numpy as np

from Forces.choose_force_function import get_force_function, MC_coefficient
from MonteCarlo import MC, random_BD
from brownian_dynamics import BrownianDynamics
from gui import gen_png_mp4
from image_processing import extract_parameters_from_frame

TIME_RE = re.compile("frame_(\d+).jpg")


def get_args():
    parser = argparse.ArgumentParser(description='Run MC simulation on 2 images')
    parser.add_argument("--frame1", dest='frame1', type=str, required=False,
                        help='The path of the first frame')
    parser.add_argument('--frame2', dest="frame2", type=str, required=False,
                        help='The path of the second frame ')
    parser.add_argument('--output', dest="output", type=str, required=True,
                        help='The path of the output')
    parser.add_argument("--mc_iterations", dest="mt_iterations", type=int, default=1,
                        help="Number of mc iterations")
    parser.add_argument("--bd_deltaT", dest="bd_deltaT", type=float, default=0.01,
                        help="The deltaT for bd in seconds ")
    parser.add_argument("--bd_num_iterations", dest="bd_num_iterations", type=int, default=100,
                        help="Number of iteration to run bd ")
    parser.add_argument("--force_function", dest="force_function", type=str, default="location_interaction",
                        help="The name of the force function")

    parser.add_argument("--max_x", dest="max_x", default=900, type=int,
                        help="Max size of x in pixels ")

    parser.add_argument("--max_y", dest="max_y", default=900, type=int,
                        help="Max size of y in pixels ")

    args = parser.parse_args()
    return args


def get_configuration_from_frame(frame):
    """
    Extract the configuration from the frame
    :param frame: The frame object
    :return: the configuration relevant from the frame 
    """
    conf = extract_parameters_from_frame.get_configuration_from_frame(frame)
    frame_time = int(TIME_RE.findall(frame)[0])
    return conf, frame_time


def run_mc_simulation(frame1, frame2, mc_iterations, bd_deltaT, output, force_function_name):
    """
    Get two frames and run the MC simulation on it.
    Output should be the parameters created the best fit and a movie of the simulation
    :param frame1: The start frame
    :param frame2: The end frame
    :param mc_iterations: Number of iterations form mc
    :param bd_deltaT: dT for the brownian dynamics
    :param output: The output folder
    :param force_function_name: The name of the force functions
    :return:
    """
    start_conf, s_time = get_configuration_from_frame(frame1)
    end_conf, e_time = get_configuration_from_frame(frame2)
    num_of_bd_iter = (e_time - s_time) / bd_deltaT
    function_conf = get_force_function(force_function_name)

    list_of_ranges = [(p.min, p.max) for p in function_conf[1]]
    list_dim = [int((p.max - p.min) / p.step) for p in function_conf[1]]
    parameters_names = [p.get_name() for p in function_conf[1]]

    best_parameters = MC.run_mc(mc_iter=mc_iterations, list_of_ranges=list_of_ranges, list_dim=list_dim,
                                start_conf=start_conf, num_of_bd_iter=num_of_bd_iter, output=output,
                                bd_deltaT=bd_deltaT, force_function=function_conf[0],
                                result_pic_xy=BrownianDynamics.BrownianDynamics.get_location_vector(
                                    end_conf.get_particles()),
                                result_pic_types=BrownianDynamics.BrownianDynamics.get_type_vector(
                                    end_conf.get_particles()),
                                parameters_names=parameters_names)

    parameters_vec = dict(zip(parameters_names, best_parameters))
    output_path = run_bd_simulation_with_params(bd_dt=bd_deltaT, bd_num_iterations=num_of_bd_iter,
                                                conf=start_conf, force_function=function_conf[0],
                                                output=output,
                                                parameters_vec=parameters_vec)

    return (output_path, parameters_vec)


def get_random_parameter(p):
    """
    Get a random value for parameter based on it's range
    :param p: The parameter
    :type p: MC_coefficient
    :return: a value for the parameter
    """
    return np.random.randint(0, int(round((p.max - p.min) / p.step))) * p.step + p.min


def run_bd_simulation_based_on_frame(frame, bd_num_iterations, bd_dt, output, force_function_name):
    """
    Run a BD simulation and output the movie
    """
    conf, _ = get_configuration_from_frame(frame)
    return run_bd_simulation(bd_dt, bd_num_iterations, conf, force_function_name, output)


def run_bd_simulation(bd_dt, bd_num_iterations, conf, force_function_name, output):
    function_conf = get_force_function(force_function_name)
    parameters_vec = {p.get_name(): get_random_parameter(p) for p in function_conf[1]}
    return run_bd_simulation_with_params(bd_dt, bd_num_iterations, conf, function_conf[0], output,
                                         parameters_vec)


def run_bd_simulation_with_params(bd_dt, bd_num_iterations, conf, force_function, output, parameters_vec):
    bd = BrownianDynamics.BrownianDynamics(membrane_borders=conf.get_borders(),
                                           force_function=force_function, parameters_vec=parameters_vec)
    bd.run_brownian_dynamics(particle_list=conf.get_particles(), num_iter=bd_num_iterations, dt=bd_dt)
    gen_png_mp4.singlePlot(arr=bd.get_snapshots(), output=output, size_x=conf.get_borders()[0],
                           size_y=conf.get_borders()[1])
    return gen_png_mp4.makeMovie(output)


def run_bd_simulation_generate_conf(max_x, max_y, bd_num_iterations, bd_dt, output, force_function_name):
    """
    Run a BD simulation and output the movie
    """
    conf= random_BD.random_BD(max_x=max_x, max_y=max_y)
    return run_bd_simulation(bd_dt, bd_num_iterations, conf, force_function_name, output)


def main():
    args = get_args()

    # Going to run MC
    if args.frame2:
        run_mc_simulation(frame1=args.frame1, frame2=args.frame2, mc_iterations=args.mt_iterations,
                          bd_deltaT=args.bd_deltaT, output=args.output,
                          force_function_name=args.force_function)

    # Going to run bd based on frame
    elif args.frame1:
        run_bd_simulation_based_on_frame(frame=args.frame1, bd_num_iterations=args.bd_num_iterations,
                                         bd_dt=args.bd_deltaT, output=args.output,
                                         force_function_name=args.force_function)

    # Going to run bd based on frame size
    elif args.max_x and args.max_y:
        run_bd_simulation_generate_conf(bd_num_iterations=args.bd_num_terations,
                                        bd_dt=args.bd_deltaT, output=args.output,
                                        force_function_name=args.force_function,
                                        max_x=args.max_x, max_y=args.max_y)

    # Something wrong
    else:
        sys.exit("Invalid arguments")


if __name__ == '__main__':
    main()
