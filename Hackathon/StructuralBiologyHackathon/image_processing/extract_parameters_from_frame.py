"""
given a frame, extract the following:
    1. grid size
    2. places of each particle

"""
import argparse

import cv2
import numpy as np

from brownian_dynamics import Particle, BDConfiguration


def parse_args():
    parser = argparse.ArgumentParser(description='Covert video to frames with time')
    parser.add_argument("--frame", dest='frame_path', type=str, required=True,
                        help='The path of the video')

    args = parser.parse_args()
    return args


def get_configuration_from_frame(frame_path):
    """
    Get a configuration object from a frame, currently contains:
        1. The borders of the frame
        2. The particles in the frames
    :param frame_path: The frame path
    :return: An BDConfiguration objects with the required data
    """
    frame = cv2.imread(frame_path)
    particles = get_particles_from_frame(frame)
    borders = get_membrane_borders(frame)
    return BDConfiguration.BDConfiguration(particles=particles, borders=borders)


def get_membrane_by_frame(frame):
    """
    This function is one of the possible ways to calculate the membrane borders.
    We say that the frame is the borders
    :param frame: A frame object
    :return: The size of the frame
    """
    return frame.shape[0], frame.shape[1]


def get_membrane_borders(frame):
    """
    Extract the membrane borders in whichever way we want
    :param frame: The frame object
    :return: The membrane borders, in whatever way the user decided to calculate
    """
    # This can be change as what required
    return get_membrane_by_frame(frame)


def get_one_pixel_particles(tcr_zip, cd45_zip):
    """
    This is just one way to do it
    Convert the green and red zip to particles, we use one pixel to each particle
    :param green_zip: A zip with all x,y location for the TCR
    :param red_zip: A zip with all x,y location for the CD45
    :return: the particles
    """
    particles = []
    for p in tcr_zip:
        particles.append(Particle.TCR(*p))

    for p in cd45_zip:
        particles.append(Particle.CD45(*p))

    return particles


def get_particles_from_frame(frame):
    """
    Get all the particles from the frame
    Currently support TCR - green, CD45 - red
    :param frame: The frame object
    :return: A list of all the particles in the form of the object
    """
    rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)

    red_mask = cv2.inRange(rgb, (100, 0, 0), (255, 60, 60))
    green_mask = cv2.inRange(rgb, (0, 200, 33), (255, 255, 123))

    green_x, green_y = np.nonzero(green_mask > 0)
    red_x, red_y = np.nonzero(red_mask > 0)
    green_zip = zip(green_x, green_y)
    red_zip = zip(red_x, red_y)

    # Convert the pixels to particles, that can be done by several ways
    particles = get_one_pixel_particles(green_zip, red_zip)

    return particles

#
# def main():
#     args = parse_args()
#     membrane_borders = get_membrane_borders(args.frame_path)
#     particles = get_particles_from_frame(args.frame_path)
#     force_matrix = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]]
#     bd1 = bd_simulation.BrownianDynamics(particle_list=particles, force_matrix=force_matrix,
#                                          membrane_borders=membrane_borders)
#     bd1.run_brownian_dynamics(10)
#
#
# if __name__ == '__main__':
#     main()
