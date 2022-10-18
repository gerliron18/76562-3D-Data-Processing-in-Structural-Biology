import os
import os.path as pth
import matplotlib.pyplot as plt
import moviepy.video.io.ImageSequenceClip
import numpy as np
import re
from PIL import Image

TIME_RE = re.compile("frame_(\d+).jpg")



def classify_points(moment):
    """
    Generate a dictionary represents the particles of one time slot of the
    simulation.
    :param moment:
    :return:
    """
    dict = {}

    locations_arr = moment[0]
    type_arr = moment[1]

    for i in range(len(type_arr)):
        if type_arr[i] not in dict:
            dict[type_arr[i]] = [[], []]
        dict[type_arr[i]][0].append(locations_arr[i][0])
        dict[type_arr[i]][1].append(locations_arr[i][1])

    return dict


def get_color(type):
    """
    Get the color of given particle type.
    :param type: The type of one particle
    :return: The color represent the particle type
    """
    if type == 'MHC':
        return 'blue'
    if type == 'TCR':
        return 'green'
    if type == 'CD45':
        return 'red'
    if type == 'LCK':
        return 'black'


def singlePlot(arr, size_x, size_y, output):
    """

    :param arr:
    :param size_x:
    :param size_y:
    :param output:
    :return:
    """
    filenum = 1

    for moment in arr:
        dict = classify_points(moment)

        x = np.ceil(size_x / 100)
        y = np.ceil(size_y / 100)
        fig = plt.figure(figsize=(x,y))
        # fig = plt.figure()
        ax = fig.add_subplot()
        for type in dict:
            pcolor = get_color(type)
            ax.plot(dict[type][0]-x/2, dict[type][1]+y/2, '.', color=pcolor, label=type, alpha=0.1)

        plt.axis('off')

        if pth.exists(os.path.join(output, "PNG",str(filenum) + ".png")):
            filenum += 1
        plt.savefig(os.path.join(output, "PNG", str(filenum) + ".png"))
        colorImage  = Image.open(os.path.join(output, "PNG", str(filenum) + ".png"))
        transposed  = colorImage.transpose(Image.ROTATE_270)
        transposed.save(os.path.join(output, "PNG", str(filenum) + ".png"))



def makeMovie(output):
    """
    Generates video using the simulation plots. Will make a mp4 format file
    under the 'MP4' folder showing frame by frame according to time.
    :param output: The path to the general output folder
    :return: The path to the generated video file
    """
    filename = "Movie"
    filenum = 1
    image_folder = os.path.join(output,"PNG")
    fps = 25

    image_files = [img.split(".")[0] for img in os.listdir(image_folder) if img.endswith(".png")]
    image_files.sort(key=int)
    new_image_files = [os.path.join(image_folder, s+".png") for s in image_files]


    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(new_image_files,
                                                                fps=fps)

    if pth.exists(os.path.join(output,'MP4', filename + str(filenum) + ".mp4")):
        filenum += 1

    clip.write_videofile(os.path.join(output,'MP4', filename + str(filenum) + ".mp4"))

    return str(os.path.join(output,'MP4', filename + str(filenum) + ".mp4"))
#
# def main():
#     makeMovie(r"H:\Study\university\Computational-Biology\Year 3\Semester B\3D_Structure\hackathon\yair_movies\frames\test")
#
#
# if __name__ == '__main__':
#
#     main()