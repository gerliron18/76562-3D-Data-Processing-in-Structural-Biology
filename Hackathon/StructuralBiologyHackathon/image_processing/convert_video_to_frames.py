import argparse
import os
import pathlib

import cv2
import pytesseract

pytesseract.pytesseract.tesseract_cmd = r'C:\Program Files\Tesseract-OCR\tesseract.exe'

import re

TIME_RE = re.compile("(\d+) seconds.*")


def parse_args():
    parser = argparse.ArgumentParser(description='Covert video to frames with time')
    parser.add_argument("--video", dest='video_path', type=str, required=True,
                        help='The path of the video')
    parser.add_argument('--output', dest="output", type=str, required=True,
                        help='The output path')

    args = parser.parse_args()
    return args


def extract_time_from_frame(image):
    """
    Extract time from the frame, return -1 if no time was found
    based on the re "\d+ seconds"
    :param image: the image as np obj
    :return: the time of the frame or -1 if none was found
    """
    text = pytesseract.image_to_string(image)
    frame_time = TIME_RE.findall(text)
    if len(frame_time) == 0:
        return -1

    return int(frame_time[0])


def get_frame_time(image, previous_time):
    """
    Try and get the time from the frame in multiple ways
    :param image: The image as np object
    :param previous_time: The previous time to validate we get something that make sense
    :return: The time of the frame or -1 if something want wrong:
        1. The time is before the previous time which make no sense
        2. We couldn't find any time in the frame and inverted frame
    """
    frame_time = extract_time_from_frame(image)

    # # Note: this didn't seems to improve and just take time
    # if frame_time == -1:
    #     converted_image = image - 255
    #     frame_time = extract_time_from_frame(converted_image)
    #     if frame_time != -1:
    #         print("yay")

    return frame_time if frame_time > previous_time else -1


def video_to_frames(video_path, output_path):
    if not os.path.exists(output_path):
        pathlib.Path.mkdir(output_path, parents=True, exist_ok=True)

    vidcap = cv2.VideoCapture(video_path)
    success, image = vidcap.read()
    frame_time = 0
    while success:
        temp_frame_time = get_frame_time(image, frame_time)
        if temp_frame_time == -1:
            success, image = vidcap.read()
            continue
        frame_time = temp_frame_time
        frame_path = os.path.join(output_path, "frame_%d.jpg" % frame_time)
        cv2.imwrite(frame_path, image)
        success, image = vidcap.read()


def main():
    args = parse_args()
    video_to_frames(args.video_path, args.output)


if __name__ == '__main__':
    main()
