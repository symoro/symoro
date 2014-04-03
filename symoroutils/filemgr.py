#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Perform file management operations for the SYMORO package."""


import os


SYMORO_ROBOTS_FOLDER = "symoro-robots"


def get_base_path():
    """
    Return the base path for storing all SYMORO robot files.

    Returns:
        A string specifying the base folder path.
    """
    home_folder = os.path.expanduser("~")
    return os.path.join(home_folder, SYMORO_ROBOTS_FOLDER)


def make_folders(folder_path):
    """
    Check if a specified folder path exists and create the folder path
    if it does not exist.

    Args:
        folder_path: The folder path (string) to check and create.
    """
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


def get_folder_path(robot_name):
    """
    Return the folder path to store the robot data. Also create the 
    folders if they are not already present.
    
    Args:
        robot_name: The name of the robot.

    Returns:
        A string specifying the folder path.
    """
    folder_path = os.path.join(get_base_path(), robot_name)
    make_folders(folder_path)
    return folder_path


def make_file_path(robot, ext=None):
    """
    Create the file path with the appropriate extension appended to 
    the file name using an underscore.

    Args:
        robot: An instance of the `Robot` class.
        ext: The extension (string) that is to be appended to the file
            name with an underscore.

    Returns:
        The file path (string) created.
    """
    if ext is None:
        fname = '%s.par' % robot.name
    else:
        fname = '%s_%s.txt' % (robot.name, ext)
    file_path = os.path.join(robot.directory, fname)
    make_folders(robot.directory)
    return file_path


