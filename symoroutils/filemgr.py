# -*- coding: utf-8 -*-


"""Perform file management operations for the SYMORO package."""


import os


SYMORO_ROBOTS_FOLDER = "symoro-robots"


def get_base_path(base_folder=SYMORO_ROBOTS_FOLDER):
    """
    Return the base path for storing all SYMORO robot files.

    Returns:
        A string specifying the base folder path.
    """
    home_folder = os.path.expanduser("~")
    return os.path.join(home_folder, base_folder)


def get_clean_name(name, char='-'):
    """
    Return a string that is lowercase and all whitespaces are replaced
    by a specified character.

    Args:
        name: The string to be cleaned up.
        char: The character to replace all whitespaces. The default
            character is "-" (hyphen).

    Returns:
        A string that is fully lowercase and all whitespaces replaced by
        the specified character.

    >>> get_clean_name('Random teXt')
    'random-text'
    >>> get_clean_name('Random teXt', '#')
    'random#text'
    """
    return name.lower().replace(' ', char)


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
        robot_name: The name of the robot (string).

    Returns:
        A string specifying the folder path.
    """
    robot_name = get_clean_name(robot_name)
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
        fname = '%s.par' % get_clean_name(robot.name)
    else:
        fname = '%s_%s.txt' % (get_clean_name(robot.name), ext)
    file_path = os.path.join(robot.directory, fname)
    make_folders(robot.directory)
    return file_path


