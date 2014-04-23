# -*- coding: utf-8 -*-


"""
This module contains the Screw6 data structure.
"""


from sympy import zeros
from sympy import ShapeError


class Screw6(object):
    """
    Data structure:
        Represent the data structure (base class) to hold a 6x6 matrix
        which in turn contains four 3x3 matrices.
    """
    def __init__(self, *args, **kwargs):
        """Constructor period."""
        self._val = zeros(6, 6)

    @property
    def val(self):
        """
        Get current value.

        Returns:
            A 6x6 Matrix with the current value
        """
        return self._val

    @val.setter
    def val(self, value):
        """
        Set the current value.

        Args:
            value: A 6x6 Matrix
        """
        if value.rows != 6 or value.cols != 6:
            raise ShapeError("Matrix size has to be 6x6.")
        self._val = value

    @property
    def topleft(self):
        """
        Get the top-left part of the 6x6 matrix.

        Returns:
            A 3x3 Matrix.
        """
        return self._val[0:3, 0:3]

    @property
    def topright(self):
        """
        Get the top-right part of the 6x6 matrix.

        Returns:
            A 3x3 Matrix.
        """
        return self._val[0:3, 3:6]

    @property
    def botleft(self):
        """
        Get the bottom-left part of the 6x6 matrix.

        Returns:
            A 3x3 Matrix.
        """
        return self._val[3:6, 0:3]

    @property
    def botright(self):
        """
        Get the bottom-right part of the 6x6 matrix.

        Returns:
            A 3x3 Matrix.
        """
        return self._val[3:6, 3:6]

    @topleft.setter
    def topleft(self, value):
        """
        Set the top-left part of the 6x6 matrix.

        Args:
            value: A 3x3 Matrix - top-left value.
        """
        if value.rows != 3 or value.cols != 3:
            raise ShapeError("Top-left value size has to be 3x3.")
        self._val[0:3, 0:3] = value

    @topright.setter
    def topright(self, value):
        """
        Set the top-right part of the 6x6 matrix.

        Args:
            value: A 3x3 Matrix - top-right value.
        """
        if value.rows != 3 or value.cols != 3:
            raise ShapeError("Top-right value size has to be 3x3.")
        self._val[0:3, 3:6] = value

    @botleft.setter
    def botleft(self, value):
        """
        Set the bottom-left part of the 6x6 matrix.

        Args:
            value: A 3x3 Matrix - bottom-left value.
        """
        if value.rows != 3 or value.cols != 3:
            raise ShapeError("Bottom-left value size has to be 3x3.")
        self._val[3:6, 0:3] = value

    @botright.setter
    def botright(self, value):
        """
        Set the bottom-right part of the 6x6 matrix.

        Args:
            value: A 3x3 Matrix - bottom-right value.
        """
        if value.rows != 3 or value.cols != 3:
            raise ShapeError("Bottom-right value size has to be 3x3.")
        self._val[3:6, 3:6] = value


