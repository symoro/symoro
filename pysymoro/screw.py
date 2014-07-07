# -*- coding: utf-8 -*-


"""
This module contains the Screw data structure.
"""


from sympy import zeros
from sympy import ShapeError


class Screw(object):
    """
    Data structure:
        Represent the screw notation - a 6x1 vector in which the first
        three rows represent the linear term and the last three rows
        represent the angular term.
        This class basically adds a lot of constraints to the Matrix
        class of sympy - atleast tries to :)
    """
    def __init__(self, lin=zeros(3, 1), ang=zeros(3, 1)):
        """
        Constructor period.

        Args:
            lin: A 3x1 matrix - linear term set to 0 by default.
            ang: A 3x1 matrix - angular term set to 0 by default.
        """
        self._val = zeros(6, 1)
        if lin.rows != 3 or lin.cols != 1:
            raise ShapeError("Linear term matrix size has to be 3x1.")
        elif ang.rows != 3 or ang.cols != 1:
            raise ShapeError("Angular term matrix size has to be 3x1.")
        self._val[0:3, 0] = lin
        self._val[3:6, 0] = ang

    def __str__(self):
        row_format = '[' + ('{:} ; ' * 5) + ('{:}') + ']'
        str_format = row_format.format(*(
            str(self._val[0]), str(self._val[1]), str(self._val[2]),
            str(self._val[3]), str(self._val[4]), str(self._val[5])
        ))
        return str_format

    def __repr__(self):
        repr_format = 'Screw({0})'.format(str(self))
        return repr_format

    @property
    def val(self):
        """
        Get the current value.

        Returns:
            A 6x1 Matrix (column vector) with the current value.
        """
        return self._val

    @val.setter
    def val(self, value):
        """
        Set the current value.

        Args:
            value: A 6x1 Matrix
        """
        if value.rows != 6 or value.cols != 1:
            raise ShapeError("Matrix size has to be 6x1.")
        self._val = value

    @property
    def lin(self):
        """
        Get the linear term value.

        Returns:
            A 3x1 Matrix with the current linear term value.
        """
        return self._val[0:3, 0]

    @lin.setter
    def lin(self, value):
        """
        Set the linear term value.

        Args:
            value: A 3x1 Matrix - linear term
        """
        if value.rows != 3 or value.cols != 1:
            raise ShapeError("Linear term matrix size has to be 3x1.")
        self._val[0:3, 0] = value

    @property
    def ang(self):
        """
        Get the angular term value.

        Returns:
            A 3x1 Matrix with the current angular term value.
        """
        return self._val[3:6, 0]

    @ang.setter
    def ang(self, value):
        """
        Set the linear term value.

        Args:
            value: A 3x1 Matrix - angular term
        """
        if value.rows != 3 or value.cols != 1:
            raise ShapeError("Angular term matrix size has to be 3x1.")
        self._val[3:6, 0] = value

    def __eq__(self, other):
        """Check equality between two instances of Screw."""
        if type(self) != type(other):
            raise ValueError(
                "Unable to compare %s with Screw type." % str(type(other))
            )
        return self.val == other.val

    def __ne__(self, other):
        """Check non-equality between two instances of Screw."""
        return not self == other


