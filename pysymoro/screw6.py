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
        """
        Constructor period.

        Usage:
        >>> # initialise to 0 by default
        Screw6()
        >>> # initialise to a given 6x6 matrix
        Screw6(<value>)
        >>> # intiialise each of the 4 sub-matrices individually
        Screw6(<top-left>, <top-right>, <bottom-left>, <bottom-right>)
        >>> # initialise using keywords
        Screw6(value=<value>)
        Screw6(
            tl=<top-left>, tr=<top-right>,
            bl=<bottom-left>, br=<bottom-right>
        )
        """
        self._val = zeros(6, 6)
        if len(args) == 1:
            self.val = args[0]
        elif len(args) == 4:
            self.topleft = args[0]
            self.topright = args[1]
            self.botleft = args[2]
            self.botright = args[3]
        elif len(args) > 0:
            raise NotImplementedError(
                """Screw6 Constructor does not accept %s positional
                arguments. See Usage.""" % (str(len(args)))
            )
        if len(kwargs) == 4:
            self.topleft = kwargs['tl']
            self.topright = kwargs['tr']
            self.botleft = kwargs['bl']
            self.botright = kwargs['br']
        elif len(kwargs) == 1:
            self.val = kwargs['value']
        elif len(kwargs) > 0:
            raise NotImplementedError(
                """Screw6 Constructor does not accept %s keyword
                arguments. See Usage.""" % (str(len(kwargs)))
            )

    def __str__(self):
        row_format = '[' + ((('{},' * 6) + ';') * 6) + ']'
        elements = list()
        for i in range(self._val.rows):
            for j in range(self._val.cols):
                elements.append(str(self._val[i, j]))
        str_format = row_format.format(*elements)
        return str_format

    def __repr__(self):
        repr_format = 'Screw6()'
        return repr_format

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

    def __eq__(self, other):
        """Check equality between two instances of Screw6."""
        if type(self) != type(other):
            raise ValueError(
                "Unable to compare %s with Screw6 type." % str(type(other))
            )
        return self.val == other.val

    def __ne__(self, other):
        """Check non-equality between two instances of Screw6."""
        return not self == other


