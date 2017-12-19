# -*- coding: utf-8 -*-


"""
This module contains the TransformationMatrix data structure.
"""


import sympy

from pysymoro.screw6 import Screw6
from symoroutils import tools


def get_transformation_matrix(gamma, b, alpha, d, theta, r):
    """
    Compute and get transformation matrix between any two frames.

    Args:
        gamma, b, alpha, d, theta, r: The geometric parameter values
        that give the relationship between any two frames that are
        defined using the modified DH-method.
    Returns:
        A 4x4 Matrix that is the homogenous transformation matrix
        between the frames.
    """
    # simplify computation
    c_gamma = sympy.cos(gamma)
    s_gamma = sympy.sin(gamma)
    c_alpha = sympy.cos(alpha)
    s_alpha = sympy.sin(alpha)
    c_theta = sympy.cos(theta)
    s_theta = sympy.sin(theta)
    # intermediate terms
    sg_ca = s_gamma * c_alpha
    sg_sa = s_gamma * s_alpha
    cg_ca = c_gamma * c_alpha
    cg_sa = c_gamma * s_alpha
    # t matrix elements
    t11 = (c_gamma * c_theta) - (sg_ca * s_theta)
    t12 = -(c_gamma * s_theta) - (sg_ca * c_theta)
    t13 = sg_sa
    t14 = (d * c_gamma) + (r * sg_sa)
    t21 = (s_gamma * c_theta) + (cg_ca * s_theta)
    t22 = -(s_gamma * s_theta) + (cg_ca * c_theta)
    t23 = -cg_sa
    t24 = (d * s_gamma) - (r * cg_sa)
    t31 = s_alpha * s_theta
    t32 = s_alpha * c_theta
    t33 = c_alpha
    t34 = (r * c_alpha) + b
    # t matrix
    tmat = sympy.Matrix([
        [t11, t12, t13, t14],
        [t21, t22, t23, t24],
        [t31, t32, t33, t34],
        [0, 0, 0, 1]
    ])
    return tmat


class TransformationMatrix(object):
    """
    Data structure:
        Represent the data structure to hold a transformation matrix
        between any two frames.
    """
    def __init__(self, **kwargs):
        """Constructor

        Usage:
            >>> # Positional arguments are not allowed.
            >>> # Only frames are specified. i, j are of type int.
            TransformationMatrix(i=<frame-i>, j=<frame-j>)
            >>> # Frames and parameter values are specified. i, j are
            >>> # of type int.
            TransformationMatrix(
                i=<frame-i>, j=<frame-j>, params=<param-dict>
            )
            >>> # Only parameter values are specified. In this case the
            >>> # <param-dict> must contain the `frame` and `ant` keys
            >>> # with relevant values.
            TransformationMatrix(params=<param-dict>)
        """
        if len(kwargs) >= 1 and len(kwargs) <= 3:
            self._init_default()
            self._compute_tmat()
            self._compute_tinv()
            self._compute_smat()
            self._compute_sinv()
            if len(kwargs) > 1:
                self._frame_i = int(kwargs['i'])
                self._frame_j = int(kwargs['j'])
                if 'params' in kwargs:
                    self.update(kwargs['params'])
            else:
                if self._has_frames(kwargs['params']):
                    self.update(kwargs['params'])
                else:
                    self._cleanup()
                    raise AttributeError(
                        """Cannot setup TransformationMatrix without
                        specifying the two frames - frame, ant. Input
                        was: %s""" % str(kwargs['params'])
                    )
        else:
            raise NotImplementedError(
                """Wrong use of TransformationMatrix constructor. See
                Usage."""
            )

    def __str__(self):
        str_format = (
            "T matrix of frame %d wrt frame %d:\n"
            "----------------------------------\n"
            "gamma=%s, b=%s, alpha=%s, d=%s, theta=%s, r=%s\n"
            "%s\n"
            "**********************************\n"
        ) % (
            self._frame_j, self._frame_i,
            str(self._gamma), str(self._b),
            str(self._alpha), str(self._d),
            str(self._theta), str(self._r),
            sympy.pretty(self._tmat)
        )
        return str_format

    def __repr__(self):
        return 'T of {} wrt {}'.format(self._frame_j, self._frame_i)

    @property
    def val(self):
        """Get the value of the transformation matrix

        Returns:
            A 4x4 Matrix.
        """
        if not hasattr(self, '_tmat'):
            raise AttributeError('tmat is yet to be computed')
        return self._tmat

    @property
    def rot(self):
        """Get the value of the rotation matrix

        Returns:
            A 3x3 Matrix.
        """
        if not hasattr(self, '_tmat'):
            raise AttributeError('tmat is yet to be computed')
        return self._tmat[0:3, 0:3]

    @property
    def trans(self):
        """Get the value of the translation vector

        Returns:
            A 3x1 Matrix.
        """
        if not hasattr(self, '_tmat'):
            raise AttributeError('tmat is yet to be computed')
        return self._tmat[0:3, 3:4]

    @property
    def inv(self):
        """Get the inverse of the transformation matrix

        Returns:
            A 4x4 Matrix.
        """
        if not hasattr(self, '_tinv'):
            raise AttributeError('tinv is yet to be computed')
        return self._tinv

    @property
    def inv_rot(self):
        """Get the inverse of the rotation matrix

        Returns:
            A 3x3 Matrix.
        """
        if not hasattr(self, '_tinv'):
            raise AttributeError('tinv is yet to be computed')
        return self._tinv[0:3, 0:3]

    @property
    def inv_trans(self):
        """Get the inverse of the translation vector

        Returns:
            A 3x1 Matrix.
        """
        if not hasattr(self, '_tinv'):
            raise AttributeError('tinv is yet to be computed')
        return self._tinv[0:3, 3:4]

    @property
    def s_j_wrt_i(self):
        """Get the screw form of the transformation matrix

        Returns:
            A 6x6 Matrix.
        """
        return self._smat.val

    @property
    def s_i_wrt_j(self):
        """Get the screw form of the inverse transformation matrix

        Returns:
            A 6x6 Matrix.
        """
        return self._sinv.val

    def update(self, params):
        """Update the parameter values and all the matrices accordingly

        Args:
            params: A dict in which the keys correspond to the list of
                parameters that are to be updated and the values
                correspond to the values with which the parameters are
                to be updated.
        """
        for key, value in params.items():
            attr = '_' + key
            if hasattr(self, attr):
                setattr(self, attr, value)
            elif key in ['sigma', 'mu']:
                continue
            elif key is 'frame':
                self._frame_j = int(value)
            elif key is 'ant':
                self._frame_i = int(value)
            else:
                raise AttributeError(
                    '{} is not a geometric parameter'.format(key))
        self._compute_tmat()
        self._compute_tinv()
        self._compute_smat()
        self._compute_sinv()

    def _compute_tmat(self):
        """Compute the transformation matrix between any two frames

        Call (proxy method) the actual function that computes the
        transformation matrix between any two frames.
        """
        self._tmat = get_transformation_matrix(
            self._gamma, self._b,
            self._alpha, self._d,
            self._theta, self._r
        )

    def _compute_tinv(self):
        """Compute the inverse of the transformation matrix"""
        if not hasattr(self, '_tmat'):
            raise AttributeError('tmat is yet to be computed')
        self._tinv = sympy.eye(4)
        rot_inv = self.rot.transpose()
        trans_inv = -rot_inv * self.trans
        self._tinv[0:3, 0:3] = rot_inv
        self._tinv[0:3, 3:4] = trans_inv

    def _compute_smat(self):
        """Compute the transformation matrix in Screw (6x6) matrix form"""
        self._smat = Screw6(
            tl=self.rot, tr=tools.skew(self.trans),
            bl=sympy.zeros(3, 3), br=self.rot
        )

    def _compute_sinv(self):
        """Compute the inverse transformation matrix in screw matrix form

        Compute the inverse transformation matrix in screw (6x6) matrix form.
        """
        self._sinv = Screw6(
            tl=self.inv_rot, tr=-(self.inv_rot * tools.skew(self.trans)),
            bl=sympy.zeros(3, 3), br=self.inv_rot
        )

    def _init_default(self):
        """Initialise to 0 by default"""
        self._gamma = 0
        self._b = 0
        self._alpha = 0
        self._d = 0
        self._theta = 0
        self._r = 0

    def _has_frames(self, params):
        """Check if the frame and its antecedent are specified in a given dict

        Args:
            params: A dict containing the geometric parameters.

        Returns:
            True if `frame` and `ant` keys are present in params. False
            otherwise.
        """
        if 'frame' in params and 'ant' in params:
            return True
        else:
            return False

    def _cleanup(self):
        """Remove attributes of the data structure"""
        del self._tinv
        del self._tmat
        del self._gamma
        del self._b
        del self._alpha
        del self._d
        del self._theta
        del self._r
