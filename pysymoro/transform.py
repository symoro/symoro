# -*- coding: utf-8 -*-


"""
This module contains the TransformationMatrix data structure.
"""


import sympy


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
    def __init__(self, frame_i, frame_j, params=None):
        """
        Constructor period.

        Args:
            frame_i: Frame with respect to which the transformation is
                computed (int). This is usually the antecedent frame.
            frame_j: Frame whose transformation is to be computed (int).
            params: A dict containing the geometric parameters.

        Usage:
        """
        self._frame_i = frame_i
        self._frame_j = frame_j
        self._gamma = 0
        self._b = 0
        self._alpha = 0
        self._d = 0
        self._theta = 0
        self._r = 0
        if params is not None:
            self.update(params)

    def update(self, params):
        """
        Update the parameter values and all the matrices accordingly.

        Args:
            params: A dict in which the keys correspond to the list of
                parameters that are to be updated and the values
                correspond to the values with which the parameters are
                to be updated.
        """
        for key, value in params.iteritems():
            attr = '_' + key
            if hasattr(self, attr):
                setattr(self, attr, value)
            else:
                raise AttributeError(
                    "%s is not a geometric parameter" % key
                )
        self._compute_tmat()
        self._compute_tinv()

    @property
    def val(self):
        """
        Get the value of the transformation matrix.

        Returns:
            A 4x4 Matrix
        """
        if not hasattr(self, '_tmat'):
            raise AttributeError("tmat is yet to be computed")
        return self._tmat

    @property
    def rot(self):
        """
        Get the value of the rotation matrix.

        Returns:
            A 3x3 Matrix
        """
        if not hasattr(self, '_tmat'):
            raise AttributeError("tmat is yet to be computed")
        return self._tmat[0:3, 0:3]

    @property
    def trans(self):
        """
        Get the value of the translation vector.

        Returns:
            A 3x1 Matrix
        """
        if not hasattr(self, '_tmat'):
            raise AttributeError("tmat is yet to be computed")
        return self._tmat[0:3, 3:4]

    @property
    def inv(self):
        """
        Get the inverse of the transformation matrix.

        Returns:
            A 4x4 Matrix
        """
        if not hasattr(self, '_tinv'):
            raise AttributeError("tinv is yet to be computed")
        return self._tinv

    @property
    def inv_rot(self):
        """
        Get the inverse of the rotation matrix.

        Returns:
            A 3x3 Matrix
        """
        if not hasattr(self, '_tinv'):
            raise AttributeError("tinv is yet to be computed")
        return self._tinv[0:3, 0:3]

    @property
    def inv_trans(self):
        """
        Get the inverse of the translation vector.

        Returns:
            A 3x1 Matrix
        """
        if not hasattr(self, '_tinv'):
            raise AttributeError("tinv is yet to be computed")
        return self._tinv

    def _compute_tmat(self):
        """
        Call (proxy method) the actual function that computes the
        transformation matrix between any two frames.
        """
        self._tmat = get_transformation_matrix(
            self._gamma, self._b,
            self._alpha, self._d,
            self._theta, self._r
        )

    def _compute_tinv(self):
        """
        Compute inverse of the transformation matrix.
        """
        if not hasattr(self, '_tmat'):
            raise AttributeError("tmat is yet to be computed")
        self._tinv = sympy.eye(4)
        rot_inv = self.rot.transpose()
        trans_inv = -rot_inv * self.trans
        self._tinv[0:3, 0:3] = rot_inv
        self._tinv[0:3, 3:4] = trans_inv


