# -*- coding: utf-8 -*-


"""
This module contains the DynParams data structure.
"""


import re

from sympy import eye, var
from sympy import Matrix

from pysymoro.screw import Screw
from pysymoro.screw6 import Screw6
from symoroutils import tools


class DynParams(object):
    """
    Data structure:
        Represent the data structure to hold the inertial parameters,
        friction parameters and the external forces for a given link. An
        instance of the inertia matrix, the spatial inertia matrix and
        the mass tensor term are also maintained.
    Note:
    Mass tensor refers to the MS 3x1 matrix which is the first moments
    of a link wrt its own frame of reference. MS = transpose([MX MY MZ])
    """
    def __init__(self, link, params=None):
        """Constructor

        Usage:
        DynParams(link=<link-number>)
        DynParams(link=<link-number>, params=<params-dict>)
        """
        self.link = link
        # Inertia matrix terms"""
        self.xx = None
        self.xy = None
        self.xz = None
        self.yy = None
        self.yz = None
        self.zz = None
        # Mass tensor terms"""
        self.msx = None
        self.msy = None
        self.msz = None
        # Link mass"""
        self.mass = None
        # Rotor inertia term"""
        self.ia = None
        # Coulomb friction parameter"""
        self.frc = None
        # Viscous friction parameter"""
        self.frv = None
        # External forces and moments"""
        self.fx_ext = None
        self.fy_ext = None
        self.fz_ext = None
        self.mx_ext = None
        self.my_ext = None
        self.mz_ext = None
        # lists to hold the string representation for the prefix of
        # different terms
        self._inertial_terms = {
            'xx': 'XX',
            'xy': 'XY',
            'xz': 'XZ',
            'yy': 'YY',
            'yz': 'YZ',
            'zz': 'ZZ'
        }
        self._ms_terms = {
            'msx': 'MX',
            'msy': 'MY',
            'msz': 'MZ',
            'mass': 'M'
        }
        self._fr_terms = {
            'ia': 'IA',
            'frc': 'FS',
            'frv': 'FV'
        }
        self._ext_force_terms = {
            'fx_ext': 'FX',
            'fy_ext': 'FY',
            'fz_ext': 'FZ',
            'mx_ext': 'CX',
            'my_ext': 'CY',
            'mz_ext': 'CZ'
        }
        # initialise the different parameters
        self._init_inertial_terms()
        self._init_ms_terms()
        self._init_fr_terms()
        self._init_ext_force_terms()
        # initialise with values if available
        if params is not None:
            self.update_params(params)

    def __str__(self):
        row_format = '\t' + ('{:^7}' * 20)
        str_format = row_format.format(*(
            str(self.link),
            str(self.xx), str(self.xy), str(self.xz),
            str(self.yy), str(self.yz), str(self.zz),
            str(self.msx), str(self.msy), str(self.msz), str(self.mass),
            str(self.ia), str(self.frc), str(self.frv),
            str(self.fx_ext), str(self.fy_ext), str(self.fz_ext),
            str(self.mx_ext), str(self.my_ext), str(self.mz_ext)
        ))
        return str_format

    def __repr__(self):
        repr_format = str(self).lstrip().rstrip()
        repr_format = re.sub(r'\s+', ', ', repr_format)
        repr_format = '(' + repr_format + ')'
        return repr_format

    def update_params(self, params):
        """Update the dynamic parameter values

        Args:
            params: A dict in which the keys correspond to the list of
                parameters that are to be updated and the values
                correspond to the values with which the parameters are
                to be updated.
        """
        for key, value in params.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(
                    '{} is not an attribute of DynParams'.format(key))

    def set_to_zero(self):
        """
        Set all the dynamic parameter values to zero.
        """
        # attributes that shouldn't be set to zero
        outliers = [
            'link', 'update_params', 'inertia', 'mass_tensor',
            'spatial_inertia', 'wrench', 'force', 'moment'
        ]
        # get all the attributes of the class ignoring system attributes
        attrs = [
            attr for attr in dir(self) \
            if not attr.startswith('_')
        ]
        for attr in attrs:
            # ignore attributes in outliers
            if attr in outliers:
                continue
            setattr(self, attr, 0)

    @property
    def inertia(self):
        """Get inertia (3x3) matrix."""
        return Matrix([
            [self.xx, self.xy, self.xz],
            [self.xy, self.yy, self.yz],
            [self.xz, self.yz, self.zz]
        ])

    @property
    def mass_tensor(self):
        """Get mass tensor (3x1) column vector."""
        return Matrix([self.msx, self.msy, self.msz])

    @property
    def spatial_inertia(self):
        """Get spatial inertia (6x6) matrix."""
        m_eye = self.mass * eye(3)
        ms_skew = tools.skew(self.mass_tensor)
        return Screw6(
            tl=m_eye, tr=ms_skew.transpose(),
            bl=ms_skew, br=self.inertia
        )

    @property
    def wrench(self):
        """Get external force (6x1) column vector (linear + angular)."""
        return Screw(
            lin=Matrix([self.fx_ext, self.fy_ext, self.fz_ext]),
            ang=Matrix([self.mx_ext, self.my_ext, self.mz_ext])
        )

    @property
    def force(self):
        """Get external force (3x1) column vector (linear)."""
        return self.wrench.lin

    @property
    def moment(self):
        """Get external moment (3x1) column vector (angular)."""
        return self.wrench.ang

    def _init_inertial_terms(self):
        """Initialise inertial terms"""
        for key, term in self._inertial_terms.items():
            if self.link != 0:
                value = var(term + str(self.link))
            else:
                value = 0
            setattr(self, key, value)

    def _init_ms_terms(self):
        """Initialise mass tensor terms and mass of the link"""
        for key, term in self._ms_terms.items():
            if self.link != 0:
                value = var(term + str(self.link))
            else:
                value = 0
            setattr(self, key, value)

    def _init_fr_terms(self):
        """Initialise rotor inertia and friction parameters."""
        for key, term in self._fr_terms.items():
            if self.link != 0:
                value = var(term + str(self.link))
            else:
                value = 0
            setattr(self, key, value)

    def _init_ext_force_terms(self):
        """Initialise external force terms."""
        for key, term in self._ext_force_terms.items():
            if self.link != 0:
                value = var(term + str(self.link))
            else:
                value = 0
            setattr(self, key, value)


