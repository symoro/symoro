# -*- coding: utf-8 -*-


"""
This module contains the FloatingRobot data structure. Additionally as a
temporary measure some other data structures are included in this module
as well.
"""


from sympy import eye, var
from sympy import Matrix

from pysymoro.screw import Screw
from pysymoro.dynparams import DynParams
from pysymoro.geoparams import GeoParams
from pysymoro import dynmodel
from symoroutils import filemgr
from symoroutils import tools


class FloatingRobot(object):
    """
    Represent the FloatingRobot data structure and provide methods that
    act as the gateway for robot modelling.
    """
    def __init__(
        self, name, links=0, joints=0, frames=0,
        is_floating=True, structure=tools.TREE, is_mobile=False,
        is_symbolic=True
    ):
        """
        Constructor period.

        Usage:
        """
        """Name of the robot."""
        self.name = name
        """Folder to store the files related to the robot."""
        self.directory = filemgr.get_folder_path(name)
        """Total number of links in the robot."""
        self.num_links = links
        """Total number of joints in the robot."""
        self.num_joints = joints
        """Total number of frames in the robot."""
        self.num_frames = frames
        """To indicate if the base is floating or not."""
        self.is_floating = is_floating
        """Type of the robot structure - simple, tree, closed-loop"""
        self.structure = structure
        """To indicate if the robot is a wheeled mobile robot"""
        self.is_mobile = is_mobile
        """To indicate if computation should be symbolic or numeric"""
        self.is_symbolic = is_symbolic
        # properties dependent on number of links
        """
        List to hold the dynamic parameters. The indices of the list
        start with 0 and it corresponds to parameters of link 0 (virtual
        link of the base).
        """
        self.dyns = [DynParams(j) for j in self.link_nums]
        # properties dependent on number of joints
        """
        To indicate if a joint is rigid or flexible. 0 for rigid and 1
        for flexible. The indices of the list start with 0 and
        corresponds to a virtual joint of the base. This joint is
        usually rigid.
        """
        self.etas = [0 for j in self.joint_nums]
        """Joint stiffness usually indicated by k."""
        self.stiffness = [0 for j in self.joint_nums]
        """Joint velocities."""
        self.qdots = [var('QP{0}'.format(j)) for j in self.joint_nums]
        """Joint accelerations."""
        self.qddots = [var('QDP{0}'.format(j)) for j in self.joint_nums]
        """Joint torques."""
        self.torques = [var('GAM{0}'.format(j)) for j in self.joint_nums]
        # properties dependent on number of frames
        """
        List to hold the geometric parameters. NOTE: This might be moved
        to a separate function.
        The indices of the list start with 0 and the first object
        corresponds to parameters of frame 0 (base) wrt its antecedent
        (some arbitary reference frame).
        """
        self.geos = [GeoParams(j) for j in self.frame_nums]
        # properties independent of number of links, joints and frames
        """Gravity vector a 3x1 Matrix."""
        self.gravity = Matrix([0, 0, var('G3')])
        # the values of properties below would be modified during
        # the computation of dynamic models.
        """Base velocity 6x1 column vector - a Screw."""
        self.base_vel = Screw()
        """Base acceleration 6x1 column vector - a Screw."""
        self.base_accel = Screw()
        """Transformation matrix of base wrt a reference frame at time 0."""
        self.base_tmat = eye(4)
        # call init methods
        self._init_maps()

    def __str__(self):
        str_format = ""
        # add robot description
        str_format = str_format + "Robot Description:\n"
        str_format = str_format + "------------------\n"
        str_format = str_format + ("\tName: %s\n" % str(self.name))
        str_format = str_format + ("\tLinks: %s\n" % str(self.num_links))
        str_format = str_format + ("\tJoints: %s\n" % str(self.num_joints))
        str_format = str_format + ("\tFrames: %s\n" % str(self.num_frames))
        str_format = str_format + ("\tFloating: %s\n" % str(self.is_floating))
        str_format = str_format + ("\tMobile: %s\n" % str(self.is_mobile))
        str_format = str_format + ("\tStructure: %s\n" % str(self.structure))
        str_format = str_format + '\n'
        # add geometric params
        str_format = str_format + "Geometric Parameters:\n"
        str_format = str_format + "---------------------\n"
        str_format = str_format + ('\t' + ('{:^8}' * 11) + '\n').format(*(
            'frame', 'ant', 'sigma', 'mu', 'gamma', 'b',
            'alpha', 'd', 'theta', 'r', 'q'
        ))
        for geo in self.geos:
            str_format = str_format + str(geo) + '\n'
        str_format = str_format + '\n'
        # add dynamic params
        str_format = str_format + "Dynamic Parameters:\n"
        str_format = str_format + "-------------------\n"
        str_format = str_format + ('\t' + ('{:^7}' * 20) + '\n').format(*(
            'link', 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ', 'MX', 'MY', 'MZ',
            'M', 'IA', 'Fc', 'Fv', 'FX', 'FY', 'FZ', 'CX', 'CY', 'CZ'
        ))
        for dyn in self.dyns:
            str_format = str_format + str(dyn) + '\n'
        str_format = str_format + '\n'
        # add joint params
        str_format = str_format + "Joint Parameters:\n"
        str_format = str_format + "-----------------\n"
        str_format = str_format + ('\t' + ('{:^9}' * 6) + '\n').format(*(
            'joint', 'eta', 'stiffness', 'qdot', 'qddot', 'torque'
        ))
        for jnt in self.joint_nums:
            jnt_str = ('\t' + ('{:^9}' * 6) + '\n').format(*(
                str(jnt), str(self.etas[jnt]), str(self.stiffness[jnt]),
                str(self.qdots[jnt]), str(self.qddots[jnt]),
                str(self.torques[jnt])
            ))
            str_format = str_format + jnt_str
        str_format = str_format + '\n'
        # end
        str_format = str_format + ('=*' * 60) + '='
        return str_format

    def __repr__(self):
        repr_format = (
            "(name=%s, links=%d, joints=%d, frames=%d, floating=%s, type=%s)"
        ) % (
            str(self.name), self.num_links, self.num_joints,
            self.num_frames, str(self.is_floating), str(self.structure)
        )
        return repr_format

    def get_val(self, idx, name):
        """
        Get the robot parameter values. The method is maninly to
        communicate with the UI.

        Args:
            idx: The joint/link/frame (index) value.
            name: The parameter name.

        Returns:
            The value corresponding to the name and index.
        """
        if name is 'Z':
            return 0
        elif name in self._dyn_params_map:
            attr = getattr(self, 'dyns')
            value = getattr(attr[idx], self._dyn_params_map[name])
            return value
        elif name in self._geo_params_map:
            attr = getattr(self, 'geos')
            value = getattr(attr[idx], self._geo_params_map[name])
            return value
        elif name in self._misc_params_map:
            attr = getattr(self, self._misc_params_map[name])
            value = attr[idx]
            return value
        elif name in self._base_params_map:
            return None

    def put_val(self, idx, name, value):
        """
        Modify the robot parameter values. This method is mainly to
        communicate with the UI. For other purposes, see update_params()
        method.

        Args:
            idx: The joint/link/frame (index) value.
            name: The parameter name.
            value: The value with which the parameter is to be modified.

        Returns:
            A `OK` if successful and `FAIL` otherwise.
        """
        if name is 'Z':
            return tools.OK
        elif name in self._dyn_params_map:
            param = {int(idx): {self._dyn_params_map[name]: value}}
            return self.update_params('dyns', param)
        elif name in self._geo_params_map:
            if name in ['frame', 'ant', 'sigma', 'mu']:
                value = int(value)
            param = {int(idx): {self._geo_params_map[name]: value}}
            return self.update_params('geos', param)
        elif name in self._misc_params_map:
            key = self._misc_params_map[name]
            param = {int(idx): {key: value}}
            return self.update_params('misc', param)
        elif name in self._base_params_map:
            return self.update_params('base', param)
        return tools.FAIL

    def update_params(self, kind, params):
        """
        Update the parameter values of the robot.

        Args:
            kind: A string metioning the paramter type. More specifically
                indicating whether the parameters correspond to geometric
                (`geos`), dynamic (`dyns`), base velocity (`base`),
                base acceleration (`base`) or others (`misc`).
            params: A nested dict containing the index values as the key
                and another dict containing the parameter names and
                values as the value. See Examples for a valid dict
                argument that is accepted by this method.

        Example 1:
            params = {
                2: {'ant': 0, 'sigma': 1, 'mu': 1},
                4: {'ant': 2, 'alpha': 0, 'd': 1, 'r': 0},
                7: {'ant': 5}
            }
        Example 2:
            params = {
                3: {'xx': 'XX3', 'msx': 'MX3', 'mass': 'M3'},
                5: {'yz': 0, 'xz': 0, 'xy': 0},
                2: {'fx_ext': 0, 'fy_ext': 0, 'mz_ext': 0},
                4: {'ia': 'IA4', 'frv': 0, 'frc': 0}
            }
        Example 3:
            params = {
                1: {'qdots': 'QP1', 'qddots': 'QDP1', 'torques': 'GAM1'},
                2: {'etas': 1, 'stiffness': 'k2'}
            }
        Example 4:
            params = {
                0: {'gravity': 'GX'},
                1: {'gravity': 'GY'},
                2: {'gravity': 'GZ'}
            }
        """
        if kind in ['dyns', 'geos']:
            self._update_dyns_geos(kind, params)
        elif kind is 'misc':
            self._update_misc(params)
        elif kind is 'base':
            raise NotImplementedError(
                "Yet to be supported."
            )
        else:
            errmsg = "`kind` can be ['dyns', 'geos', 'misc', 'base']. "
            errmsg = errmsg + ("Current input: {0}").format(kind)
            raise ValueError(errmsg)
        return tools.OK

    def set_dyns_to_zero(self, links=None):
        """
        Set all the dynamic parameter values to zero for a specified
        list of links. If no link is specified, dynamic parameters for
        all links are set to zero.

        Args:
            links: An iterable object with the link numbers.
        """
        if links == None:
            links = self.link_nums
        for link in links:
            if link in self.link_nums:
                self.dyns[link].set_to_zero()
            else:
                err_msg = "Link number {} does not belong to the robot."
                raise IndexError(err_msg.format(link))

    def compute_idym(self):
        """
        Compute the Inverse Dynamic Model of the robot using the
        recursive Newton-Euler algorithm.
        """
        self.idym = dynmodel.inverse_dynamic_model(self)

    def compute_ddym(self):
        """
        Compute the Direct Dynamic Model of the robot using the
        recursive Newton-Euler algorithm.
        """
        self.ddym = dynmodel.direct_dynamic_model(self)

    @property
    def link_nums(self):
        """
        Get the link numbers.

        Returns:
            An iteratable object with the link numbers.
        Note:
            Add 1 to number of links. This serves two purposes -
            one, it indicates a virtual link - to represent the base and
            two, it makes sure list index 1 corresponds to link 1 and so
            on.
        """
        return range(self.num_links + 1)

    @property
    def joint_nums(self):
        """
        Get the joint numbers.

        Returns:
            An iteratable object with the joint numbers.
        Note:
            Add 1 to number of joints. This serves two purposes -
            one, it indicates a virtual joint 0 to represent the base and
            two, it makes sure list index 1 corresponds to joint 1, list
            index 2 to joint 2 and so on.
        """
        return range(self.num_joints + 1)

    @property
    def frame_nums(self):
        """
        Get the frame numbers.

        Returns:
            An iteratable object with the frame numbers.
        Note:
            Add 1 to number of frames in order to include the base
            frame as well. This also makes the list index 1 correspond
            to frame 1 and so on.
        """
        return range(self.num_frames + 1)

    @property
    def q_vec(self):
        """
        Get the list of joint variables.

        Returns:
            A list containing the joint variables.
        """
        q = list()
        for j in self.joint_nums:
            if j == 0:
                continue
            q.append(self.geos[j].q)
        return q

    @property
    def q_passive(self):
        """
        Get the list of passive joint variables.

        Returns:
            A list containing the passive joint variables.
        """
        q = list()
        for j in self.joint_nums:
            if j == 0: continue
            if self.geos[j].mu == 0 and self.geos[j].sigma != 2:
                q.append(self.geos[j].q)
        return q

    @property
    def q_active(self):
        """
        Get the list of active joint variables.

        Returns:
            A list containing the active joint variables.
        """
        q = list()
        for j in self.joint_nums:
            if j == 0: continue
            if self.geos[j].mu == 1 and self.geos[j].sigma != 2:
                q.append(self.geos[j].q)
        return q

    @property
    def passive_joints(self):
        """
        Get the list of joint numbers (indices) corresponding to passive
        joints.

        Returns:
            A list containing the passive joint numbers (indices).
        """
        joints = list()
        for j in self.joint_nums:
            if j == 0: continue
            if self.geos[j].mu == 0 and self.geos[j].sigma != 2:
                joints.append(j)
        return joints

    @property
    def active_joints(self):
        """
        Get the list of joint numbers (indices) corresponding to active
        joints.

        Returns:
            A list containing the active joint numbers (indices).
        """
        joints = list()
        for j in self.joint_nums:
            if j == 0: continue
            if self.geos[j].mu == 1 and self.geos[j].sigma != 2:
                joints.append(j)
        return joints

    def _update_base(self, params):
        """
        Update the base velocity and acceleration values of the robot.
        """
        pass

    def _update_misc(self, params):
        """
        Update the miscellaneous parameter values of the robot.
        """
        for key in params:
            idx = int(key)
            curr_params = params[key]
            for attr in curr_params:
                if not hasattr(self, attr):
                    raise AttributeError((
                        "{0} is not an attribute of Robot."
                    ).format(str(attr)))
                else:
                    try:
                        param = getattr(self, attr)
                        param[idx] = curr_params[attr]
                    except IndexError:
                        raise IndexError((
                            "`{0}` doesnt have {1} index value."
                        ).format(attr, str(idx)))

    def _update_dyns_geos(self, attr, params):
        """
        Update the geometric and dynamic parameter values of the robot.
        """
        param = getattr(self, attr)
        for key in params:
            idx = int(key)
            try:
                param[idx].update_params(params[key])
            except IndexError:
                raise IndexError((
                    "`{0}` doesnt have {1} index value."
                ).format(attr, str(idx)))

    def _init_maps(self):
        """
        Initialise maps (dict) that will be used to read and write the
        different parameter values of the robot. The main purpose of
        these maps is to talk easily with the user interface.
        """
        self._dyn_params_map = dict([
            ('XX', 'xx'), ('XY', 'xy'), ('XZ', 'xz'),
            ('YY', 'yy'), ('YZ', 'yz'), ('ZZ', 'zz'),
            ('MX', 'msx'), ('MY', 'msy'), ('MZ', 'msz'), ('M', 'mass'),
            ('IA', 'ia'), ('FS', 'frc'), ('FV', 'frv'),
            ('FX', 'fx_ext'), ('FY', 'fy_ext'), ('FZ', 'fz_ext'),
            ('CX', 'mx_ext'), ('CY', 'my_ext'), ('CZ', 'mz_ext')
        ])
        self._geo_params_map = dict([
            ('j', 'frame'), ('ant', 'ant'), ('sigma', 'sigma'),
            ('mu', 'mu'), ('gamma', 'gamma'), ('b', 'b'),
            ('alpha', 'alpha'), ('d', 'd'), ('theta', 'theta'), ('r', 'r')
        ])
        self._base_params_map = dict([
            ('V0', 'base_vel'), ('VP0', 'base_accel'),
            ('W0', 'base_vel'), ('WP0', 'base_accel')
        ])
        self._misc_params_map = dict([
            ('QP', 'qdots'), ('QDP', 'qddots'), ('GAM', 'torques'),
            ('eta', 'etas'), ('k', 'stiffness'), ('G', 'gravity'),
        ])

    # methods for backward-compatability with the old Robot (for fixed
    # base) class.


