# -*- coding: utf-8 -*-


"""
This module contains the FloatingRobot data structure. Additionally as a
temporary measure some other data structures are included in this module
as well.
"""


from sympy import eye, var
from sympy import Matrix

from pysymoro.dynparams import DynParams
from pysymoro.geoparams import GeoParams
from pysymoro.screw import Screw
from symoroutils import filemgr
from symoroutils import tools


class FloatingRobot(object):
    """
    Represent the FloatingRobot data structure and provide methods that
    act as the gateway for robot modelling.
    """
    def __init__(
        self, name, links=0, joints=0, frames=0, 
        is_floating=True, structure=tools.TREE
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
        """Joint velocities."""
        self.qdots = [var('QD{0}'.format(j)) for j in self.joint_nums]
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
        self.base_acc = Screw()
        """Transformation matrix of base wrt a reference frame."""
        self.base_tmat = eye(4)

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
        str_format = str_format + ("\tStructure: %s\n" % str(self.structure))
        # add joint type - rigid or flexible
        str_format = str_format + "\tJoint Type: " + str(self.etas) + '\n'
        str_format = str_format + '\n'
        # add geometric params
        str_format = str_format + "Geometric Parameters:\n"
        str_format = str_format + "---------------------\n"
        for geo in self.geos:
            str_format = str_format + str(geo) + '\n'
        str_format = str_format + '\n'
        # add dynamic params
        str_format = str_format + "Dynamic Parameters:\n"
        str_format = str_format + "-------------------\n"
        for dyn in self.dyns:
            str_format = str_format + str(dyn) + '\n'
        str_format = str_format + '\n' + ('=*' * 60) + '='
        return str_format

    def __repr__(self):
        repr_format = (
            "(name=%s, links=%d, joints=%d, frames=%d, floating=%s, type=%s)"
        ) % (
            str(self.name), self.num_links, self.num_joints,
            self.num_frames, str(self.is_floating), str(self.structure)
        )
        return repr_format

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
        return xrange(self.num_links + 1)

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
        return xrange(self.num_joints + 1)

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
        return xrange(self.num_frames + 1)

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


