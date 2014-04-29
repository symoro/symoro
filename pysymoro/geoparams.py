# -*- coding: utf-8 -*-


"""
This module contains the GeoParams data structure.
"""


class GeoParams(object):
    """
    Data structure:
        Represent the data structure to hold the geometric parameters.
    """
    def __init__(self, frame, params=None):
        """
        Constructor period.

        Usage:
        GeoParams(frame=<frame-number>)
        GeoParams(frame=<frame-number>, params=<params-dict>)
        """
        self.frame = frame
        self.ant = 0
        self.sigma = 0
        self.mu = 0
        self.gamma = 0
        self.b = 0
        self.alpha = 0
        self.d = 0
        self.theta = 0
        self.r = 0
        # initialise with values if available
        if params is not None:
            self.update_params(params)

    def update_params(self, params):
        """
        Update the geometric parameter values.

        Args:
            params: A dict in which the keys correspond to the list of
                parameters that are to be updated and the values
                correspond to the values with which the parameters are
                to be updated.
        """
        for key, value in params.iteritems():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(
                    "%s is not an attribute of GeoParams" % key
                )

    @property
    def q(self):
        """
        Get the joint variable value.

        Returns:
            The value of theta for a revolute joint and in the case of
            prismatic joint the value of r. For a fixed joint 0 is
            returned.
        """
        if self.sigma == 2:
            return 0
        return ((1 - self.sigma) * self.theta) + (self.sigma * self.r)


