#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import math
import random
import unittest

from sympy import var

from symoroutils import tools
from pysymoro.robotf import FloatingRobot as Robot


def planar2r():
    robo = Robot('planar2r', 2, 2, 3, False, tools.SIMPLE)
    # update geometric params
    params = {
        1: {'sigma': 0, 'mu': 1, 'theta': var('th1')},
        2: {'sigma': 0, 'mu': 1, 'alpha': 0, 'd': var('L1'), 'theta': var('th2')},
        3: {'sigma': 2, 'd': var('L2')}
    }
    robo.update_params('geos', params)
    # update dynamic params
    params = {
        1: {'xx': 0, 'xy': 0, 'xz': 0, 'yy': 0, 'yz': 0, 'ia': 0}
    }
    robo.update_params('dyns', params)
    # update joint params
    params = {
        0: {'qdots': 0, 'qddots': 0, 'torques': 0}
    }
    robo.update_params('misc', params)
    # update gravity vector
    params = {
        0: {'gravity': 0},
        1: {'gravity': 0},
        2: {'gravity': var('G3')},
    }
    robo.update_params('misc', params)
    return robo


def planar2r_numerical(is_floating=False):
    # numerical values (random in some case)
    L1 = 0.5
    L2 = 0.4
    ZZ1 = 3.7
    ZZ2 = 0.35
    MX2 = 0.4
    MY2 = 0.15
    M1 = 1.2
    M2 = 0.8
    FC1 = 0.3
    FC2 = 0.25
    FV1 = 0.3
    FV2 = 0.18
    IA1 = 0.0
    IA2 = 0.0
    G3 = -9.81
    # create robot
    robo = Robot(
        'planar2r', 2, 2, 3, is_floating,
        tools.SIMPLE, is_symbolic=False
    )
    robo.set_dyns_to_zero()
    # update geometric params
    params = {
        1: {'sigma': 0, 'mu': 1},
        2: {'sigma': 0, 'mu': 1, 'd': L1},
        3: {'sigma': 2, 'd': L2}
    }
    robo.update_params('geos', params)
    # update dynamic params
    params = {
        1: {
            'zz': ZZ1, 'frc': FC1, 'frv': FV1, 'ia': IA1, 'mass': M1
        },
        2: {
            'zz': ZZ2, 'frc': FC2, 'frv': FV2, 'ia': IA2, 'mass': M2,
            'msx': MX2, 'msy': MY2
        }
    }
    robo.update_params('dyns', params)
    # update joint params
    params = {0: {'qdots': 0, 'qddots': 0, 'torques': 0}}
    robo.update_params('misc', params)
    # update gravity vector
    params = {
        0: {'gravity': 0},
        1: {'gravity': 0},
        2: {'gravity': G3},
    }
    robo.update_params('misc', params)
    return robo


def set_planar2r_joint_state(robo, q, qdot, qddot):
    """
    Set the joint states for a 2R planar robot.

    Args:
        robo: An instance of Robot class.
        q: Joint positions. A list of size 2.
        qdot: Joint velocities. A list of size 2.
        qddot: Joint accelerations. A list of size 2.

    Returns:
        The modified instance of Robot class.
    """
    q1 = q[0]
    q2 = q[1]
    q1dot = qdot[0]
    q2dot = qdot[1]
    q1ddot = qddot[0]
    q2ddot = qddot[1]
    # update joint variables
    params = {1: {'theta': q1}, 2: {'theta': q2}}
    robo.update_params('geos', params)
    # update joint params
    params = {
        1: {'qdots': q1dot, 'qddots': q1ddot},
        2: {'qdots': q2dot, 'qddots': q2ddot}
    }
    robo.update_params('misc', params)
    return robo


def set_planar2r_joint_torque(robo, qtorque):
    """
    Set the joint torques for a 2R planar robot.

    Args:
        robo: An instance of Robot class.
        qtorque: Joint torques. A list of size 2.

    Returns:
        The modified instance of Robot class.
    """
    gam1 = qtorque[0]
    gam2 = qtorque[1]
    # update torque values
    params = {1: {'torques': gam1}, 2: {'torques': gam2}}
    robo.update_params('misc', params)
    return robo


class TestDynModelPlanar2rFixed(unittest.TestCase):
    """
    Unit test for testing the inverse and direct dynamic model
    computation for floating base robots algorithm. This testing is done
    numerically.
    """
    def setUp(self):
        pass

    def test_when_zero(self):
        """
        Test the dynamic model computation when joint position,
        velocity and acceleration are set to zero.
        """
        robo = planar2r_numerical()
        # set joint state
        robo = set_planar2r_joint_state(
            robo, [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]
        )
        # compute IDyM
        robo.compute_idym()
        # set torque values for DDyM
        robo = set_planar2r_joint_torque(
            robo, [robo.idym.torques[1], robo.idym.torques[2]]
        )
        # compute DDyM
        robo.compute_ddym()
        # check if the result of IDyM (computed torques) are zero
        self.assertEqual(robo.idym.torques[1], 0.0)
        self.assertEqual(robo.idym.torques[2], 0.0)
        # check if the result of DDyM (computed qddots) are zero
        self.assertEqual(robo.ddym.qddots[1], 0.0)
        self.assertEqual(robo.ddym.qddots[2], 0.0)
        # check if input to IDyM is same as output of DDyM
        self.assertEqual(robo.ddym.qddots[1], robo.qddots[1])
        self.assertEqual(robo.ddym.qddots[2], robo.qddots[2])

    def test_when_random(self):
        """
        Test the dynamic model computation when joint position,
        velocity and acceleration are set to random meaningful values.
        """
        robo = planar2r_numerical()
        # initialise joint position, velocity, acceleration to random
        # values
        random.seed(math.pi)
        q = list(random.uniform(-math.pi, math.pi) for j in range(2))
        qdot = list(random.uniform(-math.pi, math.pi) for j in range(2))
        qddot = list(random.uniform(-math.pi, math.pi) for j in range(2))
        # set joint state
        robo = set_planar2r_joint_state(robo, q, qdot, qddot)
        # compute IDyM
        robo.compute_idym()
        print(robo.idym)
        # set torque values for DDyM
        robo = set_planar2r_joint_torque(
            robo, [robo.idym.torques[1], robo.idym.torques[2]]
        )
        # compute DDyM
        robo.compute_ddym()
        #
        print('\n')
        print(q)
        print(qdot)
        print(qddot)
        print(robo.idym.torques)
        print(robo.ddym.qddots)
        # check if input to IDyM is same as output of DDyM
        self.assertEqual(robo.ddym.qddots[1], robo.qddots[1])
        self.assertEqual(robo.ddym.qddots[2], robo.qddots[2])


class TestDynModelPlanar2rFloating(unittest.TestCase):
    """
    Unit test for testing the inverse and direct dynamic model
    computation for floating base robots algorithm. This testing is done
    numerically.
    """
    def setUp(self):
        pass

    def test_when_zero(self):
        """
        Test the dynamic model computation when joint position,
        velocity and acceleration are set to zero.
        """
        robo = planar2r_numerical(is_floating=True)
        # set joint state
        robo = set_planar2r_joint_state(
            robo, [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]
        )
        # compute IDyM
        robo.compute_idym()
        # set torque values for DDyM
        robo = set_planar2r_joint_torque(
            robo, [robo.idym.torques[1], robo.idym.torques[2]]
        )
        # compute DDyM
        robo.compute_ddym()
        # assertions
        # check if the result of IDyM (computed torques) are zero
        self.assertEqual(robo.idym.torques[1], 0.0)
        self.assertEqual(robo.idym.torques[2], 0.0)
        # check if the result of DDyM (computed qddots) are zero
        self.assertEqual(robo.ddym.qddots[1], 0.0)
        self.assertEqual(robo.ddym.qddots[2], 0.0)
        # check if input to IDyM is same as output of DDyM
        self.assertEqual(robo.ddym.qddots[1], robo.qddots[1])
        self.assertEqual(robo.ddym.qddots[2], robo.qddots[2])

    def test_when_random(self):
        """
        Test the dynamic model computation when joint position,
        velocity and acceleration are set to random meaningful values.
        """
        robo = planar2r_numerical(is_floating=True)
        # initialise joint position, velocity, acceleration to random
        # values
        random.seed(math.pi)
        q = list(random.uniform(-math.pi, math.pi) for j in range(2))
        qdot = list(random.uniform(-math.pi, math.pi) for j in range(2))
        qddot = list(random.uniform(-math.pi, math.pi) for j in range(2))
        # set joint state
        robo = set_planar2r_joint_state(robo, q, qdot, qddot)
        # compute IDyM
        robo.compute_idym()
        # set torque values for DDyM
        robo = set_planar2r_joint_torque(
            robo, [robo.idym.torques[1], robo.idym.torques[2]]
        )
        # compute DDyM
        robo.compute_ddym()
        #
        print('\n')
        print(q)
        print(qdot)
        print(qddot)
        print(robo.idym.torques)
        print(robo.ddym.qddots)
        # assertions
        # check if input to IDyM is same as output of DDyM
        self.assertEqual(robo.ddym.qddots[1], robo.qddots[1])
        self.assertEqual(robo.ddym.qddots[2], robo.qddots[2])


def run_tests():
    """Load and run the unittests"""
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(
        TestDynModelPlanar2rFixed
    )
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()


