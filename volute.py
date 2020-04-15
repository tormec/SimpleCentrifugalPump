"""Methods to calculate the volute."""

import math


def absolute_velocity_throat(cu1):
    """Calculate absolute velocity at throat section.

    :param cu1 (float): peripheral velocity component [m/s]
    :return c_thr (float): absolute velocity [m/s]
    """
    c_thr = .5 * cu1

    return c_thr


def area(flow, c_thr, theta=2*math.pi):
    """Calculate area at i-section.

    :param flow (float): flow rate [m^3/s]
    :param c_thr (float): absolute velocity [m/s]
    :param theta (float): angle at which eval. volute section [rad]
    :return a (float): area [m^2]
    """
    a = (flow / c_thr) * (theta / 2 * math.pi)

    return a


def angle_theta(n, i):
    """Calculate winding angle around volute axis.

    :param n (int): num. of divisions of the volute.
    :param i (int): section
    :return theta (float): winding angle [rad]
    """
    theta_0 = math.radians(10)
    theta = theta_0 + i * (2 * math.pi - theta_0) / (n - 1)

    return theta


def diameter(d_2):
    """Calculate internal diameter at the start winding angle.

    :param d_2 (float): diameter of the impeller at section 1 [m]
    :return r3 (float): radius [m]
    """
    d = 1.13 * d_2

    return d


def width(theta, a_thr=None, b_2=None):
    """Calculate width volute vane at different winding angles.

    :param a_thr (float): area at throat section [m^2]
    :param theta (float): angle at which eval. volute section [rad]
    :param b_2 (float): impeller width at section 1 [m]
    :return b (list): diameters at different winding angles [m]
    """
    if b_2 is not None:
        b = 1.715 * b_2
    else:
        b = (2 * a_thr * theta / math.pi**2)**.5

    return b
