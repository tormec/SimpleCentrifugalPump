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
    a = (flow / c_thr) * (theta / (2 * math.pi))

    return a


def angle_theta(n, i):
    """Calculate winding angle around volute axis.

    :param n (int): num. of divisions of the volute.
    :param i (int): section
    :param b_3 (float): width at section 3
    :return theta (float): winding angle [rad]
    """
    theta = i * (2 * math.pi) / (n - 1)

    return theta


def diameter(d_2):
    """Calculate internal diameter at the start winding angle.

    :param d_2 (float): diameter of the impeller at section 1 [m]
    :return r3 (float): radius [m]
    """
    d = 1.1 * d_2

    return d


def width_theta_min(b_2, a_thr):
    """Calculate minimum width and relative winding angle.

    :param b_2 (float): impeller width at section 2 [m]
    :param a_thr (float): area at throat section [m^2]
    """
    b = 1.8 * b_2
    theta = (math.pi * b)**2 / (2 * a_thr)

    return (b, theta)


def width(theta, a_thr, b_3):
    """Calculate width volute vane at different winding angles.

    :param theta (float): angle at which eval. volute section [rad]
    :param a_thr (float): area at throat section [m^2]
    :param b_3 (float): volute width at section 3 [m]
    :return b (list): diameters at different winding angles [m]
    """
    b = (2 * a_thr * theta / math.pi**2)**.5
    if b > b_3:
        return (b, theta)
    else:
        return None
