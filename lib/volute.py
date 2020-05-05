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

    area(flow, c_thr):
        area at throat section
    area(flow, c_thr, theta):
        area at any section located by the angle theta

    :param flow (float): flow rate [m^3/s]
    :param c_thr (float): absolute velocity [m/s]
    :param theta (float): angle at which evaluate volute section [rad]
    :return a (float): area [m^2]
    """
    a = (flow / c_thr) * (theta / (2 * math.pi))

    return a


def angle_theta(n, i):
    """Calculate angle around volute axis.

    :param n (int): num. of divisions of the volute
    :param i (int): i-section
    :param b_3 (float): width at section 3
    :return theta (float): angle [rad]
    """
    theta = i * (2 * math.pi) / (n - 1)

    return theta


def diameter(d_2):
    """Calculate internal diameter at the start angle.

    :param d_2 (float): diameter of the impeller at section 2 [m]
    :return d (float): diameter [m]
    """
    d = round(1.1 * d_2, 3)

    return d


def width_min(b_2, a_thr):
    """Calculate minimum width and relative angle.

    :param b_2 (float): impeller width at section 2 [m]
    :param a_thr (float): area at throat section [m^2]
    :return b, theta (tuple): minimum volute width [m] and relative angle [rad]
    """
    b = round(1.5 * b_2, 3)
    theta = (math.pi * b)**2 / (2 * a_thr)

    return (b, theta)


def width(theta, a_thr, b_3):
    """Calculate width volute vane at different angles.

    :param theta (float): angle at which eval. volute section [rad]
    :param a_thr (float): area at throat section [m^2]
    :param b_3 (float): volute width at section 3 [m]
    :return b, theta (tuple): volute widths [m] and relative angles [rad]
    """
    b = (2 * a_thr * theta / math.pi**2)**.5

    if b > b_3:
        return (b, theta)
    else:
        return None
