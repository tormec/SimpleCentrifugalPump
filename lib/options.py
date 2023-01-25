"""Methods to calculate different impeller dimensioning options."""

import math
import lib.constants as CN


def type_number(omega, flow, head):
    """Calculate centrifugal pump's specific speed.

    :param omega (float): angular velocity [rad/s]
    :param flow (float): flow rate [m^3/s]
    :param head (float): head [m]
    :return cappa (float): specific speed
    """
    cappa = omega * flow**0.5 / (CN.G * head)**0.75

    return cappa


def rotational_speed(np, slip, hz):
    """Calculate rotational speed at different number of poles of an AC motor.

    :param np (int): number of poles of an AC motor
    :param slip (int): slip factor [%]
    :param hz (int): utility frequency [Hz]
    :return rpm (float): rotational speed [rpm]
    """
    rpm = 120 * hz / np * (1 - slip / 100)

    return rpm


def efficency_poly(cappa):
    """Calculate efficency for a given pump's specific speed.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta         .700 .850 .900 .916 .923 .928 .931 .932 .933 .935 .932
    weights     ones(cappa)
    n           5

    :param cappa (float): specific speed
    :return eta (float): efficency
    """
    coef = [-0.171, 7.400, -19.717, 25.671, -16.239, 3.990]
    eta = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return eta


def efficency_hyd_poly(cappa):
    """Calculate hydraulic efficency for a given pump's specific speed.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta_hyd     .790 .820 .851 .870 .900 .911 .915 .918 .919 .920 .918
    weights     ones(cappa)
    n           5

    :param cappa (float): specific speed
    :return eta_hyd (float): hydraulic efficency
    """
    coef = [0.778, -0.224, 2.011, -3.295, 2.162, -0.513]
    eta_hyd = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return eta_hyd


def efficency_vol_poly(cappa):
    """Calculate volumetric efficency for a given pump's specific speed.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta_hyd     .940 .948 .953 .956 .957 .958 .959 .960 .961 .962 .963
    weights     ones(cappa)
    n           5

    :param cappa (float): specific speed
    :return eta_vol (float): volumetric efficency
    """
    coef = [0.907, 0.236, -0.433, 0.378, -0.144, 0.016]
    eta_vol = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return eta_vol


def flow_number_poly(cappa):
    """Calculate flow number for a given pump's specific speed.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    phi         .080 .093 .100 .110 .120 .130 .140 .150 .160 .165 .170
    weights     ones(cappa)
    n           5

    :param cappa (float): specific speed
    :return phi (float): flow number
    """
    coef = [0.029, 0.416, -1.090, 1.665, -1.149, 0.288]
    phi = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return phi


def head_number_poly(cappa):
    """Calculate head number for a given pump's specific speed.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    psi         .583 .575 .560 .535 .515 .489 .465 .441 .415 .395 .380
    weights     ones(cappa)
    n           5

    :param cappa (float): specific speed
    :return psi (float): head number
    """
    coef = [0.531, 0.613, -2.339, 3.255, -2.284, 0.641]
    psi = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return psi


def phi2b(d, u, phi, flow, x=1, eta_vol=1):
    """Calculate impeller width for a given flow number.

    :param d (float): diameter [m]
    :param u (float): blade velocity [m/s]
    :param phi (float): flow number
    :param flow (float): flow rate [m^3/s]
    :param x (float): blade blockage
    :param eta_vol (float): volumetric efficency
    :return b (float): impeller width [m]
    """
    b = flow / (math.pi * d * u * phi * x * eta_vol)

    return b


def psi2u(psi, head):
    """Calculate blade velocity for a given head number.

    :param psi (float): head number
    :param head (float): head [m]
    :return u (float): blade velocity [m/s]
    """
    u = (CN.G * head / psi)**0.5

    return u


def cappa2npsh(cappa, head):
    """Calculate the npsh required for a given pump's specific speed.

    :param cappa (float): specific speed
    :param head (float): head [m]
    :return npsh_req (float): neat positive suction head required [m]
    """
    npsh_req = .25 * cappa**(4/3) * head

    return npsh_req


def width0diameter(b, d):
    """Calculate rate width over diameter.

    :param b (float): impeller width [m]
    :param d (float): diameter [m]
    :return bd (float): width over diameter
    """
    bd = b / d

    return bd
