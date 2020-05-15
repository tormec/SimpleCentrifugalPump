"""Methods to calculate the impeller."""

import math
import lib.constants as CN


def type_number(omega, flow, head):
    """Calculate centrifugal pump's typical number.

    :param omega (float): angular velocity [rad/s]
    :param flow (float): flow rate [m^3/s]
    :param head (float): head [m]
    :return cappa (float): typical number
    """
    cappa = omega * flow**0.5 / (CN.G * head)**0.75

    return cappa


def rpm2omega(rpm):
    """Calculate angular velocity for a given rotationl speed.

    :param rpm (float): rotational speed [rpm]
    :return omega (float): angular velocity [rad/s]
    """
    omega = 2 * math.pi * rpm / 60

    return omega


def efficency_poly(cappa):
    """Calculate efficency for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta         .700 .850 .900 .916 .923 .928 .931 .932 .933 .935 .932
    weights     ones(cappa)
    n           5

    :param cappa (float): typical number
    :return eta (float): efficency
    """
    coef = [-0.171, 7.400, -19.717, 25.671, -16.239, 3.990]
    eta = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return eta


def efficency_hyd_poly(cappa):
    """Calculate hydraulic efficency for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta_hyd     .790 .820 .851 .870 .900 .911 .915 .918 .919 .920 .918
    weights     ones(cappa)
    n           5

    :param cappa (float): typical number
    :return eta_hyd (float): hydraulic efficency
    """
    coef = [0.778, -0.224, 2.011, -3.295, 2.162, -0.513]
    eta_hyd = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return eta_hyd


def efficency_vol_poly(cappa):
    """Calculate volumetric efficency for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta_hyd     .940 .948 .953 .956 .957 .958 .959 .960 .961 .962 .963
    weights     ones(cappa)
    n           5

    :param cappa (float): typical number
    :return eta_vol (float): volumetric efficency
    """
    coef = [0.907, 0.236, -0.433, 0.378, -0.144, 0.016]
    eta_vol = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return eta_vol


def flow_number_poly(cappa):
    """Calculate flow number for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    phi         .080 .093 .100 .110 .120 .130 .140 .150 .160 .165 .170
    weights     ones(cappa)
    n           5

    :param cappa (float): typical number
    :return phi (float): flow number
    """
    coef = [0.029, 0.416, -1.090, 1.665, -1.149, 0.288]
    phi = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return phi


def flow_number(d, b, u, flow):
    """Calculate flow number.

    :param d (float): diameter [m]
    :param b (float): impeller width [m]
    :param u (float): absolute velocity [m/s]
    :param flow (float): flow rate [m^3/s]
    :return phi (float): flow number
    """
    phi = flow / (math.pi * d * b * u)

    return phi


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


def head_number_poly(cappa):
    """Calculate head number for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    psi         .583 .575 .560 .535 .515 .489 .465 .441 .415 .395 .380
    weights     ones(cappa)
    n           5

    :param cappa (float): typical number
    :return psi (float): head number
    """
    coef = [0.531, 0.613, -2.339, 3.255, -2.284, 0.641]
    psi = sum([val * cappa**idx for idx, val in enumerate(coef)])

    return psi


def head_number(u, head):
    """Calculate head number.

    :param u (float): absolute velocity [m/s]
    :param head (float): head [m]
    :return psi (float): head number
    """
    psi = (CN.G * head) / u**2

    return psi


def theoretic_head_number(psi, eta_hyd):
    """Calculate theoretic head number.

    :param psi (float): head number
    :param eta_hyd (float): hydraulic efficency
    :return psi_th (float): theoretic head number
    """
    psi_th = psi / eta_hyd

    return psi_th


def psi2u(psi, head):
    """Calculate blade velocity for a given head number.

    :param psi (float): head number
    :param head (float): head [m]
    :return u (float): blade velocity [m/s]
    """
    u = (CN.G * head / psi)**0.5

    return u


def width0diameter(b, d):
    """Calculate rate width over diameter.

    :param b (float): impeller width [m]
    :param d (float): diameter [m]
    :return bd (float): width over diameter
    """
    bd = b / d

    return bd


def npsh_req(c, w, lm, lw):
    """Calculate the neat positive suction head required.

    :param c (float): absolute velocity [m/s]
    :param w (float): relative velocity [m/s]
    :param lm (float): loss coefficient
    :param lw (float): low-pressure peak coefficient at blades
    :return npsh_req (float): neat positive suction head required [m]
    """
    npsh_req = lw * w**2 / (2 * CN.G) + (1 + lm) * c**2 / (2 * CN.G)

    return npsh_req


def cappa2npsh(cappa, head):
    """Calculate the npsh required for a given pump's typical number.

    :param cappa (float): typical number
    :param head (float): head [m]
    :return npsh_req (float): neat positive suction head required [m]
    """
    npsh_req = .25 * cappa**(4/3) * head

    return npsh_req


def hub_blockage(d, d_hu):
    """Calculate hub blockage.

    :param d (float): diameter [m]
    :param d_hu (float): hub diameter [m]
    :return x (float): hub blockage factor
    """
    x = 1 - (d_hu / d)**2

    return x


def blade_blockage(beta_b, d, t, z):
    """Calculate blade blockage.

    :param beta_b (float): angle between relative and blade velocity [m/s]
    :param d (float): diameter [m]
    :param t (float): blade thickness [m]
    :param z (int): number of blades
    :return x (float): blade blockage factor
    """
    x = 1 - (z * t) / (math.pi * d * math.sin(beta_b))

    return x


def diameter_omega(omega, u):
    """Calculate diameter function of angular velocity.

    :param omega (float): angular velocity [rad/s]
    :param u (float): absolute velocity [m/s]
    :return d (float): diameter [m]
    """
    d_omega = 2 * u / omega

    return d_omega


def diameter_npsh(omega, x, flow, lm, lw, km, eta_vol):
    """Calculate diameter favouring min npsh required.

    :param omega (float): angular velocity [rad/s]
    :param x (float): blockage factor
    :param flow (float): flow rate [m^3/s]
    :param lm (float): loss coefficient
    :param lw (float): low-pressure peak coefficient at blades
    :param km (float): rate between circumeferential velocity cm2 and c0
    :param eta_vol (float): volumetric efficency
    :return d_npsh (float): diameter with min npsh required [m]
    """
    d_npsh = 2 * ((2 * flow**2 * km**2 * (1 + lm + lw)) /
                  (eta_vol**2 * math.pi**2 * omega**2 * x**2 * lw))**(1/6)

    return d_npsh


def diameter_efficency(omega, x, flow, km, eta_vol):
    """Calculate diameter favouring max total efficency.

    :param omega (float): angular velocity [rad/s]
    :param x (float): hub blockage
    :param flow (float): flow rate [m^3/s]
    :param km (float): rate between c_1m and c_0 velocity
    :param eta_vol (float): volumetric efficency
    :return d_eff (float): diameter with max efficency [m]
    """
    d_eff = 2 * ((2 * flow**2 * km**2) /
                 (eta_vol**2 * math.pi**2 * omega**2 * x**2))**(1/6)

    return d_eff


def diameter_flow(omega, x, flow, eta_vol):
    """Calculate diameter for a given flow rate.

    :param omega (float): angular velocity [rad/s]
    :param x (float): hub blockage
    :param flow (float): flow rate [m^3/s]
    :param eta_vol (float): volumetric efficency
    :param eta_vol (float): volumetric efficency
    :return d_flow (float): diameter as function of the flow rate [m]
    """
    d_flow = ((flow * 8 * 3.03) / (omega * math.pi * x * eta_vol))**(1/3)

    return d_flow


def average_diam(d_npsh, d_eff, d_flow):
    """Calculate diameter as average value.

    :param d_npsh (float): diameter with min NPSH_r [m]
    :param d_eff (float): diameter with max efficency [m]
    :param d_flow (float): diameter as function of the flow rate [m]
    :return d_avg (float): diameter as average value [m]
    """
    d_types = [d_npsh, d_eff, d_flow]
    d_avg = math.fsum(d_types) / len(d_types)

    return d_avg


def standard_diam(d):
    """Calculate diameter according to standard diameters.

    :param d (float): diameter [m]
    :return d_std (float): standard diameter [m]
    """
    dif_abs = []
    for i in range(len(CN.D_INT)):
        dif_abs.append(abs(CN.D_INT[i] - d))
    dif_min = min(dif_abs)
    d_std = CN.D_INT[dif_abs.index(dif_min)]

    return d_std


def curvature_rad(d_2):
    """Calculate curvature radius of the shroud at section 0.

    :param d_2 (float): diameter at section 2 [m]
    :return r_cvt (float): curvature radius [m]
    """
    r_cvt = .06 * d_2

    return r_cvt


def streamline_diam(d_hu, d_0, theta=None, r_msl=None):
    """Calculate middle streamline diameter at a given angle.

    streamline_diam(d_hu, d_0):
        streamline diameter at section 0
    streamline_diam(d_hu, d_0, theta, r_msl):
        streamline diameter at i-section located by angle theta

    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :param theta (float): angle vertical axis and streamline curv. radius [rad]
    :param r_msl (float): middle streamline curvature radius [m]
    :return d_msl (float): middle streamline diameter [m]
    """
    d_msl = (d_0 + d_hu) / 2
    if theta is not None:
        d_msl = d_msl + (r_msl * (1 - math.cos(theta))) * 2

    return d_msl


def streamline_curv_rad(d_hu, d_0, r_cvt):
    """Calculate middle streamline curvature radius.

    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :param r_cvt (float): curvature radius [m]
    :return r_msl (float): middle streamline curvature radius [m]
    """
    r_msl = (d_0 - d_hu) / 4 + r_cvt

    return r_msl


def streamline_len(r_msl, d_2=None, d_msl=None, theta=None):
    """Calculate middle streamline length.

    streamline_len(r_msl, d_2, d_msl):
        middle streamline length from section 0 to 2
    streamline_len(r_msl, theta):
        middle streamline length from section 0 to i located by angle theta

    :param r_msl (float): middle streamline curvature radius [m]
    :param d_2 (float): diameter at section 2 [m]
    :param d_msl (float): middle streamline diameter [m]
    :param theta (float): angle vertical axis and streamline curv. radius [rad]
    :return l_msl (float): middle streamline length [m]
    """
    if theta is None:
        l_msl = math.pi / 2 * r_msl + (d_2 - d_msl - 2 * r_msl) / 2
    else:
        l_msl = theta * r_msl
    return l_msl


def area(l_imsl, l_msl, d_hu, d_0, d_2, b_2, x_2):
    """Calculate impeller vane area at i-section along the middle streamline.

    :param l_imsl (float): middle streamline length at i-section [m]
    :param l_msl (float): middle streamline length [m]
    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :param d_2 (float): diameter at section 2 [m]
    :param b_2 (list): impeller width at section 2 [m]
    :param x_2 (float): blade blockage at section 2
    :return a_i (float): area at i-section
    """
    a_0 = (d_0**2 - d_hu**2) * math.pi / 4
    a_2 = math.pi * d_2 * b_2 * x_2
    a_i = a_0 + (a_2 - a_0) * l_imsl / l_msl

    return a_i


def width(d, a=None, c_m=None, flow=None, x=None, eta_vol=1):
    """Calculate impeller width.

    width(d, u, flow, z, t, beta_b):
        width for a given flow rate
    width(d, a):
        width for a given area

    :param d (float): diameter [m]
    :param a (flaot): area [m^2]
    :param c_m (float): meridional component of the absolute velocity [m/s]
    :param flow (float): flow rate [m^3/s]
    :param t (float): blade thickness [m]
    :param z (int): number of blades
    :param beta_b (float): angle between rel. and blade velocity [m/s]
    :param eta_vol (float): volumetric efficency
    :return b (float): impeller width [m]
    """
    if a is None:
        b = flow / (math.pi * d * x * c_m * eta_vol)
    else:
        b = a / (math.pi * d)

    return b


def meridional_abs_vel(u, phi):
    """Calculate meridional component of the absolute velocity.

    :param u (float): blade velocity [m/s]
    :param phi (float): flow number
    :return c_m (float): meridional component of the absolute velocity [m/s]
    """
    c_m = phi * u

    return c_m


def circumferential_abs_vel(u, c_m, beta_b):
    """Calculate circumferential component of the absolute velocity.

    :param u (float): blade velocity [m/s]
    :param c_m (float): mmeridional component of the absolute velocity [m/s]
    :param beta_b (float): angle between rel. and blade velocity [m/s]
    :return c_u (float): circumferential component of the blade vel. [m/s]
    """
    c_u = u - c_m / math.tan(beta_b)

    return c_u


def blade_vel(omega, d):
    """Calculate blade velocity for a given angular velocity.

    :param omega (float): angular velocity [rad/s]
    :param d (float): diameter [m]
    :return u (float): absolute velocity [m/s]
    """
    u = omega * d / 2

    return u


def relative_vel(c_m, beta_b):
    """Calculate relative velocity.

    :param c_m (float): meridional component of the absolute velocity [m/s]
    :param beta_b (float): angle between rel. and blade velocity [m/s]
    :return w (float): relative velocity [m/s]
    """
    w = c_m / math.sin(beta_b)

    return w


def angle_theta(n, i):
    """Calculate angle between vertical axis and middle streamline.

    :param n (int): num. of divisions of the impeller vane curved line
    :param i (int): section
    :return theta (float): angle between vertical and middle streamline [rad]
    """
    theta = i * math.pi / (2 * (n - 1))

    return theta


def angle_gamma(r_cvt, r_msl, theta):
    """Claculate blade tilt angle between horizontal and center blade.

    The angle is determined so that each streamline in the impeller vane has
    the same length for a given angle between vertical axis and center blade.

    :param r_cvt (float): curvature radius [m]
    :param r_msl (float): middle streamline curvature radius [m]
    :param theta (float): angle between vertical and middle streamline [rad]
    :return gamma (float): angle between meridional blade vel. and vert. [rad]
    """
    l_arc = (r_msl - r_cvt) * (math.pi / 2 - theta)

    delta = l_arc / r_cvt

    l_height = r_cvt * math.sin(delta)

    b_half = r_msl - r_cvt * math.cos(delta)
    alpha = math.atan(l_height / b_half)
    if delta < theta:
        gamma = math.pi / 2 - theta - alpha
    else:
        gamma = None

    return gamma


def angle_beta(u, phi, eta_vol, x, gamma=0, psi_th=0, u_sf=0):
    """Return function in blade built angle as variable.

    :param u (float): blade velocity [m/s]
    :param phi (float): flow number
    :param eta_vol (float): volumetric efficency
    :param x (float): blade blockage
    :param gamma (float): angle between meridional blade vel. and vert. [rad]
    :param psi_th (float): theoretic head number
    :param u_sf (float): slip factor [m/s]
    :return (function): function in blade built angle as variable
    """
    return lambda beta_b: beta_b - math.atan(phi * math.cos(gamma) /
                                             (x * eta_vol *
                                             (1 - u_sf / u - psi_th)))


def slip_factor(u, beta_b, z):
    """Calculate slip factor with Wiesner's formula.

    :param u (float): absolute velocity [m/s]
    :param beta_b (float): angle between rel. and blade velocity [m/s]
    :param z (int): number of blades
    :return u_sf (float): slip factor [m/s]
    """
    u_sf = u * (math.sin(beta_b))**0.5 / z**0.7

    return u_sf


def degree_reaction(phi, beta_b, z):
    """Calculate degree of reaction.

    :param phi (float): flow number
    :param beta_b (float): angle between rel. and blade velocity [m/s]
    :param z (int): number of blades
    :return epsilon_ract (float): degree of reaction
    """
    epsilon_ract = 1 - (1 - phi / math.tan(beta_b) -
                        (math.sin(beta_b))**.5 / z**.7) / 2

    return epsilon_ract
