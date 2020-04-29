"""Methods to calculate the impeller."""

import math
import constants as CN


def type_number(omega, flow, head):
    """Calculate centrifugal pump's typical number.

    :param omega (float): angular velocity [rad/s]
    :param flow (float): flow rate [m^3/s]
    :param head (float): head [m]
    :return cappa (float): typical number
    """
    cappa = omega * flow**0.5 / (CN.G * head)**0.75

    return cappa


def cappa2rpm(cappa, flow, head):
    """Calculate rotational speed for a given typical number.

    :param cappa (float): typical number
    :param flow (float): flow rate [m^3/s]
    :param head (float): head [m]
    :return rpm (float): rotational speed [rpm]
    """
    omega = cappa * (CN.G * head)**0.75 / flow**0.5
    rpm = omega * 60 / (2 * math.pi)

    return rpm


def rpm2np(rpm, slip, hz):
    """Calculate number of poles of an AC motor for a given rotational speed.

    :param rpm (float): rotational speed [rpm]
    :param slip (int): slip factor [%]
    :param hz (int): utility frequency [Hz]
    :return np (int): number of poles
    """
    np = round(120 * hz / rpm * (1 - slip / 100))

    return np


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
    eta         .600 .750 .800 .890 .910 .920 .928 .929 .930 .929 .928
    weights     ones(cappa)
    n           3

    :param cappa (float): typical number
    :return eta (float): efficency
    """
    eta = 0.237 + 2.332 * cappa - 2.569 * cappa**2 + 0.924 * cappa**3

    return eta


def efficency_hyd_poly(cappa):
    """Calculate hydraulic efficency for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta_hyd     .600 .700 .750 .875 .895 .910 .913 .914 .915 .914 .913
    weights     ones(cappa)
    n           3

    :param cappa (float): typical number
    :return eta_hyd (float): hydraulic efficency
    """
    eta_hyd = 0.268 + 1.989 * cappa - 1.986 * cappa**2 + 0.646 * cappa**3

    return eta_hyd


def efficency_vol_poly(cappa):
    """Calculate volumetric efficency for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta_hyd     .910 .940 .950 .953 .955 .958 .960 .963 .965 .968 .970
    weights     ones(cappa)
    n           3

    :param cappa (float): typical number
    :return eta_vol (float): volumetric efficency
    """
    eta_vol = 0.854 + 0.390 * cappa - 0.476 * cappa**2 + 0.195 * cappa**3

    return eta_vol


def flow_number_poly(cappa):
    """Calculate flow number for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    phi         .080 .093 .100 .110 .120 .130 .140 .150 .160 .165 .170
    weights     ones(cappa)
    n           3

    :param cappa (float): typical number
    :return phi (float): flow number
    """
    phi = 0.0675 + 0.0557 * cappa + 0.0839 * cappa**2 - 0.0490 * cappa**3

    return phi


def flow_number(d, b, u, x, flow, eta_vol):
    """Calculate flow number.

    :param d (float): diameter [m]
    :param b (float): impeller width [m]
    :param u (float): absolute velocity [m/s]
    :param x (float): blade blockage
    :param flow (float): flow rate [m^3/s]
    :param eta_vol (float): volumetric efficency
    :return phi (float): flow number
    """
    phi = flow / (math.pi * d * b * u * x * eta_vol)

    return phi


def theoretic_flow_number(phi, x, eta_vol):
    """Calculate theoretic flow number.

    :param phi (float): flow number
    :param x (float): blade blockage
    :param eta_vol (float): volumetric efficency
    :return phi_th (float): theoretic flow number
    """
    phi_th = phi / (x * eta_vol)

    return phi_th


def head_number_poly(cappa):
    """Calculate head number for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    psi         .55 .54 .53 .52 .51 .49 .45 .43 .41 .40 .38
    weights     ones(cappa)
    n           3

    :param cappa (float): typical number
    :return psi (float): head number
    """
    psi = 0.520 + 0.237 * cappa - 0.595 * cappa**2 + 0.251 * cappa**3

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


def blade_blockage(beta_c, d, thk, z):
    """Calculate blade blockage.

    :param beta_c (float): angle between relative and blade velocity [m/s]
    :param d (float): diameter [m]
    :param thk (float): blade thickness [m]
    :param z (int): number of blades
    :return x (float): blade blockage factor
    """
    x = 1 - (z * thk) / (math.pi * d * math.sin(beta_c))

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
    :return r_c (float): curvature radius [m]
    """
    r_c = .06 * d_2

    return r_c


def streamline_diam(d_hu, d_0, theta=None, r_slc=None):
    """Calculate middle streamline diameter at a given angle.

    streamline_diam(d_hu, d_0):
        streamline diameter at section 0
    streamline_diam(d_hu, d_0, theta, r_slc):
        streamline diameter at i-section located by angle theta

    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :param theta (float): angle vertical axis and streamline curv. radius [rad]
    :param r_slc (float): streamline curvature radius [m]
    :return d_sl (float): middle streamline diameter [m]
    """
    d_sl = (d_0 + d_hu) / 2
    if theta is not None:
        d_sl = d_sl + (r_slc * (1 - math.cos(theta))) * 2

    return d_sl


def streamline_curv_rad(d_hu, d_0, r_c):
    """Calculate middle streamline curvature radius.

    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :param r_c (float): curvature radius [m]
    :return r_slc (float): streamline curvature radius [m]
    """
    r_slc = (d_0 - d_hu) / 4 + r_c

    return r_slc


def streamline_len(r_slc, d_2=None, d_sl=None, theta=None):
    """Calculate middle streamline length.

    streamline_len(r_slc, d_2, d_sl):
        middle streamline length from section 0 to 2
    streamline_len(r_slc, d_2, d_sl, theta):
        middle streamline length from section 0 to i located by angle theta

    :param r_slc (float): streamline curvature radius [m]
    :param d_2 (float): diameter at section 2 [m]
    :param d_sl (float): middle streamline diameter [m]
    :param theta (float): angle vertical axis and streamline curv. radius [rad]
    :return l_mid (float): middle streamline length [m]
    """
    if theta is None:
        l_sl = math.pi / 2 * r_slc + (d_2 - d_sl - 2 * r_slc) / 2
    else:
        l_sl = theta * r_slc
    return l_sl


def area(l_isl, l_sl, d_hu, d_0, d_2, b_2, x_2):
    """Calculate impeller vane area at i-section along the middle streamline.

    :param l_isl (float): middle streamline length at i-section [m]
    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :param d_2 (float): diameter at section 2 [m]
    :param b_2 (list): impeller width at section 2 [m]
    :param x_2 (float): blade blockage at section 2
    :return a_i (float): area at i-section
    """
    a_0 = (d_0**2 - d_hu**2) * math.pi / 4
    a_2 = math.pi * d_2 * b_2 * x_2
    a_i = a_0 + (a_2 - a_0) * l_isl / l_sl

    return a_i


def width(d_isl, a_i=None, u_2=None, phi=None, flow=None, x_2=1, eta_vol=1):
    """Calculate impeller width at i-section.

    width(d_isl, u_2, phi, flow, x_2, eta_vol):
        width for a given flow rate
    width(d_isl, a_i):
        width for a given area

    :param d_isl (float): diameter at i-section along middle streamline [m]
    :param a_i (flaot): area at i-section [m^2]
    :param u_2 (float): blade velocity at section 2 [m/s]
    :param phi (float): flow coefficient
    :param flow (float): flow rate [m^3/s]
    :param x_2 (float): blade blockage at section 2
    :param eta_vol (float): volumetric efficency
    :return b_i (float): impeller width at i-section [m]
    """
    if a_i is None:
        b_i = flow / (math.pi * d_isl * u_2 * phi * x_2 * eta_vol)
    else:
        b_i = a_i / (math.pi * d_isl)

    return b_i


def meridional_abs_vel(u, phi_th):
    """Calculate meridional component of the absolute velocity.

    :param u (float): absolute velocity [m/s]
    :param phi_th (float): theoretic flow number
    :return c_m (float): meridional component of the absolute velocity [m/s]
    """
    c_m = phi_th * u

    return c_m


def circumferential_abs_vel(u, c_m, beta_c):
    """Calculate circumferential component of the absolute velocity.

    :param u (float): absolute velocity [m/s]
    :param c_m (float): mmeridional component of the absolute velocity [m/s]
    :param beta_c (float): angle between rel. and blade velocity [m/s]
    :return c_u (float): circumferential component of the abs. vel. [m/s]
    """
    c_u = u - c_m / math.tan(beta_c)

    return c_u


def blade_vel(omega, d):
    """Calculate blade velocity for a given angular velocity.

    :param omega (float): angular velocity [rad/s]
    :param d (float): diameter [m]
    :return u (float): absolute velocity [m/s]
    """
    u = omega * d / 2

    return u


def relative_vel(c_m, beta_c):
    """Calculate relative velocity.

    :param c_m (float): meridional component of the absolute velocity [m/s]
    :param beta_c (float): angle between rel. and blade velocity [m/s]
    :return w (float): relative velocity [m/s]
    """
    w = c_m / math.sin(beta_c)

    return w


def angle_theta(n, i):
    """Calculate angle between vertical axis and middle streamline.

    :param n (int): num. of divisions of the impeller vane curved line
    :param i (int): section
    :return theta (float): angle between vertical and middle streamline [rad]
    """
    theta = i * math.pi / (2 * (n - 1))

    return theta


def angle_gamma(r_c, r_slc, theta):
    """Claculate blade tilt angle between horizontal and center blade.

    The angle is determined so that each streamline in the impeller vane has
    the same length for a given angle between vertical axis and center blade.

    :param r_c (float): curvature radius [m]
    :param r_slc (float): streamline curvature radius [m]
    :param theta (float): angle between vertical and middle streamline [rad]
    :return gamma (float): angle between meridional abs. vel. and vert. [rad]
    """
    l_arc = (r_slc - r_c) * (math.pi / 2 - theta)

    delta = l_arc / r_c

    l_chord = r_c * math.sin(delta)

    b_half = r_slc - r_c * math.cos(delta)
    alpha = math.atan(l_chord / b_half)
    if delta < theta:
        gamma = math.pi / 2 - theta - alpha
    else:
        gamma = None

    return gamma


def angle_beta(u, c_m, gamma=0, psi_th=0, u_sf=0):
    """Calculate blade construction angle between rel. and blade velocity vect.

    :param u (float): blade velocity [m/s]
    :param c_m (float): meridional component of the absolute velocity [m/s]
    :param gamma (float): angle between meridional abs. vel. and vert. [rad]
    :param psi_th (float): theoretic head number
    :param u_sf (float): slip factor
    :return beta (float): angle between rel. and blade velocity vectors [rad]
    """
    beta = math.atan(c_m * math.cos(gamma) / (u * (1 - psi_th) - u_sf))

    return beta


def slip_factor(u, beta_c, z):
    """Calculate slip factor with Wiesner's formula.

    :param u (float): absolute velocity [m/s]
    :param beta_c (float): angle between rel. and blade velocity [m/s]
    :param z (int): number of blades
    :return u_sf (float): slip factor [m/s]
    """
    u_sf = u * (math.sin(beta_c))**0.5 / z**0.7

    return u_sf


def degree_reaction(phi_th, beta_c, z):
    """Calculate degree of reaction.

    :param phi_th (float): flow coefficient corrected
    :param beta_c (float): angle between rel. and blade velocity [m/s]
    :param z (int): number of blades
    :return epsilon_ract (float): degree of reaction
    """
    epsilon_ract = 1 - (1 - phi_th / math.tan(beta_c) -
                        (math.sin(beta_c))**.5 / z**.7) / 2

    return epsilon_ract
