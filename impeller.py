"""Methods to calculate the impeller."""

import math
import constants as CNST


def type_number(omega, flow, head):
    """Calculate centrifugal pump's typical number.

    :param omega (float): angular velocity [rad/s]
    :param flow (float): flow rate [m^3/s]
    :param head (float): head [m]
    :return cappa (float): typical number
    """
    cappa = omega * flow**0.5 / (CNST.G * head)**0.75

    return cappa


def cappa2rpm(cappa, flow, head):
    """Calculate rotational speed for a given typical number.

    :param cappa (float): typical number
    :param flow (float): flow rate [m^3/s]
    :param head (float): head [m]
    :return rpm (float): rotational speed [rpm]
    """
    omega = cappa * (CNST.G * head)**0.75 / flow**0.5
    rpm = omega * 60 / (2 * math.pi)

    return rpm


def rpm2pp(rpm, slip, hz):
    """Calculate polar pairs of an AC motor for a given rotational speed.

    :param rpm (float): rotational speed [rpm]
    :param slip (int): slip factor [%]
    :param hz (int): utility frequency [Hz]
    :return pp (int): polar pairs
    """
    pp = 120 * hz / rpm * (1 - slip / 100)

    return pp


def flow_number_poly(cappa):
    """Calculate flow number for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    phi         .080 .093 .100 .110 .120 .130 .140 .150 .160 .165 .170
    weights     ones(cappa)
    n           2

    :param cappa (float): typical number
    :return phi (float): flow number
    """
    phi = 0.0567636 + 0.118979 * cappa - 0.0188811 * cappa**2

    return phi


def head_number_poly(cappa):
    """Calculate head number for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    psi         .55 .54 .53 .52 .51 .49 .45 .43 .41 .40 .38
    weights     ones(cappa)
    n           2

    :param cappa (float): typical number
    :return psi (float): head number
    """
    psi = 0.5747273 - 0.0864569 * cappa - 0.0687646 * cappa**2

    return psi


def efficency_poly(cappa):
    """Calculate efficency for a given pump's typical number.

    The polynomial has been calculated applaying the curve fitting at nodes
    cappa       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
    eta         0 .650 .800 .890 .910 .920 .928 .929 .930 .929 .928
    weights     ones(cappa)
    n           2

    :param cappa (float): typical number
    :return eta_coef (float): efficency coefficients
    """
    eta = - 0.2774 + 3.0137016 * cappa - 1.7473193 * cappa**2

    return eta


def blade_vel_psi(psi, head):
    """Calculate peripheral velocity function of the head number.

    :param psi (float): head number
    :param head (float): head [m]
    :return u (float): peripheral velocity [m/s]
    """
    u_psi = (CNST.G * head / psi)**0.5

    return u_psi


def diameter_rpm(u, rpm):
    """Calculate diameter function of the rotationl speed.

    :param u (float): peripheral velocity [m/s]
    :param rpm (float): rotational speed [rpm]
    :return d (float): diameter [m]
    """
    d_rpm = (60 * u) / (math.pi * rpm)

    return d_rpm


def width_phi(u, d, phi, flow):
    """Calculate impeller width function of the flow number.

    :param u (float): peripheral velocity [m/s]
    :param d (float): diameter [m]
    :param phi (float): flow number
    :param flow (float): flow rate [m^3/s]
    :return b (float): impeller width [m]
    """
    b = flow / (math.pi * d * u * phi)

    return b


def width0diameter(b, d):
    """Calculate rate width over diameter.

    :param b (float): impeller width [m]
    :param d (float): diameter [m]
    :return bd (float): width over diameter
    """
    bd = b / d

    return bd


def npsh_req(cappa, head):
    """Calculate the neat positive suction head required.

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
    :return x (float): hub blockage
    """
    x = 1 - (d_hu / d)**2

    return x


def blade_blockage(beta_c, d, thk, z):
    """Calculate blade blockage.

    :param beta_c (float): angle between rel. and blade velocity [m/s]
    :param d (float): diameter [m]
    :param thk (float): blade thickness [m]
    :param z (int): number of blades
    :return x (float): blade blockage
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
    :param lm (float): loss coefficient at section 0
    :param lw (float): low-pressure peak coefficient at blades at section 0
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
    :param km (float): rate between circumeferential velocity cm2 and c0
    :param eta_vol (float): volumetric efficency
    :return d_eff (float): diameter with max efficency [m]
    """
    d_eff = 2 * ((2 * flow**2 * km**2) /
                 (eta_vol**2 * math.pi**2 * omega**2 * x**2))**(1/6)

    return d_eff


def diameter_flow(omega, x, flow, eta_vol):
    """Calculate diameter function of the flow rate.

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
    for i in range(len(CNST.D_INT)):
        dif_abs.append(abs(CNST.D_INT[i] - d))
    dif_min = min(dif_abs)
    d_std = CNST.D_INT[dif_abs.index(dif_min)]

    return d_std


def streamline_diam(d_hu, d_0):
    """Calculate middle streamline diameter at section 0.

    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :return d_mid (float): middle streamline diameter [m]
    """
    d_mid = (d_0 + d_hu) / 2

    return d_mid


def curvature_rad(d_1):
    """Calculate curvature radius of the shroud at section 0.

    :param d_1 (float): diameter [m]
    :return r_cvt (float): curvature radius [m]
    """
    r_cvt = .06 * d_1

    return r_cvt


def streamline_rad(d_hu, d_0, r_cvt):
    """Calculate curvature radius of the middle streamline at section 0.

    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :param r_cvt (float): curvature radius [m]
    :return r_mid (float): streamline radius [m]
    """
    r_mid = (d_0 - d_hu) / 4 + r_cvt

    return r_mid


def streamline_len(r_cvt, r_mid, d_1, d_0):
    """Calculate length middle streamline.

    :param r_cvt (float): curvature radius [m]
    :param r_mid (float): streamline radius [m]
    :param d_1 (float): diameter at section 1 [m]
    :param d_0 (float): diameter at section 0 [m]
    :return l_mid (float): middle streamline length [m]
    """
    l_mid = math.pi / 2 * r_mid + (d_1 - d_0 - 2 * r_cvt) / 2

    return l_mid


def angle_theta_2(r_cvt, r_mid, d_1, d_0, d_2):
    """Calculate angle between radius middle streamline and vertical axis
    at section 2.

    :param r_cvt (float): curvature radius [m]
    :param r_mid (float): streamline radius [m]
    :param d_1 (float): diameter at section 1 [m]
    :param d_0 (float): diameter at section 0 [m]
    :param d_2 (float): diameter at section 2 [m]
    :return theta_2 (float): angle between streamline rad. and vert. [rad]
    """
    theta_2 = math.acos((d_1 - d_0 - 2 * r_cvt - d_2) / (2 * r_mid))

    return theta_2


def width_2(a_0, a_1, r_mid, l_mid, theta_2, d_2):
    """Calculate impeller width at section 2.

    :param a_0 (float): area at section 0 [m^2]
    :param a_1 (float): area at section 1 [m^2]
    :param r_mid (float): radius mean streamline [m]
    :param l_mid (float): length mean streamline [m]
    :param theta_2 (float): angle between streamline rad. and vert. [rad]
    :param d_2 (float): diameter at section 2 [m]
    :return b_2 (float): impeller width at section 2[m]
    """
    b_2 = (a_0 + (a_1 - a_0) * (r_mid * theta_2) / l_mid) / (math.pi * d_2)

    return b_2


def width_1(d_1, u_1, phi, x_1, flow, eta_vol):
    """Calculate impeller width at section 1.

    :param d_1 (float): diameter at section 1 [m]
    :param u_1 (float): absolute velocity at section 1 [m/s]
    :param phi (float): flow coefficient
    :param flow (float): flow rate [m^3/s]
    :param eta_vol (float): volumetric efficency
    :return b_1 (float): impeller width at section 1
    """
    b_1 = flow / (math.pi * d_1 * u_1 * phi * x_1 * eta_vol)

    return b_1


def area_0(d_hu, d_0):
    """Calculate area at section 0.

    :param d_hu (float): hub diameter [m]
    :param d_0 (float): diameter at section 0 [m]
    :return a0 (float): area [m^2]
    """
    a_0 = (d_0**2 - d_hu**2) * math.pi / 4

    return a_0


def area_1(d1, b1, x1):
    """Calculate area at section 1.

    :param d1 (float): diameter [m]
    :param b1 (list): impeller width [m]
    :param x1 (float): blade blockage
    :return a1 (float): area [m^2]
    """
    a1 = math.pi * d1 * b1 * x1

    return a1


def width_impeller_vane(r_mid, l_mid, d_int, a0, a1):
    """Calculate width impeller vane as different diameters
    at equals angles in the curved zone.

    :param r_mid (float): radius mean streamline [m]
    :param l_mid (float): length mean streamline [m]
    :param d_int (float): center radius curvature front shroud [m]
    :param a0 (float): area [m^2]
    :param a1 (float): area [m^2]
    :return b_im (list): diameters at equals angles in the curved zone [m]
    """
    # in the curved zone, the mean streamline is a quarter of arc of circle
    # and along its path are considered:
    n = 11  # num. divisions: points = segments + 1
    theta = []  # cumulative angles
    length = []  # cumulative lengths
    area = []  # area at different lengths
    r = []  # radius at different angles
    b_im = []
    for i in range(n):
        theta.append((math.pi / (2 * n)) * i)
        length.append(r_mid * theta[i])
        area.append(a0 + (a1 - a0) * length[i] / l_mid)
        r.append(d_int / 2 - r_mid * math.cos(theta[i]))
        b_im.append(area[i] / (2 * math.pi * r[i]))

    return list(zip(theta, length, area, r, b_im))


def meridional_abs_vel(b, d, x, flow, eta_vol):
    """Calculate meridional velocity component of the absolute velocity.

    :param b (float): impeller width [m]
    :param d (float): diameter [m]
    :param x (float): blade blockage
    :param flow (float): flow rate [m^3/s]
    :param eta_vol (float): volumetric efficency
    :return c_m (float): meridional velocity component [m/s]
    """
    c_m = flow / (math.pi * d * b * x * eta_vol)

    return c_m


def circumferential_abs_vel(u, c_m, beta_c):
    """Calculate circumferential velocity component of the absolute
    velocity.

    :param u (float): absolute velocity [m/s]
    :param c_m (float): meridional velocity component [m/s]
    :param beta_c (float): angle between rel. and blade velocity [m/s]
    :return c_u (float): absolute velocity component [m/s]
    """
    c_u = u - c_m * 1 / math.tan(beta_c)

    return c_u


def blade_vel(omega, d):
    """Calculate blade velocity function of angular velocity.

    :param omega (float): angular velocity [rad/s]
    :param d (float): diameter [m]
    :return u (float): absolute velocity [m/s]
    """
    u = omega * d / 2

    return u


def relative_vel(c_m, beta_c):
    """Calculate relative velocity.

    :param c_m (float): meridional velocity [m/s]
    :param beta_c (float): angle between rel. and blade velocity [m/s]
    :return w (float): relative velocity [m/s]
    """
    w = c_m / math.sin(beta_c)

    return w


def angle_beta_2c(cm2, u2, gamma_2):
    """Calculate blade working angle between relative and blade
    velocity at section 2.

    :param cm2 (float): meridional velocity [m/s]
    :param u2 (float): absolute velocity [m/s]
    :param gamma_2 (int): measured angle between cm2 and vertical [deg]
    :return beta_2c (float): angle between rel. and circum. velocity [m/s]
    """
    gammar_2 = math.radians(gamma_2)
    beta_2c = math.atan((cm2 * math.cos(gammar_2)) / u2)

    return beta_2c


def angle_beta_1c(psi_th, phi_th, u_1sf, u_1):
    """Calculate blade working angle between relative and blade
    velocity at section 1.

    :param psi_th (float): theoretic head coefficient
    :param phi_th (float): flow coefficient corrected
    :param u_1sf (float): slip factor at section 1 [m/s]
    :param u_1 (float): absolute velocity at section 1 [m/s]
    :return beta_1c (float): angle between rel. and blade velocity [m/s]
    """
    beta_1c = math.atan(phi_th / (1 - psi_th - u_1sf / u_1))

    return beta_1c


def head_number(u, head):
    """Calculate head number.

    :param u (float): absolute velocity [m/s]
    :param head (float): head [m]
    :return psi (float): head coefficient
    """
    psi = (CNST.G * head) / u**2

    return psi


def flow_number(d, b, u, x, flow, eta_vol):
    """Calculate flow number.

    :param d (float): diameter [m]
    :param b (float): impeller width [m]
    :param u (float): absolute velocity [m/s]
    :param x (float): blade blockage
    :param flow (float): flow rate [m^3/s]
    :param eta_vol (float): volumetric efficency
    :return phi (float): flow coefficient
    """
    phi = flow / (math.pi * d * b * u * x * eta_vol)

    return phi


def theoretic_head_number(psi, eta_hyd):
    """Calculate theoretic head coefficient.

    :param psi (float): head coefficient
    :param eta_hyd (float): hydraulic efficency
    :return psi_th (float): theoretic head coefficient
    """
    psi_th = psi / eta_hyd

    return psi_th


def theoretic_flow_number(phi, x, eta_vol):
    """Calculate theoretic flow coefficient.

    :param phi (float): flow coefficient
    :param x (float): blade blockage
    :param eta_vol (float): volumetric efficency
    :return phi_th (float): flow coefficient corrected
    """
    phi_th = phi / (x * eta_vol)

    return phi_th


def slip_factor(u, beta_c, z):
    """Calculate slip factor with Wiesner"s formula.

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
