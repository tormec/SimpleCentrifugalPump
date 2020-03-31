#!/usr/bin/env python3
"""Calculation of geometrical dimensions for centrifugal pump."""

import math

# name sections
# 0: impeller eye
# 1: impeller blade trailing edge
# 2: impeller blade leading edge
# 3: volute start wrap angle

# physical constants
G = 9.81  # gravity acceleration [m/s^2]
RHO = 1000  # water density

# couple poles for electric motor
CPOLES = [2, 4, 6, 8]

# internal diameters for commercial pipes
D_INT = [.054, .070, .082, .107, .132, .159, .207, .260, .310, .340]  # [m]


class Pre_Values(object):
    """Methods for feasability study of the impeller."""

    def rotational_speed(self, slip, hz, cp):
        """Calculate rotational speed at diff. couple poles.

        :param slip (int): slip of electric induction motor [%]
        :param hz (int): utility frequency [Hz]
        :param cp (int): couple of poles of electric induction motor
        :return rpm (float): rotational speed [rpm]
        """
        rpm = 120 * hz / cp * (1 - slip / 100)

        return rpm

    def angular_velocity(self, rpm):
        """Calculate angular velocity pump shaft.

        :param rpm (float): rotational speed [rpm]
        :return omega (float): angular velocity [rad/s]
        """
        omega = 2 * math.pi * rpm / 60

        return omega

    def type_number(self, rpm, flow, head):
        """Calculate centrifugal pump's typical number.

        :param rpm (float): rotational speed [rpm]
        :param flow (float): flow rate [m^3/s]
        :param head (float): head [m]
        :return k_num (float): typical number
        """
        omega = self.angular_velocity(rpm)
        k_num = omega * flow**0.5 / (G * head)**0.75

        return k_num

    def k2n(self, k_num, flow, head):
        """Calculate the rotational speed for a given typical number.

        :param k_num (float): typical number
        :param flow (float): flow rate [m^3/s]
        :param head (float): head [m]
        :return rpm (float): rotational speed [rpm]
        """
        omega = k_num * (G * head)**0.75 / flow**0.5
        rpm = omega * 60 / (2 * math.pi)

        return rpm

    def n2cp(self, rpm, slip, hz):
        """Calculate the couple of poles of electric induction motor
        for a given rotational speed.

        :param rpm (float): rotational speed [rpm]
        :param slip (int): slip of the electric induction motor [%]
        :param hz (int): utility frequency [Hz]
        :return cp (int): couple of poles of electric induction motor
        """
        cp = 120 * hz / rpm * (1 - slip / 100)

        return cp

    def flow_number(self, k_num):
        """Calculate value of flow number for a given pump's typical number.

        The polynomial has been calculated applaying the curve fitting at nodes
        k_num       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
        phi         .080 .093 .100 .110 .120 .130 .140 .150 .160 .165 .170
        weights     ones(k_num)
        n           2

        :param k_num (float): typical number
        :return phi (float): flow number
        """
        phi = 0.0567636 + 0.118979 * k_num - 0.0188811 * k_num**2

        return phi

    def head_number(self, k_num):
        """Calculate value of head number for a given pump's typical number.

        The polynomial has been calculated applaying the curve fitting at nodes
        k_num       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
        psi         .55 .54 .53 .52 .51 .49 .45 .43 .41 .40 .38
        weights     ones(k_num)
        n           2

        :param k_num (float): typical number
        :return psi (float): head number
        """
        psi = 0.5747273 - 0.0864569 * k_num - 0.0687646 * k_num**2

        return psi

    def efficency(self, k_num):
        """Calculate value of efficency for a given pump's typical number.

        The polynomial has been calculated applaying the curve fitting at nodes
        k_num       .2 .3 .4 .5 .6 .7 .8 .9 1.0 1.1 1.2
        eta         0 .650 .800 .890 .910 .920 .928 .929 .930 .929 .928
        weights     ones(k_num)
        n           2

        :param k_num (float): typical number
        :return eta_coef (float): efficency coefficients
        """
        eta = - 0.2774 + 3.0137016 * k_num - 1.7473193 * k_num**2

        return eta

    def _circumferential_velocity_1(self, head, psi):
        """Calculate circumferential velocity at section 1.

        :param head (float): head [m]
        :param psi (float): head number
        :return u1 (float): circumferential velocity [m/s]
        """
        u1 = (G * head / psi)**0.5

        return u1

    def _diameter_1(self, u1, rpm):
        """Calculate diameter at section 1.

        :param u1 (float): circumferential velocity [m/s]
        :param rpm (float): rotational speed [rpm]
        :return d1 (float): diameter [m]
        """
        d1 = (60 * u1) / (math.pi * rpm)

        return d1

    def _width_1(self, u1, d1, flow, phi):
        """Calculate impeller width at section 1.

        :param u1 (float): circumferential velocity [m/s]
        :param d1 (float): diameter [m]
        :param flow (float): flow rate [m^3/s]
        :param phi (float): flow number
        :return b1 (float): impeller width [m]
        """
        b1 = flow / (math.pi * d1 * u1 * phi)

        return b1

    def width_over_diameter_1(self, b1, d1):
        """Calculate rate width over diameter at section 1.

        :param b1 (float): impeller width [m]
        :param d1 (float): diameter [m]
        :return bd1 (float): width over diameter
        """
        bd1 = b1 / d1

        return bd1

    def npsh_r(self, k_num, head):
        """Calculate the neat positive suction head required.

        :param k_num (float): typical number
        :param head (float): head [m]
        :return npsh_r (float): neat positive suction head required [m]
        """
        npsh_r = .25 * k_num**(4/3) * head

        return npsh_r


class Shaft(object):
    """Methods for pump shaft design."""

    def power(self, eta_tot, flow, head):
        """Calculate power the el. motor should provide at pump shaft.

        :param eta_t (float): total efficency
        :param flow (float): flow rate [m^3/s]
        :param head (float): head [m]
        :return powr (float): power [W]
        """
        powr = (RHO * flow * G * head) / eta_tot

        return powr

    def torque(self, powr, omega):
        """Calculate torque at pump shaft.

        :param powr (float): power [W]
        :param omega (float): angular velocity [rad/s]
        :return torq (float): torque [Nm]
        """
        torq = powr / omega

        return torq

    def diameter_shaft(self, torq, tau_adm):
        """Calculate shaft diameter.

        :param torq (float): torque [Nm]
        :param tau_adm (int): tau admissible [MPa]
        :return d_sh (float): shaft diameter [m]
        """
        d_sh = ((32 * torq) / (math.pi * (tau_adm * 10**6)))**(1/3)

        return d_sh

    def diameter_hub(self, d_sh):
        """Calculate hub diameter.

        :param d_sh (float): shaft diameter [m]
        :return d_hu (float): hub diameter [m]
        """
        d_hu = 1.5 * d_sh

        return d_hu


class Impeller(object):
    """Methods for impeller design."""

    def hub_blockage_0(self, d0, d_hu):
        """Calculate hub blockage at section 0.

        :param d0 (float): diameter [m]
        :param d_hu (float): hub diameter [m]
        :return x0 (float): hub blockage
        """
        x0 = 1 - (d_hu / d0)**2

        return x0

    def blade_blockage_1(self, beta_1c, d1, thk, z):
        """Calculate blade blockage at section 1.

        :param beta_1c (float): angle between rel. and circum. velocity [m/s]
        :param d1 (float): diameter [m]
        :param thk (float): blade thickness [m]
        :param z (int): number of blades
        :return x1 (float): blade blockage
        """
        x1 = 1 - (z * thk) / (math.pi * d1 * math.sin(beta_1c))

        return x1

    def blade_blockage_2(self, beta_2c, thk, z, d2):
        """Calculate blade blockage at section 2.

        :param beta_2c (float): angle between rel. and circum. velocity [m/s]
        :param thk (float): blade thickness
        :param z (int): number of blades
        :param d2 (float): measured diameter
        :return x2 (float): blade blockage
        """
        x2 = 1 - (z * thk) / (math.pi * d2 * math.sin(beta_2c))

        return x2

    def diameter_1(self, omega, u1):
        """Calculate diameter at section 1.

        :param omega (float): angular velocity [rad/s]
        :param u1 (float): circumferential velocity [m/s]
        :return d1 (float): diameter [m]
        """
        d1 = 2 * u1 / omega

        return d1

    def diameter_0_npsh(self, omega, x0, flow, lm, lw, km, eta_vol):
        """Calculate diameter at section 0 with min npsh_r.

        :param omega (float): angular velocity [rad/s]
        :param x0 (float): hub blockage
        :param flow (float): flow rate [m^3/s]
        :param lm (float): loss coefficient at section 0
        :param lw (float): low-pressure peak coefficient at blades at section 0
        :param km (float): rate between circumeferential velocity cm2 and c0
        :param eta_vol (float): volumetric efficency
        :return d0_npsh (float): diameter with min npsh_r [m]
        """
        d0_npsh = 2 * (
                       (2 * flow**2 * km**2 * (1 + lm + lw)) /
                       (eta_vol**2 * math.pi**2 * omega**2 * x0**2 * lw)
                      )**(1/6)

        return d0_npsh

    def diameter_0_efficency(self, omega, x0, flow, km, eta_vol):
        """Calculate diameter at section 0 with max total efficency.

        :param omega (float): angular velocity [rad/s]
        :param x0 (float): hub blockage
        :param flow (float): flow rate [m^3/s]
        :param km (float): rate between circumeferential velocity cm2 and c0
        :param eta_vol (float): volumetric efficency
        :return d0_eff (float): diameter with max efficency [m]
        """
        d0_eff = 2 * (
                      (2 * flow**2 * km**2) /
                      (eta_vol**2 * math.pi**2 * omega**2 * x0**2)
                     )**(1/6)

        return d0_eff

    def diameter_0_compromise(self, omega, x0, flow, eta_vol):
        """Calculate diameter at section 0 as compromise solution between
        min npsh_r and max total efficency.

        :param omega (float): angular velocity [rad/s]
        :param x0 (float): hub blockage
        :param flow (float): flow rate [m^3/s]
        :param eta_vol (float): volumetric efficency
        :param eta_vol (float): volumetric efficency
        :return d0_cmp (float): diameter as compromise solution [m]
        """
        d0_cmp = (
                 (flow * 8 * 3.03) /
                 (omega * math.pi * x0 * eta_vol)
                 )**(1/3)

        return d0_cmp

    def diameter_0(self, d0_npsh, d0_eff, d0_cmp):
        """Calculate diameter at section 0 according to standard diameters.

        :param d0_npsh (float): diameter with min NPSH_r [m]
        :param d0_eff (float): diameter with max efficency [m]
        :param d0_cmp (float): diameter as compromise solution [m]
        :return d0 (float): diameter [m]
        """
        d0_val = [d0_npsh, d0_eff, d0_cmp]
        d0_avg = math.fsum(d0_val) / len(d0_val)
        dif_abs = []
        for i in range(len(D_INT)):
            dif_abs.append(abs(D_INT[i] - d0_avg))
        dif_min = min(dif_abs)
        d0 = D_INT[dif_abs.index(dif_min)]

        return d0

    def diameter_mean_streamline(self, d_hu, d0):
        """Calculate mean streamline diameter.

        :param d_hu (float): hub diameter [m]
        :param d0 (float): diameter [0]
        :return d_mid (float): mean streamline diameter [m]
        """
        d_mid = (d0 + d_hu) / 2

        return d_mid

    def radius_curvature_front_shroud(self, d1):
        """Calculate radius curvature front shroud.

        :param d1 (float): diameter [m]
        :return r_cvt (float): radius curvature front shroud [m]
        """
        r_cvt = .06 * d1

        return r_cvt

    def center_mean_streamline(self, r_cvt, d0):
        """Calculate center radius curvature of the front shroud.

        :param r_cvt (float): radius curvature front shroud [m]
        :param d0 (float): diameter [0]
        :return d_int (float): center radius curvature front shroud [m]
        """
        d_int = d0 + 2 * r_cvt

        return d_int

    def angle_theta_2(self, d_int, r_mid, d2):
        """Calculate angle between radius mean streamline and vertical
        at section 2.

        :param d_int (float): center radius curvature front shroud [m]
        :param r_mid (float): radius mean streamline [m]
        :param d2 (float): measured diameter
        :return theta_2 (float): angle between radius m. stream. and vert.[rad]
        """
        theta_2 = math.acos((d_int - d2) / (2 * r_mid))

        return theta_2

    def width_2(self, a0, a1, r_mid, l_mid, theta_2, d2):
        """Calculate impeller width at section 2.

        :param a0 (float): area [m^2]
        :param a1 (float): area [m^2]
        :param r_mid (float): radius mean streamline [m]
        :param l_mid (float): length mean streamline [m]
        :param theta_2 (float): angle between radius m. stream. and vert. [rad]
        :param d2 (float): measured diameter
        :return b2 (float): impeller width [m]
        """
        b2 = (a0 + (a1 - a0) * (r_mid * theta_2) / l_mid) / (math.pi * d2)

        return b2

    def width_1(self, d1, u1, phi, x1, flow, eta_vol):
        """Calculate impeller width at section 1.

        :param d1 (float): diameter [m]
        :param u1 (float): circumferential velocity [m/s]
        :param phi (float): flow coefficient
        :param flow (float): flow rate [m^3/s]
        :param eta_vol (float): volumetric efficency
        :return x1 (float): blade blockage
        """
        b1 = flow / (math.pi * d1 * u1 * phi * x1 * eta_vol)

        return b1

    def area_0(self, d_hu, d0):
        """Calculate area at section 0.

        :param d_hu (float): hub diameter [m]
        :param d0 (float): diameter [0]
        :return a0 (float): area [m^2]
        """
        a0 = (d0**2 - d_hu**2) * math.pi / 4

        return a0

    def area_1(self, d1, b1, x1):
        """Calculate area at section 1.

        :param d1 (float): diameter [m]
        :param b1 (list): impeller width [m]
        :param x1 (float): blade blockage
        :return a1 (float): area [m^2]
        """
        a1 = math.pi * d1 * b1 * x1

        return a1

    def radius_mean_streamline(self, d_int, d_mid):
        """Calculate radius mean streamline.

        :param d_int (float): center radius curvature front shroud [m]
        :param d_mid (float): mean streamline diameter [m]
        :return r_mid (float): radius mean streamline [m]
        """
        r_mid = (d_int - d_mid) / 2

        return r_mid

    def length_mean_streamline(self, r_mid, d_int, d1):
        """Calculate length mean streamline.

        :param r_mid (float): radius mean streamline [m]
        :param d_int (float): center radius curvature front shroud [m]
        :param d1 (float): diameter [m]
        :return l_mid (float): length mean streamline [m]
        """
        l_mid = math.pi / 2 * r_mid + (d1 - d_int) / 2

        return l_mid

    def width_impeller_vane(self, r_mid, l_mid, d_int, a0, a1):
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

    def meridional_velocity_component_2(self, b2, d2, x2, flow, eta_vol):
        """Calculate meridional velocity component of the absolute velocity
        at section 2.

        :param b2 (float): impeller width [m]
        :param d2 (float): measured diameter
        :param x2 (float): blade blockage
        :param flow (float): flow rate [m^3/s]
        :param eta_vol (float): volumetric efficency
        :return cm2 (float): meridional velocity component [m/s]
        """
        cm2 = flow / (math.pi * d2 * b2 * x2 * eta_vol)

        return cm2

    def meridional_velocity_component_1(self, b1, d1, x1, flow, eta_vol):
        """Calculate meridional velocity component of the absolute velocity
        at section 1.

        :param b1 (float): impeller width [m]
        :param d1 (float): diameter [m]
        :param x1 (float): blade blockage
        :param flow (float): flow rate [m^3/s]
        :param eta_vol (float): volumetric efficency
        :return cm1 (float): meridional velocity component [m/s]
        """
        cm1 = flow / (math.pi * d1 * b1 * x1 * eta_vol)

        return cm1

    def circumferential_velocity_component_2(self, u2, cm2, beta_2c):
        """Calculate circumferential velocity component of the absolute
        velocity at section 2.

        :param u2 (float): circumferential velocity [m/s]
        :param cm2 (float): meridional velocity component [m/s]
        :param beta_2c (float): angle between rel. and circum. velocity [m/s]
        :return cu2 (float): circumferential velocity component [m/s]
        """
        cu2 = u2 - cm2 * 1 / math.tan(beta_2c)

        return cu2

    def circumferential_velocity_component_1(self, u1, cm1, beta_1c):
        """Calculate circumferential velocity component of the absolute
        velocity at section 1.

        :param u1 (float): circumferential velocity [m/s]
        :param cm1 (float): meridional velocity component [m/s]
        :param beta_1c (float): angle between rel. and circum. velocity [m/s]
        :return cu1 (float): circumferential velocity component [m/s]
        """
        cu1 = u1 - cm1 * 1 / math.tan(beta_1c)

        return cu1

    def circumferential_velocity_2(self, omega, d2):
        """Calculate circumferential velocity at section 2.

        :param omega (float): angular velocity [rad/s]
        :param d2 (float): measured diameter
        :return u2 (float): circumferential velocity [m/s]
        """
        u2 = omega * d2 / 2

        return u2

    def circumferential_velocity_1(self, omega, d1):
        """Calculate circumferential velocity at section 1.

        :param omega (float): angular velocity [rad/s]
        :param d1 (float): diameter [m]
        :return u1 (float): circumferential velocity [m/s]
        """
        u1 = omega * d1 / 2

        return u1

    def relative_velocity_2(self, cm2, beta_2c):
        """Calculate relative velocity at section 2.

        :param cm2 (float): meridional velocity [m/s]
        :param beta_2c (float): angle between rel. and circum. velocity [m/s]
        :return w2 (float): relative velocity [m/s]
        """
        w2 = cm2 / math.sin(beta_2c)

        return w2

    def relative_velocity_1(self, cm1, beta_1c):
        """Calculate relative velocity at section 1.

        :param cm1 (float): meridional velocity [m/s]
        :param beta_1c (float): angle between rel. and circum. velocity [m/s]
        :return w1 (float): relative velocity [m/s]
        """
        w1 = cm1 / math.sin(beta_1c)

        return w1

    def angle_beta_2c(self, cm2, u2, gamma_2):
        """Calculate blade working angle between relative and circumferential
        velocity at section 2.

        :param cm2 (float): meridional velocity [m/s]
        :param u2 (float): circumferential velocity [m/s]
        :param gamma_2 (int): measured angle between cm2 and vertical [deg]
        :return beta_2c (float): angle between rel. and circum. velocity [m/s]
        """
        gammar_2 = math.radians(gamma_2)
        beta_2c = math.atan((cm2 * math.cos(gammar_2)) / u2)

        return beta_2c

    def angle_beta_1c(self, psi_th, phi_th, u1_sf, u1):
        """Calculate blade working angle between relative and circumferential
        velocity at section 1.

        :param psi_th (float): theoretic head coefficient
        :param phi_th (float): flow coefficient corrected
        :param u1_sf (float): slip factor [m/s]
        :param u1 (float): circumferential velocity [m/s]
        :return beta_1c (float): angle between rel. and circum. velocity [m/s]
        """
        beta_1c = math.atan(phi_th / (1 - psi_th - u1_sf / u1))

        return beta_1c

    def head_coefficient(self, u1, head):
        """Calculate head coefficient at section 1.

        :param u1 (float): circumferential velocity [m/s]
        :param head (float): head [m]
        :return psi (float): head coefficient
        """
        psi = (G * head) / u1**2

        return psi

    def flow_coefficient(self, d1, b1, u1, x1, flow, eta_vol):
        """Calculate flow coefficient at section 1.

        :param d1 (float): diameter [m]
        :param b1 (float): impeller width [m]
        :param u1 (float): circumferential velocity [m/s]
        :param x1 (float): blade blockage
        :param flow (float): flow rate [m^3/s]
        :param eta_vol (float): volumetric efficency
        :return phi (float): flow coefficient
        """
        phi = flow / (math.pi * d1 * b1 * u1 * x1 * eta_vol)

        return phi

    def theoretic_head_coefficient(self, psi, eta_idr):
        """Calculate theoretic head coefficient at section 1.

        :param psi (float): head coefficient
        :param eta_idr (float): idraulic efficency
        :return psi_th (float): theoretic head coefficient
        """
        psi_th = psi / eta_idr

        return psi_th

    def theoretic_flow_coefficient(self, phi, x1, eta_vol):
        """Calculate theoretic flow coefficient at section 1.

        :param phi (float): flow coefficient
        :param x1 (float): blade blockage
        :param eta_vol (float): volumetric efficency
        :return phi_th (float): flow coefficient corrected
        """
        phi_th = phi / (x1 * eta_vol)

        return phi_th

    def slip_factor(self, u1, beta_1c, z):
        """Calculate slip factor at section 1 with Wiesner's formula.

        :param u1 (float): circumferential velocity [m/s]
        :param beta_1c (float): angle between rel. and circum. velocity [m/s]
        :param z (int): number of blades
        :return u1_sf (float): slip factor [m/s]
        """
        u1_sf = u1 * (math.sin(beta_1c))**0.5 / z**0.7

        return u1_sf

    def degree_reaction(self, phi_th, beta_1c, z):
        """Calculate degree of reaction.

        :param phi_th (float): flow coefficient corrected
        :param beta_1c (float): angle between rel. and circum. velocity [m/s]
        :param z (int): number of blades
        :return epsilon_r (float): degree of reaction
        """
        epsilon_r = 1 - .5 * (1 - phi_th / math.tan(beta_1c) -
                              (math.sin(beta_1c))**.5 / z**.7)

        return epsilon_r


class Volute(object):
    """Methods for impeller design."""

    def absolute_velocity_throat(self, cu1):
        """Calculate absolute velocity at throat section.

        :param cu1 (float): circumferential velocity component [m/s]
        :return c_thr (float): absolute velocity [m/s]
        """
        c_thr = .5 * cu1

        return c_thr

    def area_throat(self, flow, c_thr):
        """Calculate area at throat section.

        :param flow (float): flow rate [m^3/s]
        :param c_thr (float): absolute velocity [m/s]
        :return a_thr (float): area [m^2]
        """
        a_thr = flow / c_thr

        return a_thr

    def radius_start(self, d1):
        """Calculate internal radius at the start wrap angle.

        :param d1 (float): diameter [m]
        :return r3 (float): radius [m]
        """
        r1 = d1 / 2
        r3 = 1.1 * r1

        return r3

    def width_start(self, b1):
        """Calculate width at the start wrap angle.

        :param b1 (float): impeller width [m]
        :return b3 (float): width [m]
        """
        b3 = 1.715 * b1

        return b3

    def width_volute_vane(self, a_thr, theta_3):
        """Calculate width volute vane at different wrap angles.

        :param a_thr (float): area [m^2]
        :param theta_3 (int): start wrap angle [deg]
        :return b_vl (list): diameters at different wrap angles [m]
        """
        n = 8  # num. divisions circumference
        theta = []  # cumulative angles
        b_vl = []
        theta_3 = math.radians(theta_3)
        for i in range(n + 1):  # 0, .., n
            theta_i = ((2 * math.pi / n) * i)
            if theta_i > theta_3:
                theta.append(theta_i)
                b_vl.append((2 * a_thr * (theta_i - theta_3) /
                            math.pi**2)**.5)

        return list(zip(theta, b_vl))


class Project(Pre_Values, Shaft, Impeller, Volute):
    """Execute the project of a centrifugal pump."""

    def __init__(self, **kwargs):
        """Take input variables and execute the project.

        :param flow (float): flow rate [m^3/s]
        :param head (float) head [m]
        :param psi_coef (list): head coefficients
        :param phi_coef (list): head coefficients
        :param eta_coef (list): total efficency
        :param slip (int): slip for electric induction motor [%]
        :param hz (int):utility frequency [Hz]
        :param tau_adm (int): tau admissible [MPa]
        :param thk (float): blade thickness [m]
        :param lm (float): loss coefficient at section 0
        :param lw (float): low-pressure peak coefficient at blades at section 0
        :param km (float): rate between circumeferential velocity cm2 and c0
        :param eta_vol (float): volumetric efficency
        :param eta_idr (float): idraulic efficency
        :param d2 (float): measured diameter
        :param gamma_2 (int): measured angle between cm2 and vertical [deg]
        :param z (int): number of blades
        :param theta_3 (int): start wrap angle [deg]
        """
        flow = kwargs["flow"]
        head = kwargs["head"]
        slip = kwargs["slip"]
        hz = kwargs["hz"]
        tau_adm = kwargs["tau_adm"]
        thk = kwargs["thk"]
        lm = kwargs["lm"]
        lw = kwargs["lw"]
        km = kwargs["km"]
        eta_vol = kwargs["eta_vol"]
        eta_idr = kwargs["eta_idr"]
        d2 = kwargs["d2"]
        gamma_2 = kwargs["gamma_2"]
        z = kwargs["z"]
        theta_3 = kwargs["theta_3"]

        # feasability study
        fs_rpm = []
        fs_k_num = []
        fs_phi = []
        fs_psi = []
        fs_eta = []
        fs_u1 = []
        fs_d1 = []
        fs_b1 = []
        fs_bd1 = []
        fs_npsh_r = []
        for i, cp in enumerate(CPOLES):
            rpm = self.rotational_speed(slip, hz, cp)
            k_num = self.type_number(rpm, flow, head)
            # allow typical numbers only in the domain of centrifugal pumps
            if 0.2 <= k_num <= 1.2:
                fs_rpm.append(rpm)
                fs_k_num.append(k_num)
                fs_phi.append(self.flow_number(k_num))
                fs_psi.append(self.head_number(k_num))
                fs_eta.append(self.efficency(k_num))
                fs_u1.append(self._circumferential_velocity_1(head, fs_psi[i]))
                fs_d1.append(self._diameter_1(fs_u1[i], fs_rpm[i]))
                fs_b1.append(self._width_1(fs_u1[i], fs_d1[i], flow,
                                           fs_phi[i]))
                fs_bd1.append(self.width_over_diameter_1(fs_b1[i], fs_d1[i]))
                fs_npsh_r.append(self.npsh_r(k_num, head))

        # choose the best typical number
        # avoid k_num > .55 becuase it requires double curvature blades
        k_num = []
        for k in fs_k_num:
            if k <= .55:
                k_num.append(k)
        if len(k_num) > 0:
            k_num = max(k_num)

        # project values
        rpm = self.k2n(k_num, flow, head)
        cp = self.n2cp(rpm, slip, hz)
        phi = self.flow_number(k_num)
        psi = self.head_number(k_num)
        eta = self.efficency(k_num)
        u1 = self._circumferential_velocity_1(head, psi)
        d1 = self._diameter_1(u1, rpm)
        b1 = self._width_1(u1, d1, flow, phi)
        bd1 = self.width_over_diameter_1(b1, d1)
        npsh_r = self.npsh_r(k_num, head)

        # pump shaft design
        omega = self.angular_velocity(rpm)
        powr = self.power(eta, flow, head)
        torq = self.torque(powr, omega)
        d_sh = self.diameter_shaft(torq, tau_adm)
        d_hu = self.diameter_hub(d_sh)

        # impeller design
        dif = 1
        err = .001
        u1 = [u1]
        while dif > err:
            d1 = self.diameter_1(omega, u1[-1])
            d1 = round(d1 * 1000) / 1000
            u1.append(self.circumferential_velocity_1(omega, d1))
            dif = abs(u1[-1] - u1[-2])
        u1 = u1[-1]

        psi = self.head_coefficient(u1, head)

        dif = 1
        err = .001
        x0 = [1]
        while dif > err:
            d0_npsh = self.diameter_0_npsh(omega, x0[-1], flow, lm, lw, km,
                                           eta_vol)
            d0_eff = self.diameter_0_efficency(omega, x0[-1], flow, km,
                                               eta_vol)
            d0_cmp = self.diameter_0_compromise(omega, x0[-1], flow, eta_vol)
            d0 = self.diameter_0(d0_npsh, d0_eff, d0_cmp)
            x0.append(self.hub_blockage_0(d0, d_hu))
            if len(x0) > 2:
                dif = abs(x0[-1] - x0[-2])
        x0 = x0[-1]

        d_mid = self.diameter_mean_streamline(d_hu, d0)
        r_cvt = self.radius_curvature_front_shroud(d1)
        d_int = self.center_mean_streamline(r_cvt, d0)
        r_mid = self.radius_mean_streamline(d_int, d_mid)
        l_mid = self.length_mean_streamline(r_mid, d_int, d1)
        a0 = self.area_0(d_hu, d0)

        dif = 1
        err = .001
        x1 = [1]
        u1_sf = 0
        while dif > err:
            b1 = self.width_1(d1, u1, phi, x1[-1], flow, eta_vol)
            b1 = round(b1 * 1000) / 1000
            phi = self.flow_coefficient(d1, b1, u1, x1[-1], flow, eta_vol)
            a1 = self.area_1(d1, b1, x1[-1])
            b_im = self.width_impeller_vane(r_mid, l_mid, d_int, a0, a1)
            psi_th = self.theoretic_head_coefficient(psi, eta_idr)
            phi_th = self.theoretic_flow_coefficient(phi, x1[-1], eta_vol)
            beta_1c = self.angle_beta_1c(psi_th, phi_th, u1_sf, u1)
            epsilon_r = self.degree_reaction(phi_th, beta_1c, z)
            u1_sf = self.slip_factor(u1, beta_1c, z)
            x1.append(self.blade_blockage_1(beta_1c, d1, thk, z))
            cm1 = self.meridional_velocity_component_1(b1, d1, x1[-1], flow,
                                                       eta_vol)
            cu1 = self.circumferential_velocity_component_1(u1, cm1,  beta_1c)
            w1 = self.relative_velocity_1(cm1, beta_1c)
            if len(x1) > 2:
                dif = abs(x1[-1] - x1[-2])
        x1 = x1[-1]

        dif = 1
        err = .001
        x2 = [1]
        while dif > err:
            theta_2 = self.angle_theta_2(d_int, r_mid, d2)
            b2 = self.width_2(a0, a1, r_mid, l_mid, theta_2, d2)
            u2 = self.circumferential_velocity_2(omega, d2)
            cm2 = self.meridional_velocity_component_2(b2, d2, x2[-1], flow,
                                                       eta_vol)
            beta_2c = self.angle_beta_2c(cm2, u2, gamma_2)
            cu2 = self.circumferential_velocity_component_2(u2, cm2, beta_2c)
            w2 = self.relative_velocity_2(cm2, beta_2c)
            x2.append(self.blade_blockage_2(beta_2c, thk, z, d2))
            if len(x2) > 2:
                dif = abs(x2[-1] - x2[-2])
        x2 = x2[-1]

        # volute design
        r3 = self.radius_start(d1)
        b3 = self.width_start(b1)
        c_thr = self.absolute_velocity_throat(cu1)
        a_thr = self.area_throat(flow, c_thr)
        b_vl = self.width_volute_vane(a_thr, theta_3)

        self.results = locals()


def main(**kwargs):
    """Print results."""
    prj = Project(**kwargs)

    for k, v in sorted(list(prj.results.items())):
        if type(v) in (float, int, list):
            print(k, ' ', v)


if __name__ == '__main__':
    main(flow=.011, head=25,  # [m^3/s], [m]
         slip=3, hz=50,  # slip and utility frequency for electric motor
         tau_adm=30,  # tau admissible for C40 steel [MPa]
         thk=.003,  # blade thickness [m]
         lm=.04,  # loss coefficient at section 0
         lw=.50,  # low-pressure peak coefficient at blades at section 0
         km=1.2,  # rate between circumeferential velocity cm2 and c0
         eta_vol=.940,  # volumetric efficency
         eta_idr=.88,  # idraulic efficency
         d2=.112,  # measured diameter [m]
         gamma_2=5,  # measured angle between cm2 and vertical [deg]
         z=7,  # number of blades
         theta_3=10)  # start wrap angle for volute [deg]
