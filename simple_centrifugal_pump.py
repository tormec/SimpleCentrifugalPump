#!/usr/bin/env python3
"""Calculation of geometrical dimensions for centrifugal pump."""

import math

# name sections
# 0: impeller eye
# 1: impeller blade trailing edge
# 2: impeller blade leading edge

# physical constants
G = 9.81  # gravity acceleration [m/s^2]
RHO = 1000  # water density

# couple ploes for electric motor
CPOLES = [2, 4, 6, 8]

# internal diameters for commercial pipes
D_INT = [.054, .070, .082, .107, .132, .159, .207, .260, .310, .340]  # [m]


class Pre_Values(object):
    """Methods for feasability study of the impeller."""

    def rotational_speed(self, slip, hz):
        """Calculate rotational speed at diff. couple poles.

        :param slip (int): slip for electric induction motor [%]
        :param hz (int): utility frequency [Hz]
        :return rpm (list): rotational speed [rpm]
        """
        rpm = []
        for cp in CPOLES:
            rpm.append(120 * hz / cp * (1 - slip / 100))

        return rpm

    def type_number(self, rpm, flow, head):
        """Calculate centrifugal pump's type number at diff. couple poles.

        :param rpm (list): rotational speed [rpm]
        :param flow (float): flow rate [m^3/s]
        :param head (float): head [m]
        :return k_num (list): type number
        """
        k_num = []
        for n in rpm:
            omega = 2 * math.pi * n / 60
            k_num.append(omega * flow**0.5 / (G * head)**0.75)

        return k_num

    def _circumferential_velocity_1(self, head, psi_coef):
        """Calculate circumferential velocity at section 1
        at diff. couple poles.

        :param head (float): head [m]
        :param psi_coef (list): head coefficients
        :return u1 (list): circumferential velocity [m/s]
        """
        u1 = []
        for psi in psi_coef:
            u1.append((G * head / psi)**0.5)

        return u1

    def _diameter_1(self, u1, rpm):
        """Calculate diameter at section 1 at diff. couple poles.

        :param u1 (list): circumferential velocity [m/s]
        :param rpm (list): rotational speed [rpm]
        :return d1 (list): diameter [m]
        """
        d1 = []
        zipped = zip(u1, rpm)
        for u, n in zipped:
            d1.append((60 * u) / (math.pi * n))

        return d1

    def _width_1(self, u1, d1, flow, phi_coef):
        """Calculate impeller width at section 1 at diff. couple poles.

        :param u1 (list): circumferential velocity [m/s]
        :param d1 (list): diameter [m]
        :param flow (float): flow rate [m^3/s]
        :param phi_coef (list): flow coefficients
        :return b1 (list): impeller width [m]
        """
        b1 = []
        zipped = zip(u1, d1, phi_coef)
        for u, d, phi in zipped:
            b1.append(flow / (math.pi * d * u * phi))

        return b1

    def width_over_diameter_1(self, b1, d1):
        """Calculate rate width over diameter at section 1
        at diff. couple poles.

        :param b1 (list): impeller width [m]
        :param d1 (list): diameter [m]
        :return bd1 (list): width over diameter
        """
        bd1 = []
        zipped = zip(b1, d1)
        for b, d in zipped:
            bd1.append(b / d)

        return bd1

    def npsh_r(self, k_num, head):
        """Calculate the neat positive suction head required
        at diff. couple poles.

        :param k_num (list): type number
        :param head (float): head [m]
        :return npsh_r (list): neat positive suction head required [m]
        """
        npsh_r = []
        for k in k_num:
            npsh_r.append(.25 * k**(4/3) * head)

        return npsh_r


class Shaft(object):
    """Methods for pump shaft design."""

    def angular_velocity(self, rpm):
        """Calculate angular velocity pump shaft.

        :param rpm (float): rotational speed [rpm]
        :return omega (float): angular velocity [rad/s]
        """
        omega = 2 * math.pi * rpm / 60

        return omega

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
        :return r_theta_st (float): radius [m]
        """
        r1 = d1 / 2
        r_theta_st = 1.1 * r1

        return r_theta_st

    def width_start(self, b1):
        """Calculate width at the start wrap angle.

        :param b1 (float): impeller width [m]
        :return b_theta_0 (float): width [m]
        """
        b_theta_st = 1.715 * b1

        return b_theta_st

    def width_volute_vane(self, a_thr, theta_st):
        """Calculate width volute vane at different wrap angles.

        :param a_thr (float): area [m^2]
        :param theta_st (float): start wrap angle [rad]
        :return b_vl (list): diameters at different wrap angles [m]
        """
        n = 8  # num. divisions circumference
        theta = []  # cumulative angles
        b_vl = []
        theta.append(math.radians(theta_st))
        b_vl.append((2 * a_thr * theta[-1] / math.pi**2)**.5)
        for i in range(1, n + 1):
            theta.append((2 * math.pi / n) * i)
            b_vl.append((2 * a_thr * theta[-1] / math.pi**2)**.5)

        return list(zip(theta, b_vl))


class Test(Pre_Values, Shaft, Impeller, Volute):
    """Test methods."""

    def __init__(self, **kwargs):
        """Take input variables and execute test.

        :param flow (float): flow rate [m^3/s]
        :param head (float) head [m]
        :param psi_coef (list): head coefficients
        :param phi_coef (list): head coefficients
        :param eta_coef (list): total efficency
        :param slip (int): slip for electric induction motor [%]
        :param hz (int):utility frequency [Hz]
        :param cps (int): the chosen couple poles for electric motor
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
        """
        self.flow = kwargs["flow"]
        self.head = kwargs["head"]
        self.fs_psi = kwargs["psi_coef"]
        self.fs_phi = kwargs["phi_coef"]
        self.fs_eta = kwargs["eta_coef"]
        self.slip = kwargs["slip"]
        self.hz = kwargs["hz"]
        self.cps = kwargs["cps"]
        self.tau_adm = kwargs["tau_adm"]
        self.thk = kwargs["thk"]
        self.lm = kwargs["lm"]
        self.lw = kwargs["lw"]
        self.km = kwargs["km"]
        self.eta_vol = kwargs["eta_vol"]
        self.eta_idr = kwargs["eta_idr"]
        self.d2 = kwargs["d2"]
        self.gamma_2 = kwargs["gamma_2"]
        self.z = kwargs["z"]
        self.theta_st = kwargs["theta_st"]

        # feasability study
        self.fs_rpm = self.rotational_speed(self.slip, self.hz)
        self.fs_k_num = self.type_number(self.fs_rpm, self.flow, self.head)
        self.fs_u1 = self._circumferential_velocity_1(self.head, self.fs_psi)
        self.fs_d1 = self._diameter_1(self.fs_u1, self.fs_rpm)
        self.fs_b1 = self._width_1(self.fs_u1, self.fs_d1, self.flow,
                                   self.fs_phi)
        self.fs_bd1 = self.width_over_diameter_1(self.fs_b1, self.fs_d1)
        self.fs_npsh_r = self.npsh_r(self.fs_k_num, self.head)

        # project values
        val = CPOLES.index(self.cps)
        self.rpm = self.rotational_speed(self.slip, self.hz)[val]
        self.k_num = self.type_number(self.fs_rpm, self.flow, self.head)[val]
        self.phi = self.fs_phi[val]
        self.psi = self.fs_psi[val]
        self.eta_tot = self.fs_eta[val]
        self.u1 = self._circumferential_velocity_1(self.head, self.fs_psi)[val]
        self.d1 = self._diameter_1(self.fs_u1, self.fs_rpm)[val]
        self.b1 = self._width_1(self.fs_u1, self.fs_d1, self.flow,
                                self.fs_phi)[val]
        self.npsh_r = self.npsh_r(self.fs_npsh_r, self.head)[val]

        # pump shaft design
        self.omega = self.angular_velocity(self.rpm)
        self.powr = self.power(self.eta_tot, self.flow, self.head)
        self.torq = self.torque(self.powr, self.omega)
        self.d_sh = self.diameter_shaft(self.torq, self.tau_adm)
        self.d_hu = self.diameter_hub(self.d_sh)

        # impeller design
        dif = 1
        err = .001
        u1 = [self.u1]
        while dif > err:
            self.d1 = self.diameter_1(self.omega, u1[-1])
            self.d1 = round(self.d1 * 1000) / 1000
            u1.append(self.circumferential_velocity_1(self.omega, self.d1))
            dif = abs(u1[-1] - u1[-2])
        self.u1 = u1[-1]

        self.psi = self.head_coefficient(self.u1, self.head)

        dif = 1
        err = .001
        x0 = [1]
        while dif > err:
            self.d0_npsh = self.diameter_0_npsh(self.omega, x0[-1], self.flow,
                                                self.lm, self.lw, self.km,
                                                self.eta_vol)
            self.d0_eff = self.diameter_0_efficency(self.omega, x0[-1],
                                                    self.flow, self.km,
                                                    self.eta_vol)
            self.d0_cmp = self.diameter_0_compromise(self.omega, x0[-1],
                                                     self.flow, self.eta_vol)
            self.d0 = self.diameter_0(self.d0_npsh, self.d0_eff, self.d0_cmp)
            x0.append(self.hub_blockage_0(self.d0, self.d_hu))
            if len(x0) > 2:
                dif = abs(x0[-1] - x0[-2])
        self.x0 = x0[-1]

        self.d_mid = self.diameter_mean_streamline(self.d_hu, self.d0)
        self.r_cvt = self.radius_curvature_front_shroud(self.d1)
        self.d_int = self.center_mean_streamline(self.r_cvt, self.d0)
        self.r_mid = self.radius_mean_streamline(self.d_int, self.d_mid)
        self.l_mid = self.length_mean_streamline(self.r_mid, self.d_int,
                                                 self.d1)
        self.a0 = self.area_0(self.d_hu, self.d0)

        dif = 1
        err = .001
        x1 = [1]
        self.u1_sf = 0  # initialization
        while dif > err:
            self.b1 = self.width_1(self.d1, self.u1, self.phi, x1[-1],
                                   self.flow, self.eta_vol)
            self.b1 = round(self.b1 * 1000) / 1000
            self.phi = self.flow_coefficient(self.d1, self.b1, self.u1,
                                             x1[-1], self.flow, self.eta_vol)
            self.a1 = self.area_1(self.d1, self.b1, x1[-1])
            self.b_im = self.width_impeller_vane(self.r_mid, self.l_mid,
                                                 self.d_int, self.a0, self.a1)
            self.psi_th = self.theoretic_head_coefficient(self.psi,
                                                          self.eta_idr)
            self.phi_th = self.theoretic_flow_coefficient(self.phi, x1[-1],
                                                          self.eta_vol)
            self.beta_1c = self.angle_beta_1c(self.psi_th, self.phi_th,
                                              self.u1_sf, self.u1)
            self.epsilon_r = self.degree_reaction(self.phi_th, self.beta_1c,
                                                  self.z)
            self.u1_sf = self.slip_factor(self.u1, self.beta_1c, self.z)
            x1.append(self.blade_blockage_1(self.beta_1c, self.d1, self.thk,
                                            self.z))
            self.cm1 = self.meridional_velocity_component_1(self.b1, self.d1,
                                                            x1[-1], self.flow,
                                                            self.eta_vol)
            self.cu1 = self.circumferential_velocity_component_1(self.u1,
                                                                 self.cm1,
                                                                 self.beta_1c)
            self.w1 = self.relative_velocity_1(self.cm1, self.beta_1c)
            if len(x1) > 2:
                dif = abs(x1[-1] - x1[-2])
        self.x1 = x1[-1]

        dif = 1
        err = .001
        x2 = [1]
        while dif > err:
            self.theta_2 = self.angle_theta_2(self.d_int, self.r_mid, self.d2)
            self.b2 = self.width_2(self.a0, self.a1, self.r_mid, self.l_mid,
                                   self.theta_2, self.d2)
            self.u2 = self.circumferential_velocity_2(self.omega, self.d2)
            self.cm2 = self.meridional_velocity_component_2(self.b2, self.d2,
                                                            x2[-1], self.flow,
                                                            self.eta_vol)
            self.beta_2c = self.angle_beta_2c(self.cm2, self.u2, self.gamma_2)
            self.cu2 = self.circumferential_velocity_component_2(self.u2,
                                                                 self.cm2,
                                                                 self.beta_2c)
            self.w2 = self.relative_velocity_2(self.cm2, self.beta_2c)
            x2.append(self.blade_blockage_2(self.beta_2c, self.thk, self.z,
                                            self.d2))
            if len(x2) > 2:
                dif = abs(x2[-1] - x2[-2])
        self.x2 = x2[-1]

        # volute design
        self.r_theta_st = self.radius_start(self.d1)
        self.b_theta_st = self.width_start(self.b1)
        self.c_thr = self.absolute_velocity_throat(self.cu1)
        self.a_thr = self.area_throat(self.flow, self.c_thr)
        self.b_vl = self.width_volute_vane(self.a_thr, self.theta_st)


def main(**kwargs):
    """Execute test and print results."""
    test = Test(**kwargs)

    for k, v in sorted(list(globals().items())):
        if type(v) in (float, int, list):
            print(k, ' ', v)

    l = test.__dict__.items()  # list of tuples [(k, v)]
    for i in sorted(l, key=lambda l: l[0]):  # sort by k
        if type(i[1]) is list:
            t = []
            for j in i[1]:
                if type(j) is tuple:
                    t.append(tuple(round(k, 5) for k in j))
                elif type(j) is float:
                    t.append(round(j, 5))
            print(i[0], ' ', t)
        else:
            print(i[0], ' ', round(i[1], 5))


if __name__ == '__main__':
    main(flow=.018, head=28,
         slip=3, hz=50,  # slip and utility frequency for electric motor
         psi_coef=[.50, .55, .69],  # head coef. for diff. couple poles
         phi_coef=[.13, .09, .08],  # flow coef. for diff. couple poles
         eta_coef=[.88, .76, .73],  # efficency for diff. couple poles
         cps=4,  # the chosen couple poles for electric motor
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
         theta_st=10)  # start wrap angle for volute [deg]