#!/usr/bin/env python3
"""Calculation of geometrical dimensions for centrifugal pump."""

import math

# name sections
# 0: impeller eye
# 1: impeller blade trailing edge
# 2: impeller blade leading edge

# given data
FLOW = 0.018  # flow rate [m^3/s]
HEAD = 28  # head [m]

# physical constants
G = 9.81  # gravity acceleration [m/s^2]
RHO = 1000  # water density

# caracteristics of the electric motor
CPOLES = [2, 4, 6]  # number of couple poles
SLIP = 3  # slip for electric induction motor [%]
HZ = 50  # frequency [Hz]

# steel chosen for the pump shaft
TAU = 30  # tau admissible for C40 steel [MPa]

# blades
S = .003  # width [m]

# coefficients due to loss
LM = 0.04  # loss coefficient at section 0
LW = 0.50  # low-pressure peak coefficient at blades at section 0
KM = 1.2  # relation between circumeferential velocity cm2 and c0
ETA_V = 0.940  # volumetric efficency
ETA_IDR = .87  # idraulic efficency

# internal diameters for commercial pipes
D_INT = [.054, .070, .082, .107, .132, .159, .207, .260, .310, .340]  # [m]

# TO CHANGE
FHT_2 = [.13, .50, .88]  # phi, psi, eta for el. motor with 2 couple poles
FHT_4 = [.09, .55, .76]  # phi, psi, eta for el. motor with 4 couple poles
FHT_6 = [.08, .69, .73]  # phi, psi, eta for el. motor with 6 couple poles
CPS = 4  # chosen el. motor with 4 couple poles
D2 = .112  # measured diameter [m]
GAMMA = 5  # measured angle between cm2 and the vertical at blade edge [deg]
Z = 7  # number of blades


class Pre_Values(object):
    """Feasibility study for choosing the principal dimensions."""

    def __init__(self):
        self.fs_rpm = self.rotational_speed()
        self.fs_k_num = self.type_number(self.fs_rpm)
        self.fs_u1 = Pre_Values.circumferential_velocity_1(self)
        self.fs_d1 = Pre_Values.diameter_1(self, self.fs_u1, self.fs_rpm)
        self.fs_b1 = Pre_Values.width_1(self, self.fs_u1, self.fs_d1)
        self.fs_bd1 = self.width_over_diameter_1(self.fs_b1, self.fs_d1)
        self.fs_npsh_r = self.npsh_r(self.fs_k_num)

    def rotational_speed(self):
        """Calculate rotational speed at diff. couple poles.

        :return rpm (list): rotational speed [rpm]
        """
        rpm = []
        for cp in CPOLES:
            rpm.append(120 * HZ / cp * (1 - SLIP / 100))

        return rpm

    def type_number(self, rpm):
        """Calculate centrifugal pump's type number at diff. couple poles.

        :param rpm (list): rotational speed [rpm]
        :return k_num (list): type number
        """
        k_num = []
        for n in rpm:
            omega = 2 * math.pi * n / 60
            k_num.append(omega * FLOW**0.5 / (G * HEAD)**0.75)

        return k_num

    def circumferential_velocity_1(self):
        """Calculate circumferential velocity at section 1
        at diff. couple poles.

        :return u1 (list): circumferential velocity [m/s]
        """
        u1 = []
        zipped = list(zip(FHT_2, FHT_4, FHT_6))
        head_coef = zipped[1]
        for psi in head_coef:
            u1.append((G * HEAD / psi)**0.5)

        return u1

    def diameter_1(self, u1, rpm):
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

    def width_1(self, u1, d1):
        """Calculate impeller width at section 1 at diff. couple poles.

        :param u1 (list): circumferential velocity [m/s]
        :param d1 (list): diameter [m]
        :return b1 (list): impeller width [m]
        """
        b1 = []
        zipped = list(zip(FHT_2, FHT_4, FHT_6))
        flow_coef = zipped[0]
        zipped = zip(u1, d1, flow_coef)
        for u, d, phi in zipped:
            b1.append(FLOW / (math.pi * d * u * phi))

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

    def npsh_r(self, k_num):
        """Calculate the neat positive suction head required
        at diff. couple poles.

        :param k_num (list): type number
        :return npsh_r (list): neat positive suction head required [m]
        """
        npsh_r = []
        for k in k_num:
            npsh_r.append(.25 * k**(4/3) * HEAD)

        return npsh_r


class Project_Values(Pre_Values):
    """Project values chosen from the feasibility study."""

    def __init__(self):
        Pre_Values.__init__(self)  # initialize base class's init()

        val = CPOLES.index(CPS)
        self.rpm = self.rotational_speed()[val]
        self.k_num = self.type_number(self.fs_rpm)[val]
        self.phi = list(zip(FHT_2, FHT_4, FHT_6))[0][val]
        self.psi = list(zip(FHT_2, FHT_4, FHT_6))[1][val]
        self.eta_tot = list(zip(FHT_2, FHT_4, FHT_6))[2][val]
        self.u1 = Pre_Values.circumferential_velocity_1(self)[val]
        self.d1 = Pre_Values.diameter_1(self, self.fs_u1, self.fs_rpm)[val]
        self.b1 = Pre_Values.width_1(self, self.fs_u1, self.fs_d1)[val]
        self.npsh_r = self.npsh_r(self.fs_npsh_r)[val]


class Shaft(Project_Values):
    """Dimensioning pump shaft."""

    def __init__(self):
        Project_Values.__init__(self)  # initialize base class's init()

        self.omega = self.angular_velocity(self.rpm)
        self.powr = self.power(self.eta_tot)
        self.torq = self.torque(self.powr, self.omega)
        self.d_sh = self.diameter_shaft(self.torq)
        self.d_hu = self.diameter_hub(self.d_sh)

    def angular_velocity(self, rpm):
        """Calculate angular velocity pump shaft.

        :param rpm (float): rotational speed [rpm]
        :return omega (float): angular velocity [rad/s]
        """
        omega = 2 * math.pi * rpm / 60

        return omega

    def power(self, eta_tot):
        """Calculate power the el. motor should provide at pump shaft.

        :param eta_t (float): total efficency
        :return powr (float): power [W]
        """
        powr = (RHO * FLOW * G * HEAD) / eta_tot

        return powr

    def torque(self, powr, omega):
        """Calculate torque at pump shaft.

        :param powr (float): power [W]
        :param omega (float): angular velocity [rad/s]
        :return torq (float): torque [Nm]
        """
        torq = powr / omega

        return torq

    def diameter_shaft(self, torq):
        """Calculate shaft diameter.

        :param torq (float): torque [Nm]
        return d_sh (float): shaft diameter [m]
        """
        d_sh = ((32 * torq) / (math.pi * (TAU * 10**6)))**(1/3)

        return d_sh

    def diameter_hub(self, d_sh):
        """Calculate hub diameter.

        :param d_sh (float): shaft diameter [m]
        :return d_hu (float): hub diameter [m]
        """
        d_hu = 1.5 * d_sh

        return d_hu


class Impeller(Shaft):
    """Dimensioning impeller."""

    def __init__(self):
        Shaft.__init__(self)  # initialize base class's init()

        x0 = [1]
        x1 = [1]
        x2 = [1]
        u1_sf = [0]
        df_x0 = 1
        df_x1 = 1
        df_x2 = 1
        df_u1_sf = 1
        er = .001
        while (df_x0 > er and df_x1 > er and df_x2 > er and df_u1_sf > er):
            # calculate diameter at section 1
            self.d1 = self.diameter_1(self.omega, self.u1)
            # round diameter at section 1
            self.d1 = round(self.d1 * 1000) / 1000
            # calculate circumferential velocity at section 1
            self.u1 = self.circumferential_velocity_1(self.omega, self.d1)
            # calculate head coefficient
            self.psi = self.head_coefficient(self.u1)
            # calculate  impeller width at section 1
            self.b1 = self.width_1(self.d1, self.u1, self.phi, x1[-1])
            # round impeller width at section 1
            self.b1 = round(self.b1 * 1000) / 1000
            # calculate flow coefficient
            self.phi = self.flow_coefficient(self.d1, self.b1, self.u1)
            # calculate diameter at section 0 with min npsh_r
            self.d0_npsh = self.diameter_0_npsh(self.omega, x0[-1])
            # calculate diameter at section 0 with max tot efficency
            self.d0_eff = self.diameter_0_efficency(self.omega, x0[-1])
            # calculate diameter at section 0 as compromise solution
            self.d0_cmp = self.diameter_0_compromise(self.omega, x0[-1])
            # calculate diameter at section 0 as average of diameters
            self.d0 = self.diameter_0(self.d0_npsh, self.d0_eff, self.d0_cmp)
            # calculate hub blockage
            x0.append(self.hub_blockage_0(self.d0, self.d_hu))
            # calculate diameter mean streamline
            self.d_mid = self.diameter_mean_streamline(self.d_hu, self.d0)
            # calculate radius curvature front shroud
            self.r_cvt = self.radius_curvature_front_shroud(self.d1)
            # calculate center radius curvature front shroud
            self.d_int = self.center_mean_streamline(self.r_cvt, self.d0)
            # calculate radius mean streamline
            self.r_mid = self.radius_mean_streamline(self.d_int, self.d_mid)
            # calculate length mean streamline
            self.l_mid = self.length_mean_streamline(self.r_mid, self.d_int,
                                                     self.d1)
            # calculate area at section 0
            self.a0 = self.area_0(self.d_hu, self.d0)
            # calculate area at section 1
            self.a1 = self.area_1(self.d1, self.b1, x1[-1])
            # calculate diameters at equals angles in the curved zone
            self.b_vn = self.width_vane(self.r_mid, self.l_mid, self.d_int,
                                        self.a0, self.a1)
            # calculate angle between radius m. stream. and vert. at sec. 2
            self.theta_2 = self.angle_theta_2(self.d_int, self.r_mid)
            # calculate impeller width at section 2
            self.b2 = self.width_2(self.a0, self.a1, self.r_mid, self.l_mid,
                                   self.theta_2)
            # calculate circumferential velocity at section 2
            self.u2 = self.circumferential_velocity_2(self.omega)
            # calculate meridional velocity at section 2
            self.cm2 = self.meridional_velocity_2(self.b2, x2[-1])
            # calculate blade working angle at section 2
            self.beta_2c = self.angle_beta_2c(self.cm2, self.u2)
            # calculate blade blockage at section 2
            x2.append(self.blade_blockage_2(self.beta_2c))
            # calculate relative velocity at section 2
            self.w2 = self.relative_velocity_2(self.cm2, self.beta_2c)
            # calculate theoretic head coefficient
            self.psi_th = self.theoretic_head_coefficient(self.psi)
            # calculate theoretic flow coefficient
            self.phi_th = self.theoretic_flow_coefficient(self.phi, x1[-1])
            # calculate blade working angle at section 1
            self.beta_1c = self.angle_beta_1c(self.psi_th, self.phi_th,
                                              u1_sf[-1], self.u1)
            # calculate slip factor at section 1
            u1_sf.append(self.slip_factor(self.u1, self.beta_1c))
            # calculate blade blockage at section 1
            x1.append(self.blade_blockage_1(self.beta_1c, self.d1))
            # calculate meridional velocity at section 1
            self.cm1 = self.meridional_velocity_1(self.b1, self.d1, x1[-1])
            # calculate relative velocity at section 1
            self.w1 = self.relative_velocity_1(self.cm1, self.beta_1c)
            # calculate degree of reaction
            self.epsilon_r = self.degree_reaction(self.phi_th, self.beta_1c)

            if len(x0) > 1:
                    df_x0 = x0[-1] - x0[-2]
            if len(x1) > 1:
                    df_x1 = x1[-1] - x1[-2]
            if len(x2) > 1:
                    df_x2 = x2[-1] - x2[-2]
            if len(u1_sf) > 1:
                    df_u1_sf = u1_sf[-1] - u1_sf[-2]

        self.x0 = x0[-1]
        self.x1 = x1[-1]
        self.x2 = x2[-1]
        self.u1_sf = u1_sf[-1]

        print(math.degrees(self.beta_2c))
        print(math.degrees(self.beta_1c))

    def hub_blockage_0(self, d0, d_hu):
        """Calculate hub blockage at section 0.

        :param d0 (float): diameter [m]
        :param d_hu (float): hub diameter [m]
        :return x0 (float): hub blockage
        """
        x0 = 1 - (d_hu / d0)**2

        return x0

    def blade_blockage_1(self, beta_1c, d1):
        """Calculate blade blockage at section 1.

        :param beta_1c (float): angle between rel. and circum. velocity [m/s]
        :param d1 (float): diameter
        :return x1 (float): blade blockage
        """
        x1 = 1 - (Z * S) / (math.pi * d1 * math.sin(beta_1c))

        return x1

    def blade_blockage_2(self, beta_2c):
        """Calculate blade blockage at section 2.

        :param beta_2c (float): angle between rel. and circum. velocity [m/s]
        :return x2 (float): blade blockage
        """
        x2 = 1 - (Z * S) / (math.pi * D2 * math.sin(beta_2c))

        return x2

    def diameter_1(self, omega, u1):
        """Calculate diameter at section 1.

        :param omega (float): angular velocity [rad/s]
        :param u1 (float): circumferential velocity [m/s]
        :return d1 (float): diameter [m]
        """
        d1 = 2 * u1 / omega

        return d1

    def diameter_0_npsh(self, omega, x0):
        """Calculate diameter at section 0 with min npsh_r.

        :param omega (float): angular velocity [rad/s]
        :param x0 (float): hub blockage
        :return d0_npsh (float): diameter with min npsh_r [m]
        """
        d0_npsh = 2 * (
                       (2 * FLOW**2 * KM**2 * (1 + LM + LW)) /
                       (ETA_V**2 * math.pi**2 * omega**2 * x0**2 * LW)
                      )**(1/6)

        return d0_npsh

    def diameter_0_efficency(self, omega, x0):
        """Calculate diameter at section 0 with max total efficency.

        :param omega (float): angular velocity [rad/s]
        :param x0 (float): hub blockage
        :return d0_eff (float): diameter with max efficency [m]
        """
        d0_eff = 2 * (
                      (2 * FLOW**2 * KM**2) /
                      (ETA_V**2 * math.pi**2 * omega**2 * x0**2)
                     )**(1/6)

        return d0_eff

    def diameter_0_compromise(self, omega, x0):
        """Calculate diameter at section 0 as compromise solution between
        min npsh_r and max total efficency.

        :param omega (float): angular velocity [rad/s]
        :param x0 (float): hub blockage
        :return d0_cmp (float): diameter as compromise solution [m]
        """
        d0_cmp = (
                 (FLOW * 8 * 3.03) /
                 (omega * math.pi * x0 * ETA_V)
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

    def angle_theta_2(self, d_int, r_mid):
        """Calculate angle between radius mean streamline and vertical
        at section 2.

        :param d_int (float): center radius curvature front shroud [m]
        :param r_mid (float): radius mean streamline [m]
        :return theta_2 (float): angle between radius m. stream. and vert.[rad]
        """
        theta_2 = math.acos((d_int - D2) / (2 * r_mid))

        return theta_2

    def width_2(self, a0, a1, r_mid, l_mid, t2):
        """Calculate impeller width at section 2.

        :param a0 (float): area [m^2]
        :param a1 (float): area [m^2]
        :param r_mid (float): radius mean streamline [m]
        :param l_mid (float): length mean streamline [m]
        :param theta_2 (float): angle between radius m. stream. and vert. [rad]
        :return b2 (float): impeller width [m]
        """
        b2 = (a0 + (a1 - a0) * (r_mid * t2) / l_mid) / (math.pi * D2)

        return b2

    def width_1(self, d1, u1, phi, x1):
        """Calculate impeller width at section 1.

        :param d1 (float): diameter [m]
        :param u1 (float): circumferential velocity [m/s]
        :param phi (float): flow coefficient
        :x1 (float): blade blockage
        """
        b1 = FLOW / (math.pi * d1 * u1 * phi * x1 * ETA_V)

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

    def width_vane(self, r_mid, l_mid, d_int, a0, a1):
        """Calculate width impeller vane as different diameters
        at equals angles in the curved zone.

        :param r_mid (float): radius mean streamline [m]
        :param l_mid (float): length mean streamline [m]
        :param d_int (float): center radius curvature front shroud [m]
        :param a0 (float): area [m^2]
        :param a1 (float): area [m^2]
        :return b_vn (list): diameters at equals angles in the curved zone [m]
        """
        # in the curved zone, the mean streamline is a quarter of arc of circle
        # and along its path are considered:
        n = 11  # num. divisions: points = segments + 1
        theta = []  # cumulative angles
        length = []  # cumulative lengths
        area = []  # area at different lengths
        r = []  # radius at different angles
        b_vn = []
        for i in range(n):
            theta.append((math.pi / (2 * n)) * i)
            length.append(r_mid * theta[i])
            area.append(a0 + (a1 - a0) * length[i] / l_mid)
            r.append(d_int / 2 - r_mid * math.cos(theta[i]))
            b_vn.append(area[i] / (2 * math.pi * r[i]))

        return list(zip(theta, length, b_vn, area))

    def meridional_velocity_2(self, b2, x2):
        """Calculate meridional velocity at section 2.

        :param b2 (float): impeller width [m]
        :param x2 (float): blade blockage
        :return cm2 (float): meridional velocity [m/s]
        """
        cm2 = FLOW / (math.pi * D2 * b2 * x2 * ETA_V)

        return cm2

    def meridional_velocity_1(self, b1, d1, x1):
        """Calculate meridional velocity at section 1.

        :param b1 (float): impeller width [m]
        :param d1 (float): diameter [m]
        :param x1 (float): blade blockage
        :return cm1 (float): meridional velocity [m/s]
        """
        cm1 = FLOW / (math.pi * d1 * b1 * x1 * ETA_V)

        return cm1

    def circumferential_velocity_2(self, omega):
        """Calculate circumferential velocity at section 2.

        :param omega (float): angular velocity [rad/s]
        :return u2 (float): circumferential velocity [m/s]
        """
        u2 = omega * D2 / 2

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

    def angle_beta_2c(self, cm2, u2):
        """Calculate blade working angle between relative and circumferential
        velocity at section 2.

        :param cm2 (float): meridional velocity [m/s]
        :param u2 (float): circumferential velocity [m/s]
        :return beta_2c (float): angle between rel. and circum. velocity [m/s]
        """
        gammar = math.radians(GAMMA)
        beta_2c = math.atan((cm2 * math.cos(gammar)) / u2)

        return beta_2c

    def angle_beta_1c(self, psi_th, phi_c, u1_sf, u1):
        """Calculate blade working angle between relative and circumferential
        velocity at section 1.

        :param psi_th (float): theoretic head coefficient
        :param phi_c (float): flow coefficient corrected
        :param u1_sf (float): slip factor [m/s]
        :param u1 (float): circumferential velocity [m/s]
        :return beta_1c (float): angle between rel. and circum. velocity [m/s]
        """
        beta_1c = math.atan(phi_c / (1 - psi_th - u1_sf / u1))

        return beta_1c

    def head_coefficient(self, u1):
        """Calculate head coefficient at section 1.

        : param u1 (float): circumferential velocity [m/s]
        :return psi (float): head coefficient
        """
        psi = (G * HEAD) / u1**2

        return psi

    def flow_coefficient(self, d1, b1, u1):
        """Calculate flow coefficient at section 1.

        :param d1 (float): diameter [m]
        :param b1 (float): impeller width [m]
        :param u1 (float): circumferential velocity [m/s]
        :return phi (float): flow coefficient
        """
        phi = FLOW / (math.pi * d1 * b1 * u1)

        return phi

    def theoretic_head_coefficient(self, psi):
        """Calculate theoretic head coefficient at section 1.

        :param psi (float): head coefficient
        :return psi_th (float): theoretic head coefficient
        """
        psi_th = psi / ETA_IDR

        return psi_th

    def theoretic_flow_coefficient(self, phi, x1):
        """Calculate theoretic flow coefficient at section 1.

        :param phi (float): flow coefficient
        :param x1 (float): blade blockage
        :return phi_c (float): flow coefficient corrected
        """
        phi_c = phi / (x1 * ETA_V)

        return phi_c

    def slip_factor(self, u1, beta_1c):
        """Calculate slip factor at section 1 with Wiesner's formula.

        :param u1 (float): circumferential velocity [m/s]
        :param beta_1c (float): angle between rel. and circum. velocity [m/s]
        :return u1_sf (float): slip factor [m/s]
        """
        u1_sf = u1 * (math.sin(beta_1c))**0.5 / Z**0.7

        return u1_sf

    def degree_reaction(self, phi_c, beta_1c):
        """Calculate degree of reaction.

        :param phi_c (float): flow coefficient corrected
        :param eta_1c (float): angle between rel. and circum. velocity [m/s]
        :return epsilon_r (float): degree of reaction
        """
        epsilon_r = 1 - 0.5 * (1 - phi_c / math.tan(beta_1c) -
                               (math.sin(beta_1c))**0.5 / Z**0.7)

        return epsilon_r


def main():
    mp = Impeller()

    for k, v in sorted(list(globals().items())):
        if type(v) in (float, int, list):
            print(k, ' ', v)

    l = mp.__dict__.items()  # list of tuples [(k, v)]
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
    main()
