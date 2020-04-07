import math
import constants as CNST


class Impeller(object):
    """Methods to calculate the impeller."""

    def hub_blockage(self, d, d_hu):
        """Calculate hub blockage.

        :param d (float): diameter [m]
        :param d_hu (float): hub diameter [m]
        :return x (float): hub blockage
        """
        x = 1 - (d_hu / d)**2

        return x

    def blade_blockage_1(self, beta_c, d, thk, z):
        """Calculate blade blockage.

        :param beta_c (float): angle between rel. and blade velocity [m/s]
        :param d (float): diameter [m]
        :param thk (float): blade thickness [m]
        :param z (int): number of blades
        :return x (float): blade blockage
        """
        x = 1 - (z * thk) / (math.pi * d * math.sin(beta_c))

        return x

    def diameter_omega(self, omega, u):
        """Calculate diameter function of angular velocity.

        :param omega (float): angular velocity [rad/s]
        :param u (float): absolute velocity [m/s]
        :return d (float): diameter [m]
        """
        d = 2 * u / omega

        return d

    def diameter_npsh(self, omega, x, flow, lm, lw, km, eta_vol):
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

    def diameter_efficency(self, omega, x, flow, km, eta_vol):
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

    def diameter_flow(self, omega, x, flow, eta_vol):
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

    def diameter_avg(self, d_npsh, d_eff, d_flow):
        """Calculate diameter as average value.

        :param d_npsh (float): diameter with min NPSH_r [m]
        :param d_eff (float): diameter with max efficency [m]
        :param d_flow (float): diameter as function of the flow rate [m]
        :return d_avg (float): diameter as average value [m]
        """
        d_vals = [d_npsh, d_eff, d_flow]
        d_avg = math.fsum(d_vals) / len(d_vals)

        return d_avg

    def diameter_std(self, d):
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

    def streamline_diam(self, d_hu, d_0):
        """Calculate middle streamline diameter at section 0.

        :param d_hu (float): hub diameter [m]
        :param d (float): diameter [m]
        :return d_mid (float): middle streamline diameter [m]
        """
        d_mid = (d_0 + d_hu) / 2

        return d_mid

    def curvature_rad_shroud(self, d_1):
        """Calculate curvature radius of the shroud at section 0.

        :param d1 (float): diameter [m]
        :return r_cvt (float): curvature radius [m]
        """
        r_cvt = .06 * d_1

        return r_cvt

    def streamline_len(self, r_cvt, d_hu, d_1, d_0):
        """Calculate length mean streamline.

        :param r_cvt (float): curvature radius [m]
        :param d_hub (float): hub diameter [m]
        :param d_1 (float): diameter at section 1 [m]
        :param d_0 (float): diameter at section 0 [m]
        :return l_mid (float): middle streamline length [m]
        """
        l_mid = math.pi / 2 * ((d_0 - d_hu) / 4 + r_cvt) +\
            (d_1 - d_0 - 2 * r_cvt) / 2

        return l_mid

    def angle_theta_2(self, r_cvt, d_hu, d_1, d_0, d_2):
        """Calculate angle between radius middle streamline and vertical axis
        at section 2.

        :param r_cvt (float): curvature radius [m]
        :param d_hub (float): hub diameter [m]
        :param d_1 (float): diameter at section 1 [m]
        :param d_0 (float): diameter at section 0 [m]
        :param d_2 (float): diameter at section 2 [m]
        :return theta_2 (float): angle between streamline rad. and vert. [rad]
        """
        theta_2 = math.acos((d_1 - d_0 - 2 * r_cvt - d_2) /
                            (2 * ((d_0 - d_hu) / 4 + r_cvt)))

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

    def width_1(self, d_1, u_1, phi, x_1, flow, eta_vol):
        """Calculate impeller width at section 1.

        :param d_1 (float): diameter at section 1 [m]
        :param u_1 (float): absolute velocity at section 1 [m/s]
        :param phi (float): flow coefficient
        :param flow (float): flow rate [m^3/s]
        :param eta_vol (float): volumetric efficency
        :return b_1 (float): impeller width at section 1
        """
        b1 = flow / (math.pi * d_1 * u_1 * phi * x_1 * eta_vol)

        return b1

    def area_0(self, d_hu, d_0):
        """Calculate area at section 0.

        :param d_hu (float): hub diameter [m]
        :param d_0 (float): diameter at section 0 [m]
        :return a0 (float): area [m^2]
        """
        a_0 = (d_0**2 - d_hu**2) * math.pi / 4

        return a_0

    def area_1(self, d1, b1, x1):
        """Calculate area at section 1.

        :param d1 (float): diameter [m]
        :param b1 (list): impeller width [m]
        :param x1 (float): blade blockage
        :return a1 (float): area [m^2]
        """
        a1 = math.pi * d1 * b1 * x1

        return a1

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

    def meridional_abs_vel(self, b, d, x, flow, eta_vol):
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

    def circumeferential_abs_vel(self, u, c_m, beta_c):
        """Calculate circumferential velocity component of the absolute
        velocity.

        :param u (float): absolute velocity [m/s]
        :param c_m (float): meridional velocity component [m/s]
        :param beta_c (float): angle between rel. and blade velocity [m/s]
        :return c_u (float): absolute velocity component [m/s]
        """
        c_u = u - c_m * 1 / math.tan(beta_c)

        return c_u

    def blade_vel(self, omega, d):
        """Calculate blade velocity function of angular velocity.

        :param omega (float): angular velocity [rad/s]
        :param d (float): diameter [m]
        :return u (float): absolute velocity [m/s]
        """
        u = omega * d / 2

        return u

    def relative_vel(self, c_m, beta_c):
        """Calculate relative velocity.

        :param c_m (float): meridional velocity [m/s]
        :param beta_c (float): angle between rel. and blade velocity [m/s]
        :return w (float): relative velocity [m/s]
        """
        w = c_m / math.sin(beta_c)

        return w

    def angle_beta_2c(self, cm2, u2, gamma_2):
        """Calculate blade working angle between relative and peripheral
        velocity at section 2.

        :param cm2 (float): meridional velocity [m/s]
        :param u2 (float): absolute velocity [m/s]
        :param gamma_2 (int): measured angle between cm2 and vertical [deg]
        :return beta_2c (float): angle between rel. and circum. velocity [m/s]
        """
        gammar_2 = math.radians(gamma_2)
        beta_2c = math.atan((cm2 * math.cos(gammar_2)) / u2)

        return beta_2c

    def angle_beta_1c(self, psi_th, phi_th, u_1sf, u_1):
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

    def head_number(self, u, head):
        """Calculate head number.

        :param u (float): absolute velocity [m/s]
        :param head (float): head [m]
        :return psi (float): head coefficient
        """
        psi = (CNST.G * head) / u**2

        return psi

    def flow_number(self, d, b, u, x, flow, eta_vol):
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

    def theoretic_head_number(self, psi, eta_hyd):
        """Calculate theoretic head coefficient.

        :param psi (float): head coefficient
        :param eta_hyd (float): hydraulic efficency
        :return psi_th (float): theoretic head coefficient
        """
        psi_th = psi / eta_hyd

        return psi_th

    def theoretic_flow_number(self, phi, x, eta_vol):
        """Calculate theoretic flow coefficient.

        :param phi (float): flow coefficient
        :param x (float): blade blockage
        :param eta_vol (float): volumetric efficency
        :return phi_th (float): flow coefficient corrected
        """
        phi_th = phi / (x * eta_vol)

        return phi_th

    def slip_factor(self, u, beta_c, z):
        """Calculate slip factor with Wiesner"s formula.

        :param u (float): absolute velocity [m/s]
        :param beta_c (float): angle between rel. and blade velocity [m/s]
        :param z (int): number of blades
        :return u_sf (float): slip factor [m/s]
        """
        u_sf = u * (math.sin(beta_c))**0.5 / z**0.7

        return u_sf

    def degree_reaction(self, phi_th, beta_c, z):
        """Calculate degree of reaction.

        :param phi_th (float): flow coefficient corrected
        :param beta_c (float): angle between rel. and blade velocity [m/s]
        :param z (int): number of blades
        :return epsilon_ract (float): degree of reaction
        """
        epsilon_ract = 1 - (1 - phi_th / math.tan(beta_c) -
                            (math.sin(beta_c))**.5 / z**.7) / 2

        return epsilon_ract
