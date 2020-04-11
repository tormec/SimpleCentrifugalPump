#!/usr/bin/env python3
"""Calculation of geometrical dimensions for centrifugal pump."""

import calc
import constants as CN
import shaft as sh
import impeller as im
import volute as vl

# name sections
# 0: impeller eye
# 1: impeller blade trailing edge
# 2: impeller blade leading edge
# 3: volute start wrap angle


class Project(object):
    """Execute the project of a centrifugal pump."""

    def __init__(self, **kwargs):
        """Take input variables and execute the project.

        :param flow (float): flow rate [m^3/s]
        :param head (float) head [m]
        :param slip (int): slip of electric induction motor [%]
        :param hz (int):utility frequency [Hz]
        :param tau_adm (int): tau admissible [MPa]
        :param thk (float): blade thickness [m]
        :param lm (float): loss coefficient at section 0
        :param lw (float): low-pressure peak coefficient at blades at section 0
        :param km (float): rate between circumeferential velocity cm2 and c0
        :param eta_vol (float): volumetric efficency
        :param eta_hyd (float): idraulic efficency
        :param z (int): number of blades
        :param theta_3 (int): start wrap angle [deg]
        """
        self.flow = kwargs["flow"]
        self.head = kwargs["head"]
        # suggested values
        self.slip = 3  # slip factor for AC motor
        self.hz = 50  # utility frequency for AC motor
        self.tau_adm = 30  # admissible stress for C40 steel [MPa]
        self.thk = .003  # blade thickness [m]
        self.lm = .04  # loss coefficient at section 0
        self.lw = .50  # low-pressure peak coefficient at blades at section 0
        self.km = 1.2  # rate between peripheral velocity cm2 and c0
        self.eta_vol = .940  # volumetric efficency
        self.eta_hyd = .880  # hydraulic efficency
        self.z = 7  # number of blades
        self.theta_3 = 10  # start wrap angle for volute [deg]

        options = self.calc_options()
        choice = self.chose_option(**options)
        shaft = self.calc_shaft(**choice)
        impeller = self.calc_impeller(**{**choice, **shaft})
        # volute = self.calc_volute(**impeller)
        self.results = [options, choice, shaft, impeller] #, volute]

    def calc_options(self):
        """Calculate several design options of an impeller."""
        part = "---options---"
        pp = []
        rpm = []
        cappa = []
        phi = []
        psi = []
        eta = []
        u_1 = []
        d_1 = []
        b_1 = []
        bd_1 = []
        npsh_req = []
        for i, p in enumerate(CN.PPAIRS):
            n = sh.rotational_speed(p, self.slip, self.hz)
            k = im.type_number(sh.angular_velocity(n), self.flow, self.head)
            # only typical numbers in the domain of centrifugal pumps
            if 0.2 <= k <= 1.2:
                pp.append(p)
                rpm.append(n)
                cappa.append(k)
                phi.append(im.flow_number_poly(k))
                psi.append(im.head_number_poly(k))
                eta.append(im.efficency_poly(k))
                u_1.append(im.psi2u(psi[i], self.head))
                d_1.append(im.diameter_omega(im.rpm2omega(n), u_1[i]))
                b_1.append(im.width(d_1[i], None, u_1[i], phi[i], self.flow))
                bd_1.append(im.width0diameter(b_1[i], d_1[i]))
                npsh_req.append(im.npsh_req(k, self.head))

        results = {}
        for i in ["part", "pp", "rpm", "cappa", "phi", "psi", "eta", "u_1",
                  "d_1", "b_1", "bd_1", "npsh_req"]:
            results[i] = locals()[i]

        return results

    def chose_option(self, **kwargs):
        """Select an option of design according to a criteria."""
        cappa = kwargs["cappa"]

        part = "---chosen option---"
        for k in cappa:
            if k > .55:
                # avoid, because it requires double curvature blades
                cappa.remove(k)
        if len(cappa) > 0:
            cappa = max(cappa)
        rpm = im.cappa2rpm(cappa, self.flow, self.head)
        pp = im.rpm2pp(rpm, self.slip, self.hz)
        phi = im.flow_number_poly(cappa)
        psi = im.head_number_poly(cappa)
        eta = im.efficency_poly(cappa)
        u_1 = im.psi2u(psi, self.head)
        d_1 = im.diameter_omega(im.rpm2omega(rpm), u_1)
        b_1 = im.width(d_1, None, u_1, phi, self.flow)
        bd_1 = im.width0diameter(b_1, d_1)
        npsh_req = im.npsh_req(cappa, self.head)

        results = {}
        for i in ["part", "pp", "rpm", "cappa", "phi", "psi", "eta", "u_1",
                  "d_1", "b_1", "bd_1", "npsh_req"]:
            results[i] = locals()[i]

        return results

    def calc_shaft(self, **kwargs):
        """Calculate the pump shaft. """
        rpm = kwargs["rpm"]
        eta = kwargs["eta"]

        part = "---pump shaft---"
        omega = sh.angular_velocity(rpm)
        power = sh.power(eta, self.flow, self.head)
        torque = sh.torque(power, omega)
        d_sh = sh.shaft_diameter(torque, self.tau_adm)
        d_sh = round(d_sh, 3)
        d_hu = calc.bisect(lambda d_hu, d_sh=d_sh:
                           (d_hu**4 - d_sh**4) / d_hu - d_sh**3,
                           d_sh - 1, d_sh + 1, .001)
        d_hu = round(d_hu, 3)

        results = {}
        for i in ["part", "omega", "power", "torque", "d_sh", "d_hu"]:
            results[i] = locals()[i]

        return results

    def calc_impeller(self, **kwargs):
        """Calculate the impeller."""
        phi = kwargs["phi"]
        u_1 = kwargs["u_1"]
        omega = kwargs["omega"]
        d_hu = kwargs["d_hu"]

        part = "---impeller---"
        u_1 = [u_1]
        dif = 1
        err = .001
        while dif > err:
            d_1 = im.diameter_omega(omega, u_1[-1])
            d_1 = round(d_1, 3)
            u_1.append(im.blade_vel(omega, d_1))
            dif = abs(u_1[-1] - u_1[-2])
        u_1 = u_1[-1]

        psi = im.head_number(u_1, self.head)

        x_0 = [1]
        dif = 1
        err = .001
        while dif > err:
            d_0npsh = im.diameter_npsh(omega, x_0[-1], self.flow, self.lm,
                                       self.lw, self.km, self.eta_vol)
            d_0eff = im.diameter_efficency(omega, x_0[-1], self.flow, self.km,
                                           self.eta_vol)
            d_0flow = im.diameter_flow(omega, x_0[-1], self.flow, self.eta_vol)
            d_0avg = im.average_diam(d_0npsh, d_0eff, d_0flow)
            d_0 = im.standard_diam(d_0avg)
            x_0.append(im.hub_blockage(d_0, d_hu))
            dif = abs(x_0[-1] - x_0[-2])
        x_0 = x_0[-1]

        d_sl = im.streamline_diam(d_hu, d_0)
        r_c = im.curvature_rad(d_1)
        r_slc = im.streamline_curv_rad(d_hu, d_0, r_c)
        l_sl = im.streamline_len(r_slc, d_1, d_sl)

        psi = im.head_number(u_1, self.head)
        psi_th = im.theoretic_head_number(psi, self.eta_hyd)

        x_1 = 1
        u_1sf = 0
        dif = 1
        err = .001
        c_1u = [im.psi_th2c_u(psi_th, u_1)]
        while dif > err:
            b_1 = im.width(d_1, None, u_1, phi, self.flow, x_1, self.eta_vol)
            b_1 = round(b_1, 3)
            phi = im.flow_number(d_1, b_1, u_1, x_1, self.flow, self.eta_vol)
            phi_th = im.theoretic_flow_number(phi, x_1, self.eta_vol)
            c_1m = im.meridional_abs_vel(b_1, d_1, x_1, self.flow,
                                         self.eta_vol)
            beta_1 = im.angle_beta(c_1m, u_1, 0, c_1u[-1], u_1sf)
            epsilon_ract = im.degree_reaction(phi_th, beta_1, self.z)
            u_1sf = im.slip_factor(u_1, beta_1, self.z)
            x_1 = im.blade_blockage(beta_1, d_1, self.thk, self.z)
            c_1u.append(im.circumferential_abs_vel(u_1, c_1m,  beta_1))
            w_1 = im.relative_vel(c_1m, beta_1)
            dif = abs(c_1u[-1] - c_1u[-2])
        c_1u = c_1u[-1]

        n = 11
        theta_i = []
        l_isl = []
        b_i = []
        c_im = []
        for i in range(n):
            x_i = [1]
            dif = 1
            err = .001
            theta_i.append(im.angle_theta(n, i))
            d_isl = im.streamline_diam(d_hu, d_0, theta_i[-1], r_slc)
            l_isl.append(im.streamline_len(r_slc, theta=theta_i[-1]))
            gamma_i = im.angle_gamma(r_c, r_slc, theta_i[-1])
            if gamma_i is None:
                continue
            u_i = im.blade_vel(omega, d_isl)
            phi_i = phi
            while dif > err:
                a_i = im.area(l_isl[-1], l_sl, d_hu, d_0, d_1, b_1, x_1)
                b = im.width(d_isl, a_i)
                phi_i = im.flow_number(d_isl, b, u_i, x_i[-1], self.flow,
                                       self.eta_vol)
                c_m = im.meridional_abs_vel(b, d_isl, x_i[-1], self.flow,
                                            self.eta_vol)
                beta_i = im.angle_beta(c_m, u_i, gamma_i)
                x_i.append(im.blade_blockage(beta_i, d_isl, self.thk, self.z))
                dif = abs(x_i[-1] - x_i[-2])
            b_i.append(b)
            c_im.append(c_m)

        # n = 6
        # for i in range(n):
        #     x_i = [1]
        #     #c_iu = c_iu
        #     dif = 1
        #     err = .001
        #     d_isl = im.streamline_diam(d_hu, d_0)
        #     l_isl.append(im.streamline_len(r_slc, d_1, d_isl))
        #     u_i = im.blade_vel(omega, d_isl)
        #     #phi_i = phi_i
        #     psi_i = im.head_number(u_i, self.head)
        #     psi_ith = im.theoretic_head_number(psi_i, self.eta_hyd)
        #     while dif > err:
        #         b = im.width(d_isl, u_i, phi_i, self.flow, x_i[-1],
        #                      self.eta_vol)
        #         phi_i = im.flow_number(d_isl, b, u_i, x_i[-1], self.flow,
        #                                self.eta_vol)
        #         c_m = im.meridional_abs_vel(b, d_isl, x_i[-1], self.flow,
        #                                     self.eta_vol)
        #         beta_i = im.angle_beta(c_m, u_i, 0, c_iu)
        #         x_i.append(im.blade_blockage(beta_i, d_isl, self.thk, self.z))
        #         c_iu = im.circumferential_abs_vel(u_i, c_m, beta_i)
        #         dif = abs(x_i[-1] - x_i[-2])
        #     b_i.append(b)
        #     c_im.append(c_m)

        # dif = 1
        # err = .001
        # x_2 = [1]
        # while dif > err:
        #     theta_2 = im.diameter2theta(r_cvt, r_mid, d_1, d_0, d_2)
        #     b_2 = im.width_2(a_0, a_1, r_mid, l_mid, theta_2, d_2)
        #     u_2 = im.blade_vel(omega, d_2)
        #     c_2_m = im.meridional_abs_vel(b_2, d_2, x_2[-1], self.flow,
        #                                    self.eta_vol)
        #     beta_2_c = im.angle_beta_2c(c_2_m, u_2, self.gamma_2)
        #     c_2_u = im.circumeferential_abs_vel(u_2, c_2_m, beta_2_c)
        #     w_2 = im.relative_vel(c_2_m, beta_2_c)
        #     x_2.append(im.blade_blockage(beta_2_c, d_2, self.thk, self.z))
        #     if len(x_2) > 2:
        #         dif = abs(x_2[-1] - x_2[-2])
        # x_2 = x_2[-1]

        results = {}
        for i in ["part", "d_1", "u_1", "psi",
                  "d_0npsh", "d_0eff", "d_0flow", "d_0avg", "d_0", "x_0",
                  "d_sl", "r_c", "r_slc", "l_sl",
                  "b_1", "phi", "psi_th", "phi_th", "beta_1",
                  "epsilon_ract", "u_1sf", "x_1", "c_1m", "c_1u", "w_1",
                  "theta_i", "l_isl", "b_i", "c_im"]:
            results[i] = locals()[i]

        return results

    # def calc_volute(self, **kwargs):
    #     """Calculate the volute."""
    #     d_1 = kwargs["d_1"]
    #     b_1 = kwargs["b_1"]
    #     c_1_u = kwargs["c_1_u"]

    #     part = "---volute---"
    #     r_3 = vl.radius_start(d_1)
    #     b_3 = vl.width_start(b_1)
    #     c_thr = vl.absolute_velocity_throat(c_1_u)
    #     a_thr = vl.area_throat(self.flow, c_thr)
    #     b_vl = vl.width_volute_vane(a_thr, self.theta_3)

    #     results = {}
    #     for i in ["part", "r_3", "b_3", "c_thr", "a_thr", "b_vl"]:
    #         results[i] = locals()[i]

    #     return results


def main(**kwargs):
    """Print results."""
    prj = Project(**kwargs)

    for result in prj.results:
        for key, val in result.items():
            if type(val) == list:
                for k, v in enumerate(val):
                    if type(v) == tuple:
                        val[k] = tuple(round(t, 6) for t in v)
                    else:
                        val[k] = round(v, 6)
                print(key, " ", val)
            elif type(val) in (float, int):
                print(key, " ", round(val, 6))
            else:
                print(val)


if __name__ == "__main__":
    main(flow=.011,  # [m^3/s]
         head=25)  # [m]
