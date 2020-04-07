#!/usr/bin/env python3
"""Calculation of geometrical dimensions for centrifugal pump."""

import calc
import constants as CNST
import options as opt
import shaft as shf
import impeller as imp
import volute as vlt

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
        :param d2 (float): measured diameter
        :param gamma_2 (int): measured angle between cm2 and vertical [deg]
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
        self.d2 = .112  # measured diameter [m]
        self.gamma_2 = 5  # measured angle between cm2 and vertical [deg]
        self.z = 7  # number of blades
        self.theta_3 = 10  # start wrap angle for volute [deg]

        options = self.calc_options()
        choice = self.chose_option(**options)
        shaft = self.calc_shaft(**choice)
        impeller = self.calc_impeller(**{**choice, **shaft})
        volute = self.calc_volute(**impeller)
        self.results = [options, choice, shaft, impeller, volute]

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
        npsh_r = []
        for i, p in enumerate(CNST.PPAIRS):
            n = opt.rotational_speed(p, self.slip, self.hz)
            k = opt.type_number(n, self.flow, self.head)
            # only typical numbers in the domain of centrifugal pumps
            if 0.2 <= k <= 1.2:
                pp.append(p)
                rpm.append(n)
                cappa.append(k)
                phi.append(opt.flow_number(k))
                psi.append(opt.head_number(k))
                eta.append(opt.efficency(k))
                u_1.append(opt.peripheral_velocity(psi[i], self.head))
                d_1.append(opt.diameter(u_1[i], rpm[i]))
                b_1.append(opt.width(u_1[i], d_1[i], phi[i], self.flow))
                bd_1.append(opt.width0diameter(b_1[i], d_1[i]))
                npsh_r.append(opt.npsh_r(k, self.head))

        results = {}
        for i in ["part", "pp", "rpm", "cappa", "phi", "psi", "eta", "u_1",
                  "d_1", "b_1", "bd_1", "npsh_r"]:
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
        rpm = opt.cappa2rpm(cappa, self.flow, self.head)
        pp = opt.rpm2pp(rpm, self.slip, self.hz)
        phi = opt.flow_number(cappa)
        psi = opt.head_number(cappa)
        eta = opt.efficency(cappa)
        u_1 = opt.peripheral_velocity(psi, self.head)
        d_1 = opt.diameter(u_1, rpm)
        b_1 = opt.width(u_1, d_1, phi, self.flow)
        bd_1 = opt.width0diameter(b_1, d_1)
        npsh_r = opt.npsh_r(cappa, self.head)

        results = {}
        for i in ["part", "pp", "rpm", "cappa", "phi", "psi", "eta", "u_1",
                  "d_1", "b_1", "bd_1", "npsh_r"]:
            results[i] = locals()[i]

        return results

    def calc_shaft(self, **kwargs):
        """Calculate the pump shaft. """
        rpm = kwargs["rpm"]
        eta = kwargs["eta"]

        part = "---pump shaft---"
        omega = shf.angular_velocity(rpm)
        power = shf.power(eta, self.flow, self.head)
        torque = shf.torque(power, omega)
        d_sh = shf.shaft_diameter(torque, self.tau_adm)
        d_hu = calc.bisect(lambda d_hu, d_sh=d_sh:
                           (d_hu**4 - d_sh**4) / d_hu - d_sh**3,
                           d_sh - 1, d_sh + 1, .001)

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
        dif = 1
        err = .001
        u_1 = [u_1]
        while dif > err:
            d_1 = shf.diameter_omega(omega, u_1[-1])
            d_1 = round(d_1 * 1000) / 1000
            u_1.append(shf.blade_vel(omega, d_1))
            dif = abs(u_1[-1] - u_1[-2])
        u_1 = u_1[-1]

        psi = shf.head_number(u_1, self.head)

        dif = 1
        err = .001
        x_0 = [1]
        while dif > err:
            d_0_npsh = shf.diameter_npsh(omega, x_0[-1], self.flow,
                                         self.lm, self.lw, self.km,
                                         self.eta_vol)
            d_0_eff = shf.diameter_efficency(omega, x_0[-1], self.flow,
                                             self.km, self.eta_vol)
            d_0_flow = shf.diameter_flow(omega, x_0[-1], self.flow,
                                         self.eta_vol)
            d_0_avg = shf.average_diam(d_0_npsh, d_0_eff, d_0_flow)
            d_0 = shf.standard_diam(d_0_avg)
            x_0.append(self.hub_blockage(d_0, d_hu))
            if len(x_0) > 2:
                dif = abs(x_0[-1] - x_0[-2])
        x_0 = x_0[-1]

        d_mid = shf.streamline_diam(d_hu, d_0)
        r_cvt = shf.curvature_rad(d_1)
        r_mid = shf.streamline_rad(d_hu, d_0, r_cvt)
        l_mid = shf.streamline_len(r_cvt, r_mid, d_1, d_0)
        a_0 = shf.area_0(d_hu, d_0)

        dif = 1
        err = .001
        x_1 = [1]
        u_1_sf = 0
        while dif > err:
            b_1 = shf.width_1(d_1, u_1, phi, x_1[-1], self.flow, self.eta_vol)
            b_1 = round(b_1 * 1000) / 1000
            phi = shf.flow_coefficient(d_1, b_1, u_1, x_1[-1], self.flow,
                                       self.eta_vol)
            a_1 = shf.area_1(d_1, b_1, x_1[-1])
            psi_th = shf.theoretic_head_number(psi, self.eta_hyd)
            phi_th = shf.theoretic_flow_number(phi, x_1[-1], self.eta_vol)
            beta_1_c = shf.angle_beta_1c(psi_th, phi_th, u_1_sf, u_1)
            epsilon_ract = shf.degree_reaction(phi_th, beta_1_c, self.z)
            u_1_sf = shf.slip_factor(u_1, beta_1_c, self.z)
            x_1.append(shf.blade_blockage(beta_1_c, d_1, self.thk, self.z))
            c_1_m = shf.meridional_abs_vel(b_1, d_1, x_1[-1], self.flow,
                                           self.eta_vol)
            c_1_u = shf.circumferential_abs_vel(u_1, c_1_m,  beta_1_c)
            w_1 = shf.relative_velocity_1(c_1_m, beta_1_c)
            if len(x_1) > 2:
                dif = abs(x_1[-1] - x_1[-2])
        x_1 = x_1[-1]

        dif = 1
        err = .001
        x_2 = [1]
        while dif > err:
            theta_2 = shf.angle_theta_2(r_cvt, r_mid, d_1, d_0, d_2)
            b_2 = shf.width_2(a_0, a_1, r_mid, l_mid, theta_2, d_2)
            u_2 = shf.blade_vel(omega, d_2)
            c_2_m = shf.meridional_abs_vel(b_2, d_2, x_2[-1], self.flow,
                                           self.eta_vol)
            beta_2_c = shf.angle_beta_2c(c_2_m, u_2, self.gamma_2)
            c_2_u = shf.circumeferential_abs_vel(u_2, c_2_m, beta_2_c)
            w_2 = shf.relative_vel(c_2_m, beta_2_c)
            x_2.append(shf.blade_blockage(beta_2_c, d_2, self.thk, self.z))
            if len(x_2) > 2:
                dif = abs(x_2[-1] - x_2[-2])
        x_2 = x_2[-1]

        results = {}
        for i in ["part", "d_1", "u_1", "psi",
                  "d_0_npsh", "d_0_eff", "d_0_flow", "d_0_avg", "d_0", "x_0",
                  "d_mid", "r_cvt", "l_mid", "a_0",
                  "b_1", "phi", "a_1", "psi_th", "phi_th", "beta_1_c",
                  "epsilon_ract", "u_1_sf", "x_1", "c_1_m", "c_1_u", "w_1",
                  "theta_2", "b_2", "u_2", "c_2_m", "beta_2_c", "c_2_u", "w_2",
                  "x_2"]:
            results[i] = locals()[i]

        return results

    def calc_volute(self, **kwargs):
        """Calculate the volute."""
        d1 = kwargs["d1"]
        b1 = kwargs["b1"]
        cu1 = kwargs["cu1"]

        part = "---volute---"
        r3 = self.radius_start(d1)
        b3 = self.width_start(b1)
        c_thr = self.absolute_velocity_throat(cu1)
        a_thr = self.area_throat(self.flow, c_thr)
        b_vl = self.width_volute_vane(a_thr, self.theta_3)

        results = {}
        for i in ["part", "r3", "b3", "c_thr", "a_thr", "b_vl"]:
            results[i] = locals()[i]

        return results


def main(**kwargs):
    """Print results."""
    prj = Project(**kwargs)

    for result in prj.results:
        for key, val in result.items():
            if type(val) == list:
                for k, v in enumerate(val):
                    if type(v) == tuple:
                        val[k] = tuple(round(t, 3) for t in v)
                    else:
                        val[k] = round(v, 3)
                print(key, " ", val)
            elif type(val) in (float, int):
                print(key, " ", round(val, 3))
            else:
                print(val)


if __name__ == "__main__":
    main(flow=.011,  # [m^3/s]
         head=25)  # [m]
