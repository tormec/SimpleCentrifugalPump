#!/usr/bin/env python3
"""Calculation of geometrical dimensions for centrifugal pump."""

import calc
import constants as CN
import shaft as sh
import impeller as im
import volute as vl


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
        self.km = 1.2  # rate between peripheral velocity cm_1 and c_0
        self.z = 7  # number of blades

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
        u_2 = []
        d_2 = []
        b_2 = []
        bd_2 = []
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
                u_2.append(im.psi2u(psi[i], self.head))
                d_2.append(im.diameter_omega(im.rpm2omega(n), u_2[i]))
                b_2.append(im.width(d_2[i], None, u_2[i], phi[i], self.flow))
                bd_2.append(im.width0diameter(b_2[i], d_2[i]))
                npsh_req.append(im.cappa2npsh(k, self.head))

        results = {}
        for i in ["part", "pp", "rpm", "cappa", "phi", "psi", "eta", "u_2",
                  "d_2", "b_2", "bd_2", "npsh_req"]:
            results[i] = locals()[i]

        return results

    def chose_option(self, **kwargs):
        """Select an option of design according to a criteria."""
        cappa = kwargs["cappa"]

        part = "---chosen option---"
        for k in cappa:
            if k > .55:
                # avoid it, because it requires double curvature blades
                cappa.remove(k)
        if len(cappa) > 0:
            cappa = max(cappa)
        rpm = im.cappa2rpm(cappa, self.flow, self.head)
        pp = im.rpm2pp(rpm, self.slip, self.hz)
        phi = im.flow_number_poly(cappa)
        psi = im.head_number_poly(cappa)
        eta = im.efficency_poly(cappa)
        u_2 = im.psi2u(psi, self.head)
        d_2 = im.diameter_omega(im.rpm2omega(rpm), u_2)
        b_2 = im.width(d_2, None, u_2, phi, self.flow)
        bd_2 = im.width0diameter(b_2, d_2)
        npsh_req = im.cappa2npsh(cappa, self.head)

        results = {}
        for i in ["part", "pp", "rpm", "cappa", "phi", "psi", "eta", "u_2",
                  "d_2", "b_2", "bd_2", "npsh_req"]:
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
        cappa = kwargs["cappa"]
        phi = kwargs["phi"]
        u_2 = kwargs["u_2"]
        omega = kwargs["omega"]
        d_hu = kwargs["d_hu"]

        part = "---impeller---"
        eta_hyd = im.efficency_hyd_poly(cappa)
        eta_vol = im.efficency_vol_poly(cappa)
        u_2 = [u_2]
        dif = 1
        err = .001
        while dif > err:
            d_2 = im.diameter_omega(omega, u_2[-1])
            d_2 = round(d_2, 3)
            u_2.append(im.blade_vel(omega, d_2))
            dif = abs(u_2[-1] - u_2[-2])
        u_2 = u_2[-1]

        psi = im.head_number(u_2, self.head)
        psi_th = im.theoretic_head_number(psi, eta_hyd)

        x_0 = [1]
        dif = 1
        err = .001
        while dif > err:
            d_0npsh = im.diameter_npsh(omega, x_0[-1], self.flow, self.lm,
                                       self.lw, self.km, eta_vol)
            d_0eff = im.diameter_efficency(omega, x_0[-1], self.flow, self.km,
                                           eta_vol)
            d_0flow = im.diameter_flow(omega, x_0[-1], self.flow, eta_vol)
            d_0avg = im.average_diam(d_0npsh, d_0eff, d_0flow)
            d_0 = im.standard_diam(d_0avg)
            x_0.append(im.hub_blockage(d_0, d_hu))
            dif = abs(x_0[-1] - x_0[-2])
        x_0 = x_0[-1]

        d_sl = im.streamline_diam(d_hu, d_0)
        r_c = im.curvature_rad(d_2)
        r_slc = im.streamline_curv_rad(d_hu, d_0, r_c)
        l_sl = im.streamline_len(r_slc, d_2, d_sl)

        x_2 = [1]
        u_2sf = 0
        dif = 1
        err = .001
        while dif > err:
            b_2 = im.width(d_2, None, u_2, phi, self.flow, x_2[-1],
                           eta_vol)
            b_2 = round(b_2, 3)
            phi = im.flow_number(d_2, b_2, u_2, x_2[-1], self.flow,
                                 eta_vol)
            phi_th = im.theoretic_flow_number(phi, x_2[-1], eta_vol)
            c_2m = im.meridional_abs_vel(u_2, phi_th)
            beta_2 = im.angle_beta(u_2, c_2m, 0, psi_th, u_2sf)
            epsilon_ract = im.degree_reaction(phi_th, beta_2, self.z)
            u_2sf = im.slip_factor(u_2, beta_2, self.z)
            x_2.append(im.blade_blockage(beta_2, d_2, self.thk, self.z))
            c_2u = im.circumferential_abs_vel(u_2, c_2m,  beta_2)
            w_2 = im.relative_vel(c_2m, beta_2)
            dif = abs(x_2[-1] - x_2[-2])
        x_2 = x_2[-1]

        theta = []
        b = []
        x = []
        t = 0
        n = 16
        for i in range(n):
            x_i = [1]
            dif = 1
            err = .001
            theta_i = im.angle_theta(n, i)
            d_isl = im.streamline_diam(d_hu, d_0, theta_i, r_slc)
            l_isl = im.streamline_len(r_slc, theta=theta_i)
            gamma_i = im.angle_gamma(r_c, r_slc, theta_i)
            if gamma_i is None:
                theta.append(theta_i)
                a_i = im.area(l_isl, l_sl, d_hu, d_0, d_2, b_2, x_2)
                b.append(im.width(d_isl, a_i))
                x.append(0)
                t = i
                continue
            u_i = im.blade_vel(omega, d_isl)
            while dif > err:
                a_i = im.area(l_isl, l_sl, d_hu, d_0, d_2, b_2, x_2)
                b_i = im.width(d_isl, a_i)
                phi_i = im.flow_number(d_isl, b_i, u_i, x_i[-1], self.flow,
                                       eta_vol)
                phi_ith = im.theoretic_flow_number(phi_i, x_i[-1], eta_vol)
                c_im = im.meridional_abs_vel(u_i, phi_ith)
                beta_i = im.angle_beta(u_i, c_im, gamma_i)
                w_i = im.relative_vel(c_im, beta_i)
                x_i.append(im.blade_blockage(beta_i, d_isl, self.thk, self.z))
                dif = abs(x_i[-1] - x_i[-2])
            theta.append(theta_i)
            b.append(b_i)
            x.append(x_i[-1])

        theta_1 = theta[t + 1]
        b_1 = b[t + 1]
        x_1 = x[t + 1]
        d_1sl = im.streamline_diam(d_hu, d_0, theta_1, r_slc)
        gamma_1 = im.angle_gamma(r_c, r_slc, theta_1)
        u_1 = im.blade_vel(omega, d_1sl)
        phi_1 = im.flow_number(d_1sl, b_1, u_1, x_1, self.flow, eta_vol)
        phi_1th = im.theoretic_flow_number(phi_1, x_1, eta_vol)
        c_1m = im.meridional_abs_vel(u_1, phi_1th)
        beta_1 = im.angle_beta(u_1, c_1m, gamma_1)
        w_1 = im.relative_vel(c_1m, beta_1)

        npsh_req = im.npsh_req(c_1m, w_1, self.lm, self.lw)

        results = {}
        for i in ["part",
                  "d_0npsh", "d_0eff", "d_0flow", "d_0avg", "d_0", "x_0",
                  "d_sl", "r_c", "r_slc", "l_sl",
                  "theta_1", "gamma_1", "beta_1", "b_1", "d_1sl", "x_1",
                  "u_1", "c_1m", "w_1",
                  "d_2", "u_2", "b_2",  "beta_2", "theta", "b",
                  "c_2m", "c_2u", "w_2", "u_2sf", "x_2",
                  "psi", "psi_th", "phi",  "phi_th", "epsilon_ract",
                  "npsh_req"]:
            if i in ["beta_2", "theta", "theta_1", "gamma_1", "beta_1"]:
                results[i] = calc.rad2deg(locals()[i])
            else:
                results[i] = locals()[i]

        return results

    def calc_volute(self, **kwargs):
        """Calculate the volute."""
        d_2 = kwargs["d_2"]
        b_2 = kwargs["b_2"]
        c_2u = kwargs["c_2u"]

        part = "---volute---"
        d_3 = vl.diameter(d_2)
        c_thr = vl.absolute_velocity_throat(c_2u)
        a_thr = vl.area(self.flow, c_thr)
        b_3, theta_3 = vl.width_min(b_2, a_thr)

        n = 9
        theta = []
        b = []
        for i in range(n):
            dim = vl.width(vl.angle_theta(n, i), a_thr, b_3)
            if dim is not None:
                b.append(dim[0])
                theta.append(dim[1])

        results = {}
        for i in ["part", "d_3", "b_3", "theta_3", "c_thr", "a_thr",
                  "theta", "b"]:
            if i in ["theta_3", "theta"]:
                results[i] = calc.rad2deg(locals()[i])
            else:
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
