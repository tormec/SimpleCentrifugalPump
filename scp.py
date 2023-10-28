#!/usr/bin/env python3
"""Calculation of geometrical dimensions for centrifugal pump."""

import argparse
import lib.constants as CN
import lib.options as op
import lib.calc as cl
import lib.shaft as sh
import lib.impeller as im
import lib.volute as vl


class Project(object):
    """Execute the project of a centrifugal pump."""

    def __init__(self, **kwargs):
        """Take input variables and execute the project.

        :param flow (float): flow rate [m^3/s]
        :param head (float) head [m]
        :param hz (int): frequency of alternating current [Hz]
        :param t (float): blade thickness [m]
        """
        self.flow = kwargs["flow"]
        self.head = kwargs["head"]
        self.hz = kwargs["hz"]
        self.t = kwargs["t"]
        # suggested values
        self.slip = 3  # slip factor for AC motor [%]
        self.tau_adm = 30  # shearing stress in pump shaft [MPa]
        self.lm = .04  # loss coefficient at section 0
        self.lw = .50  # low-pressure peak coefficient at blades at section 0
        self.km = 1.2  # rate between circumeferential velocity cm2 and c0
        self.z = [6, 7, 8]  # range number of blades
        self.beta_b = [cl.deg2rad(15), cl.deg2rad(75)]  # min, max beta_b [rad]

        options = self.calc_options()
        choice = self.chose_option(**options, **kwargs)
        shaft = self.calc_shaft(**choice)
        impeller = self.calc_impeller(**{**choice, **shaft})
        volute = self.calc_volute(**impeller)
        self.results = [options, choice, shaft, impeller, volute]

    def calc_options(self):
        """Calculate several design options of an impeller."""
        part = "---options---"
        np = []
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
        for i, p in enumerate(CN.NPOLES):
            n = op.rotational_speed(p, self.slip, self.hz)
            k = op.specific_speed(sh.angular_velocity(n), self.flow, self.head)
            # only specific speeds in the domain of centrifugal pumps
            if 0.2 <= k <= 1.2:
                np.append(p)
                rpm.append(n)
                cappa.append(k)
                phi.append(op.flow_number_poly(k))
                psi.append(op.head_number_poly(k))
                eta.append(op.efficency_poly(k))
                u_2.append(op.psi2u(psi[-1], self.head))
                d_2.append(im.diameter_omega(sh.angular_velocity(n), u_2[-1]))
                b_2.append(op.phi2b(d_2[-1], u_2[-1], phi[-1], self.flow))
                bd_2.append(op.width0diameter(b_2[-1], d_2[-1]))
                npsh_req.append(op.cappa2npsh(k, self.head))

        results = {}
        for i in ["part", "np", "rpm", "cappa", "phi", "psi", "eta", "u_2",
                  "d_2", "b_2", "bd_2", "npsh_req"]:
            results[i] = locals()[i]

        return results

    def chose_option(self, **kwargs):
        """Select an option of design according to a criteria."""
        cappa = kwargs["cappa"]
        np = kwargs["np"]

        part = "---chosen option---"
        if "fnp" in kwargs and kwargs["fnp"] in np:
            idx = np.index(kwargs["fnp"])
        else:
            # avoid cappa over .55 because it requires double curvature blades
            idx = cappa.index(max([val for val in cappa if val <= .55]))
        cappa = kwargs["cappa"][idx]
        np = kwargs["np"][idx]
        rpm = kwargs["rpm"][idx]
        phi = kwargs["phi"][idx]
        psi = kwargs["psi"][idx]
        eta = kwargs["eta"][idx]
        u_2 = kwargs["u_2"][idx]
        d_2 = kwargs["d_2"][idx]
        b_2 = kwargs["b_2"][idx]
        bd_2 = kwargs["bd_2"][idx]
        npsh_req = kwargs["npsh_req"][idx]
        eta_hyd = op.efficency_hyd_poly(cappa)
        eta_vol = op.efficency_vol_poly(cappa)

        results = {}
        for i in ["part", "np", "rpm", "cappa", "phi", "psi", "eta",
                  "u_2", "d_2", "b_2", "bd_2",
                  "npsh_req", "eta_hyd", "eta_vol"]:
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
        d_sh = sh.shaft_diameter(torque, self.tau_adm, coef=2)
        d_sh = round(d_sh, 3)
        d_hu = cl.bisect(sh.hub_diameter(d_sh), d_sh - 1, d_sh + 1, .001)
        d_hu = round(d_hu, 3)

        results = {}
        for i in ["part", "omega", "power", "torque", "d_sh", "d_hu"]:
            results[i] = locals()[i]

        return results

    def suction_eye(self, **kwargs):
        """Calculate suction eye of the impeller."""
        omega = kwargs["omega"]
        d_hu = kwargs["d_hu"]
        eta_vol = kwargs["eta_vol"]

        part_0 = "---impeller suction eye---"
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
            x_0.append(im.hub_blockage(d_0avg, d_hu))
            dif = abs(x_0[-1] - x_0[-2])
        d_0 = im.standard_diam(d_0avg)
        x_0 = im.hub_blockage(d_0, d_hu)

        results = {}
        for i in ["part_0", "d_0npsh", "d_0eff", "d_0flow", "d_0avg", "x_0",
                  "d_0"]:
            results[i] = locals()[i]

        return results

    def blade_trailing_edge(self, **kwargs):
        """Calculate blade trailing edge of the impeller."""
        omega = kwargs["omega"]
        eta_vol = kwargs["eta_vol"]
        d_hu = kwargs["d_hu"]
        d_0 = kwargs["d_0"]
        x_0 = kwargs["x_0"]
        d_2 = kwargs["d_2"]
        b_2 = kwargs["b_2"]
        x_2 = kwargs["x_2"]
        r_cvt = kwargs["r_cvt"]
        r_msl = kwargs["r_msl"]
        l_msl = kwargs["l_msl"]
        z = kwargs["z"]

        part_1 = "---impeller blade trailing edge---"
        beta_1b = None
        theta = []
        b = []
        n = 16
        for i in range(n):
            theta.append(im.angle_theta(n, i))
            d_imsl = im.streamline_diam(d_hu, d_0, theta[-1], r_msl)
            l_imsl = im.streamline_len(r_msl, theta=theta[-1])
            a_i = im.area(l_imsl, l_msl, d_0, x_0, d_2, b_2, x_2)
            b.append(im.width(d_imsl, a_i))
            if beta_1b is None:
                if "theta_1" in kwargs:
                    theta_1 = cl.deg2rad(kwargs["theta_1"])
                else:
                    theta_1 = theta[-1]
                gamma_1 = im.angle_gamma(r_cvt, r_msl, theta_1)
                if gamma_1 is None:
                    continue
                b_1 = b[-1]
                u_1 = im.blade_vel(omega, d_imsl)
                phi_1 = im.flow_number(d_imsl, b_1, u_1, self.flow)
                x_1 = [1]
                dif = 1
                err = .001
                while dif > err:
                    beta_1b = cl.bisect(im.angle_beta(u_1, phi_1, eta_vol,
                                                      x_1[-1], gamma_1),
                                        self.beta_b[0], self.beta_b[-1], .001)
                    x_1.append(im.blade_blockage(beta_1b, d_imsl, self.t, z))
                    dif = abs(x_1[-1] - x_1[-2])
                x_1 = x_1[-1]

        d_1msl = im.streamline_diam(d_hu, d_0, theta_1, r_msl)
        c_1m = im.meridional_abs_vel(u_1, phi_1)
        w_1 = im.relative_vel(c_1m, beta_1b)
        npsh_req = im.npsh_req(c_1m, w_1, self.lm, self.lw)

        results = {}
        for i in ["part_1", "theta_1", "b_1", "gamma_1", "beta_1b", "d_1msl",
                  "x_1", "u_1", "c_1m", "w_1", "npsh_req", "theta", "b"]:
            if i in ["theta", "theta_1", "gamma_1", "beta_1b"]:
                results[i] = cl.rad2deg(locals()[i])
            else:
                results[i] = locals()[i]

        return results

    def blade_leading_edge(self, **kwargs):
        """Calculate blade leading edge of the impeller."""
        omega = kwargs["omega"]
        eta_hyd = kwargs["eta_hyd"]
        eta_vol = kwargs["eta_vol"]
        phi = kwargs["phi"]
        u_2 = kwargs["u_2"]
        d_hu = kwargs["d_hu"]
        d_0 = kwargs["d_0"]
        z = kwargs["z"]

        part_2 = "---impeller blade leading edge---"
        u_2 = [u_2]
        dif = 1
        err = .001
        while dif > err:
            d_2 = im.diameter_omega(omega, u_2[-1])
            d_2 = round(d_2, 3)
            u_2.append(im.blade_vel(omega, d_2))
            dif = abs(u_2[-1] - u_2[-2])
        u_2 = u_2[-1]

        c_2m = im.meridional_abs_vel(u_2, phi)
        psi = im.head_number(u_2, self.head)
        psi_th = im.theoretic_head_number(psi, eta_hyd)

        d_msl = im.streamline_diam(d_hu, d_0)
        r_cvt = im.curvature_rad(d_2)
        r_msl = im.streamline_curv_rad(d_hu, d_0, r_cvt)
        l_msl = im.streamline_len(r_msl, d_2, d_msl)

        x_2 = [1]
        u_2sf = 0
        dif = 1
        err = .001
        while dif > err:
            beta_2b = cl.bisect(im.angle_beta(u_2, phi, eta_vol, x_2[-1], 0,
                                              psi_th, u_2sf),
                                self.beta_b[0], self.beta_b[-1], .001)
            b_2 = im.width(d_2, None, c_2m, self.flow, x_2[-1], eta_vol)
            b_2 = round(b_2, 3)
            phi = im.flow_number(d_2, b_2, u_2, self.flow)
            u_2sf = im.slip_factor(u_2, beta_2b, z)
            x_2.append(im.blade_blockage(beta_2b, d_2, self.t, z))
            dif = abs(x_2[-1] - x_2[-2])
        x_2 = x_2[-1]

        c_2m = im.meridional_abs_vel(u_2, phi)
        c_2u = im.circumferential_abs_vel(u_2, c_2m,  beta_2b)
        w_2 = im.relative_vel(c_2m, beta_2b)
        epsilon_ract = im.degree_reaction(phi, beta_2b, z)

        results = {}
        for i in ["part_2", "d_msl", "r_cvt", "r_msl", "l_msl",
                  "d_2", "u_2", "b_2",  "beta_2b",
                  "c_2m", "c_2u", "w_2", "u_2sf", "x_2",
                  "psi", "psi_th", "phi",  "epsilon_ract"]:
            if i == "beta_2b":
                results[i] = cl.rad2deg(locals()[i])
            else:
                results[i] = locals()[i]

        return results

    def calc_impeller(self, **kwargs):
        """Calculate the impeller."""
        suction = self.suction_eye(**kwargs)

        # attempt to reduce dif between beta_b angles changing z
        # or changing d_1 (through theta_1)
        for z in self.z:
            leading = self.blade_leading_edge(**{**kwargs, **{"z": z},
                                              **suction})
            trailing = self.blade_trailing_edge(**{**kwargs, **{"z": z},
                                                **suction, **leading})
            if abs(leading["beta_2b"] - trailing["beta_1b"]) <= 8:
                break
            else:
                if z == self.z[-1]:
                    skip = trailing["theta"].index(trailing["theta_1"]) + 1
                    for t in trailing["theta"][skip:-1]:
                        trailing = self.blade_trailing_edge(**{**kwargs,
                                                            **{"z": z},
                                                            **{"theta_1": t},
                                                            **suction,
                                                            **leading})
                        if abs(leading["beta_2b"] - trailing["beta_1b"]) <= 8:
                            break

        results = {**suction, **trailing, **leading, **{"z": z}}

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
                results[i] = cl.rad2deg(locals()[i])
            else:
                results[i] = locals()[i]

        return results


def main(**kwargs):
    """Print results."""
    prj = Project(**kwargs)

    for result in prj.results:
        for key, val in result.items():
            if isinstance(val, list):
                for k, v in enumerate(val):
                    if isinstance(v, tuple):
                        val[k] = tuple(round(t, 6) for t in v)
                    else:
                        val[k] = round(v, 6)
                print(key, " ", val)
            elif isinstance(val, (float, int)):
                print(key, " ", round(val, 6))
            else:
                print(val)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("flow", type=float, nargs="?", default=.011,
                        help="flow rate in [m^3/s]")
    parser.add_argument("head", type=float, nargs="?", default=25,
                        help="head in [m]")
    parser.add_argument("--hz", type=int, default=50,
                        help="frequency of alternating current [Hz]")
    parser.add_argument("--t", type=float, default=.003,
                        help="blade thickness [m]")
    parser.add_argument("--fnp", type=int, default=argparse.SUPPRESS,
                        help="force a different solution choosing a number of \
                        poles 'np' of the AC motor among the options")

    args = parser.parse_args()

    main(**vars(args))
