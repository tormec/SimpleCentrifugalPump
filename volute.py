import math


class Volute(object):
    """Methods to calculate the volute."""

    def absolute_velocity_throat(self, cu1):
        """Calculate absolute velocity at throat section.

        :param cu1 (float): peripheral velocity component [m/s]
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
