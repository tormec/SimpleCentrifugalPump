import math
import constants as CNST


class Options(object):
    """Methods to calculate several design options of an impeller."""

    def rotational_speed(self, pp, slip, hz):
        """Calculate rotational speed at different pole pairs of an AC motor.

        :param pp (int): pole pairs of an AC motor
        :param slip (int): slip factor [%]
        :param hz (int): utility frequency [Hz]
        :return rpm (float): rotational speed [rpm]
        """
        rpm = 120 * hz / pp * (1 - slip / 100)

        return rpm

    def angular_velocity(self, rpm):
        """Calculate angular velocity of a shaft.

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
        :return cappa (float): typical number
        """
        omega = self.angular_velocity(rpm)
        cappa = omega * flow**0.5 / (CNST.G * head)**0.75

        return cappa

    def cappa2rpm(self, cappa, flow, head):
        """Calculate rotational speed for a given typical number.

        :param cappa (float): typical number
        :param flow (float): flow rate [m^3/s]
        :param head (float): head [m]
        :return rpm (float): rotational speed [rpm]
        """
        omega = cappa * (CNST.G * head)**0.75 / flow**0.5
        rpm = omega * 60 / (2 * math.pi)

        return rpm

    def rpm2pp(self, rpm, slip, hz):
        """Calculate polar pairs of an AC motor for a given rotational speed.

        :param rpm (float): rotational speed [rpm]
        :param slip (int): slip factor [%]
        :param hz (int): utility frequency [Hz]
        :return pp (int): polar pairs
        """
        pp = 120 * hz / rpm * (1 - slip / 100)

        return pp

    def flow_number(self, cappa):
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

    def head_number(self, cappa):
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

    def efficency(self, cappa):
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

    def peripheral_velocity(self, psi, head):
        """Calculate peripheral velocity function of the head number.

        :param psi (float): head number
        :param head (float): head [m]
        :return u (float): peripheral velocity [m/s]
        """
        u = (CNST.G * head / psi)**0.5

        return u

    def diameter(self, u, rpm):
        """Calculate diameter function of the rotationl speed.

        :param u (float): peripheral velocity [m/s]
        :param rpm (float): rotational speed [rpm]
        :return d (float): diameter [m]
        """
        d = (60 * u) / (math.pi * rpm)

        return d

    def width(self, u, d, phi, flow):
        """Calculate impeller width function of the flow number.

        :param u (float): peripheral velocity [m/s]
        :param d (float): diameter [m]
        :param phi (float): flow number
        :param flow (float): flow rate [m^3/s]
        :return b (float): impeller width [m]
        """
        b = flow / (math.pi * d * u * phi)

        return b

    def width0diameter(self, b, d):
        """Calculate rate width over diameter.

        :param b (float): impeller width [m]
        :param d (float): diameter [m]
        :return bd (float): width over diameter
        """
        bd = b / d

        return bd

    def npsh_r(self, cappa, head):
        """Calculate the neat positive suction head required.

        :param cappa (float): typical number
        :param head (float): head [m]
        :return npsh_r (float): neat positive suction head required [m]
        """
        npsh_r = .25 * cappa**(4/3) * head

        return npsh_r
