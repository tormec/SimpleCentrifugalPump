import math
import constants as CNST


class Shaft(object):
    """Methods to calculate the pump shaft."""

    def power(self, eta_tot, flow, head):
        """Calculate power the el. motor should provide at pump shaft.

        :param eta_t (float): total efficency
        :param flow (float): flow rate [m^3/s]
        :param head (float): head [m]
        :return power (float): power [W]
        """
        power = (CNST.RHO * flow * CNST.G * head) / eta_tot

        return power

    def torque(self, power, omega):
        """Calculate torque at pump shaft.

        :param power (float): power [W]
        :param omega (float): angular velocity [rad/s]
        :return torque (float): torque [Nm]
        """
        torque = power / omega

        return torque

    def shaft_diameter(self, torque, tau_adm):
        """Calculate shaft diameter.

        :param torque (float): torque [Nm]
        :param tau_adm (int): tau admissible [MPa]
        :return d_sh (float): shaft diameter [m]
        """
        d_sh = ((32 * torque) / (math.pi * (tau_adm * 10**6)))**(1/3)

        return d_sh
