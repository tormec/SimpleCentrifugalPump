"""Methods to calculate the pump shaft."""

import math
import lib.constants as CN


def rotational_speed(np, slip, hz):
    """Calculate rotational speed at different number of poles of an AC motor.

    :param np (int): number of poles of an AC motor
    :param slip (int): slip factor [%]
    :param hz (int): utility frequency [Hz]
    :return rpm (float): rotational speed [rpm]
    """
    rpm = 120 * hz / np * (1 - slip / 100)

    return rpm


def angular_velocity(rpm):
    """Calculate angular velocity.

    :param rpm (float): rotational speed [rpm]
    :return omega (float): angular velocity [rad/s]
    """
    omega = 2 * math.pi * rpm / 60

    return omega


def power(eta_tot, flow, head):
    """Calculate power the el. motor should provide at pump shaft.

    :param eta_t (float): total efficency
    :param flow (float): flow rate [m^3/s]
    :param head (float): head [m]
    :return power (float): power [W]
    """
    power = (CN.RHO * flow * CN.G * head) / eta_tot

    return power


def torque(power, omega):
    """Calculate torque at pump shaft.

    :param power (float): power [W]
    :param omega (float): angular velocity [rad/s]
    :return torque (float): torque [Nm]
    """
    torque = power / omega

    return torque


def shaft_diameter(torque, tau_adm, coef=1):
    """Calculate shaft diameter.

    :param torque (float): torque [Nm]
    :param tau_adm (float): tau admissible [MPa]
    :param coef (int): security coeffient
    :return d_sh (float): shaft diameter [m]
    """
    d_sh = ((16 * coef * torque) / (math.pi * (tau_adm * 10**6)))**(1/3)

    return d_sh
