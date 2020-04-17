import math


def bisect(f, a, b, err):
    """Implement bisection method to find the root of f(x) in [a, b].

    f(x) must be a continuos function in [a, b] and f(a) * f(b) < 0

    :param f (function): function to evaluate as def f(x): return 2 * x + 1
    :param a (float): left end
    :param b (float): right end
    :param err (float): tollerance
    """
    n_max = 1e2
    n = 0
    amp = err + 1
    fa = f(a)
    while amp >= err and n < n_max:
        n = n + 1
        amp = abs(b - a)
        x = a + amp * .5  # middle point between a and b
        fx = f(x)
        if fa * fx < 0:
            b = x
        elif fa * fx > 0:
            a = x
            fa = fx
        else:
            amp = 0

    return x


def rad2deg(rad):
    """Convert radians to degrees.

    :param rad (float/list): angle(s) in radians [rad]
    :return deg (float/list): angle(s) in degrees [deg]
    """
    if type(rad) == list:
        deg = [round(math.degrees(i), 2) for i in rad]
    else:
        deg = round(math.degrees(rad), 2)

    return deg


class Test(object):
    """Test methods"""
    def f(self, x):
        return 2 * x + 1

    def test_bisect(self, f, a, b, err):
        assert bisect(f, a, b, err) == -.5

    def test_rad2deg(self, rad):
        if type(rad) == list:
            assert rad2deg(rad) == [0, 30, 45, 60, 90]
        else:
            assert rad2deg(rad) == 45


if __name__ == "__main__":
    test = Test()
    test.test_bisect(test.f, -2, 2, .001)
    test.test_rad2deg([0, math.pi / 6, math.pi / 4, math.pi / 3, math.pi / 2])
    test.test_rad2deg(math.pi / 4)
