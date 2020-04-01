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


class Test(object):
    """Test methods"""
    def f(self, x):
        return 2 * x + 1

    def test_bisect(self, *args):
        assert bisect(*args) == -.5


if __name__ == "__main__":
    test = Test()
    test.test_bisect(test.f, -2, 2, .001)
