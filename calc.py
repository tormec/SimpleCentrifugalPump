import numpy as np


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


def polyfit(x, fx, w, n):
    """Calculate the coeffients of the polynomial that best approximate (x, y).

    :param x (list): indipendent variable
    :param fx (list): dipendent variable
    :param w (list): weight of each node
    :param n (int): degree of the polynomial
    :return a (list): coeffientes of the polynomial
    """
    A = []
    for i in range(len(x)):  # rows
        r = []
        for j in range(n, 0, -1):  # columns
            r.append(w[i]**.5 * x[i]**(n-j))
        A.append(r)

    b = []
    for i in range(len(x)):  # columns
        b.append([w[i]**.5 * fx[i]])

    A = np.array(A)
    b = np.array(b)

    a = np.linalg.lstsq(np.dot(A.T, A), np.dot(A.T, b))  # (A' * A) \ (A' * b)
    a = a[0].T.tolist()[0]

    return a


class Test(object):
    """Test methods"""
    def f(self, x):
        return 2 * x + 1

    def test_bisect(self, *args):
        assert bisect(*args) == -.5

    def test_polyfit(self, *args):
        a = [round(v, 2) for i, v in enumerate(polyfit(*args))]
        assert a == [-0.65, 3.13]


if __name__ == "__main__":
    test = Test()
    test.test_bisect(test.f, -2, 2, .001)
    test.test_polyfit([-1, 2, 5, 6], [-3, 5, 12, 21], [1, 1, 1, 1], 2)
