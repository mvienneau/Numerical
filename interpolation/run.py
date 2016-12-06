#!/usr/bin/env python3
import numpy

def piecewise_linear_interp_and_eval(x, xx):
    """
    Essentially a modified version of the cubic spline interp from the notes, but linear
    """
    def splineInterp(y):
        sortedIndex = numpy.argsort(x)
        sort_x, sort_y = x[sortedIndex], y[sortedIndex]

        #Get the change in X and change in Y
        changeX = sort_x[1:] - sort_x[:-1]
        changeY = sort_y[1:] - sort_y[:-1]

        #Create the spline matrix
        splineMatrix = numpy.zeros((len(x) - 1, 3))
        splineMatrix[:, 0] = sort_x[:-1]
        splineMatrix[:, 1] = changeY/changeX
        splineMatrix[:, 2] = sort_y[:-1]

        return splineMatrix

    def splineEval(splineMat, xx):
        leftVal = splineMat[:, 0].searchsorted(xx) - 1
        leftVal[leftVal < 0] = 0  # Use the leftmost interval even if xx<=x
        f = numpy.zeros_like(xx)
        for i, t in enumerate(xx):
            cur_Spline = splineMat[leftVal[i]]
            X = xx[i] - cur_Spline[0]
            f[i] = numpy.polyval(cur_Spline[1:], X)
        return f
    # Build an explicit matrix for the spline fit evaluated at xx
    A = numpy.zeros((len(xx), len(x)))
    for i, e in enumerate(numpy.eye(len(x), len(x))):
        S = splineInterp(e)
        A[:, i] = splineEval(S, xx)
    return A


def cosspace(a, b, n=50):
    return (a + b)/2 + (b - a)/2 * (numpy.cos(numpy.linspace(0, numpy.pi, n)))


def vander_chebyshev(x, n=None):
    if n is None:
        n = len(x)
    T = numpy.ones((len(x), n))
    if n > 1:
        T[:, 1] = x
    for k in range(1, n - 1):
        T[:, k + 1] = 2 * x * T[:, k] - T[:, k - 1]
    return T


def chebyshev_interp_and_eval(x, xx):
    """Matrix mapping from values at points x to values
    of Chebyshev interpolating polynomial at points xx"""
    A = vander_chebyshev(x)
    B = vander_chebyshev(xx, len(x))
    return B.dot(numpy.linalg.inv(A))

def f(x):
    return numpy.exp(-16*x**2)

def g(x):
    return abs(x)

def maxerror(interp_and_eval, f, xspace, interval, npoints, eps):
    for n in npoints:
        x = xspace(*interval, n)
        xx = numpy.linspace(*interval, 300)
        A = interp_and_eval(x, xx)
        max_error = numpy.linalg.norm(A.dot(f(x)) - f(xx), numpy.inf)
        if max_error < eps:
            return n, max_error
    return n, max_error


if __name__ == '__main__':
    """NORM TESTING"""
    x = numpy.linspace(-1, 1, 20)
    xx = numpy.linspace(-1, 1, 100)
    print("piecewise_linear linspace -> linspace: norm=%e" % numpy.linalg.norm(piecewise_linear_interp_and_eval(x, xx)))

    x = cosspace(-1, 1, 20)
    print("piecewise_linear cosspace -> linspace: norm=%e" % numpy.linalg.norm(piecewise_linear_interp_and_eval(x, xx)))

    x = numpy.linspace(-1, 1, 20)
    xx = numpy.linspace(-1, 1, 100)
    print("chebyshev linspace -> linspace: norm=%e" % numpy.linalg.norm(chebyshev_interp_and_eval(x, xx)))

    x = cosspace(-1, 1, 20)
    print("chebyshev cosspace -> linspace: norm=%e" % numpy.linalg.norm(chebyshev_interp_and_eval(x, xx)))

    """MAX ERROR"""
    epsilon = 1e-4
    num_points = numpy.arange(2, 101)
    n, max_error = maxerror(piecewise_linear_interp_and_eval, f, cosspace, (-1, 1), num_points, epsilon)
    print("piecewise_linear f() cosspace -> linspace: n=%d max_error=%e" % (n, max_error))

    n, max_error = maxerror(piecewise_linear_interp_and_eval, f, numpy.linspace, (-1, 1), num_points, epsilon)
    print("piecewise_linear f() linspace -> linspace: n=%d max_error=%e" % (n, max_error))

    epsilon = 1e-2
    n, max_error = maxerror(piecewise_linear_interp_and_eval, g, cosspace, (-1, 1), num_points, epsilon)
    print("piecewise_linear g() cosspace -> linspace: n=%d max_error=%e" % (n, max_error))

    n, max_error = maxerror(piecewise_linear_interp_and_eval, g, numpy.linspace, (-1, 1), num_points, epsilon)
    print("piecewise_linear g() linspace -> linspace: n=%d max_error=%e" % (n, max_error))
