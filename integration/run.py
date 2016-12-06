#! /usr/bin/env python3
import numpy

tests = []
@tests.append
class f_expx:
    def F(x):
        return numpy.exp(2 * x) / (1 + x**2)
    def f(x):
        return 2 * numpy.exp(2 * x) / (1 + x**2) - numpy.exp(2 * x) / (1 + x**2)**2 * 2 * x
@tests.append
class f_quadratic:
    def F(x):
        return x**3 / 3 - x**2 + 3 * x - 1
    def f(x):
        return x**2 - 2 * x + 3
@tests.append
class f_tanh:
    def F(x):
        return numpy.tanh(x)
    def f(x):
        return numpy.cosh(x)**(-2)
@tests.append
class f_cos:
    def F(x):
        return -20*(x-numpy.cos(x))
    def f(x):
        return -20*(numpy.sin(x) + 1)
"""
This function 'technicaly' breaks the
recursive extrapolation method.
It causes the maximum recursion depth to be exceeded.
This is due to the function breaking newton's method by diverging.
Taken from my Hw1 solution

@tests.append
class f_arctan:
    def F(x):
        return numpy.arctan(10*x) + x
    def f(x):
        return (10.0/((100.0*x)**2 + 1)) + 1
"""

# From the notebook
def fint_trapezoid(f, a, b, n=20):
    dx = (b - a) / (n - 1)     # We evaluate endpoints of n-1 intervals
    x = numpy.linspace(a, b, n)
    fx = f(x)
    fx[[0, -1]] *= .5
    # fx[0] *= .5; fx[-1] *= .5
    return numpy.sum(fx) * dx

def fint_midpoint(f, a, b, n=20):
    dx = (b - a)/n     # Width of each interval
    x = numpy.linspace(a+dx/2, b-dx/2, n)
    return numpy.sum(f(x))*dx

#From the notebook
def rec_extrap(t, a, b, n):
    I_h = fint_trapezoid(t.f, a, b, n)
    I_2h = fint_trapezoid(t.f, a, b, n // 2)
    I_extrap = I_h + (I_h - I_2h) / 3
    I_exact = t.F(b) - t.F(a)
    error = I_extrap - I_exact
    if abs(error) < 1e-4:  # Recurse until the errror is less than 1e-4
        print ('{:12s} {:12s} {: 4d}, {: 10e}'.format(
            t.__name__, "Recursive Extrapolation", n, error))
        return n
    else:
        return rec_extrap(t, a, b, (n + 10))

#From the notebook
def legendre_eval(x, n):
    P = numpy.ones((len(x), n))
    dP = numpy.ones((len(x), n))
    if n > 1:
        P[:, 1] = x
    for k in range(1, n - 1):
        #Calculate the legrendre polynomials and derivatives
        P[:, k + 1] = ((2 * k + 1) * x * P[:, k] - k * P[:, k - 1]) / (k + 1)
        dP[:, k + 1] = (2 * k + 1) * P[:, k] + dP[:, k - 1]
    return P, dP

#using newtons method to find the roots of the polynomials
def legendre_roots(n):
    x = numpy.cos(numpy.linspace(.5 / (n - 1), 1 - .5 / (n - 1), n) * numpy.pi)
    f, df = legendre_eval(x, n)
    rootsList = numpy.zeros(n)
    for i in range(0, n - 1):
        for j in range(100):
            fx = f[int(x[i])][n - 1]
            dfx = df[int(x[i])][n - 1]
            if numpy.abs(f[int(x[i])][n - 1]) < 1e-8:
                rootsList[i] = x[i]
                break
            else:
                x[i] -= fx/dfx
                f, df = legendre_eval(x, n)
    return rootsList


def fint_legendre(t, a, b, n):
    beta = legendre_roots(n)
    T = numpy.diag(beta, -1) + numpy.diag(beta, 1)
    D, V = numpy.linalg.eig(T)
    w = V[0, :]**2 * (b - a)
    x = (a + b) / 2 + (b - a) / 2 * D
    I_calc = w.dot(t.f(x))
    I_exact = t.F(b) - t.F(a)
    error = I_calc - I_exact
    print ('{:12s} {:12s} {: 4d}, {: 10e}'.format(
        t.__name__, "Gauss-Legendre Method  ", n, error))
    return
    
for t in tests:
    a, b = -2, 2
    n = 10
    fint_legendre(t, a, b, n)
for t in tests:
    a, b = -2,2
    n = 10
    rec_extrap(t, a, b, n)
