#! /usr/bin/env python3

import numpy

tests = []
@tests.append
def f0(x):
	return x*x - 2, 2*x
@tests.append
def f1(x):
	return numpy.cos(x) - x, -numpy.sin(x) - 1

@tests.append
def f3(x):
	return numpy.arctan(10*x) + x, (10.0/((100.0*x)**2 + 1)) + 1


def newton(f, x):
	for i in range(100):
		s = 2
		fx, dfx = f(x)
		if numpy.abs(fx) < 1e-12:
			return i, x, fx
		print (f.__name__, i, x, fx, dfx)
		newX = x - (fx / dfx)
		while (numpy.isnan(newX) or abs(f(newX)[0]) > abs(f(x)[0])):
			newX = x - (f(x)[0]/s*f(x)[1])
			s *= 2
		x = newX









if __name__ == "__main__":
	for f in tests:
		print(f.__name__, newton(f, 1))
