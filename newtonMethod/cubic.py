#! /usr/bin/env python3

import numpy

tests = []
@tests.append
def f0(x):
	return x*x - 2, 2*x, 2
@tests.append
def f1(x):
	return numpy.cos(x) - x, -numpy.sin(x) - 1, -numpy.cos(x)

@tests.append
def f3(x):
	return numpy.arctan(10*x) + x, (10.0/((100.0*x)**2 + 1)) + 1, -(2000*x)/((100*(x**2)+1)**2)

@tests.append
def f4(x):
	return numpy.sin(x+3) + 1.1*(x+3), numpy.cos(x+3) + 1.1, -numpy.sin(x+3)



def cubic(f, x):
	for i in range(100):
		fx, dfx, ddfx = f(x)
		sqrt = ((dfx**2)-(2*(fx*ddfx)))**(1/2)
		if numpy.abs(fx) < 1e-12:
			return i, x, fx
		print (f.__name__, i, x, fx, dfx)
		posRoot = x + ((1 / ddfx * (-dfx+sqrt)))
		negRoot = x + ((1 / ddfx * (-dfx-sqrt)))
		if (abs(f(posRoot)[0]) < abs(f(negRoot)[0])):
			x = posRoot
		else:
			x = negRoot




if __name__ == "__main__":
	for f in tests:
		print(f.__name__, cubic(f, 1))
