#! /usr/bin/env python3

import numpy

tests = []
@tests.append
def f0(x):
	return x*x - 2, 2*x
@tests.append
def f1(x):
	return numpy.cos(x) - x, -numpy.sin(x) - 1


#----Diverging Test------
#abs(x)^(1/3) * sin(x)
@tests.append
def f3(x):
	return numpy.arctan(10*x) + x, (10.0/((100.0*x)**2 + 1)) + 1

#-----Newton Converges, Cubic Diverges-----
@tests.append
def f4(x):
	return numpy.sin(x+3) + 1.1*(x+3), numpy.cos(x+3) + 1.1



def newton(f, x):
	for i in range(100):
		fx, dfx = f(x)
		if numpy.abs(fx) < 1e-12:
			return i, x, fx, fx/dfx
		print (f.__name__, i, x, fx, dfx)
		x -= fx / dfx

'''
def newton(f, x):
	fx = f(x)[0]
	i = 0
	while (numpy.abs(fx) > 1e-12):
		fx, dfx = f(x)
		print (f.__name__, i, x, fx, dfx)
		i += 1
		x -= fx / dfx
'''


if __name__ == "__main__":
	for f in tests:
		print(f.__name__, newton(f, 1))
