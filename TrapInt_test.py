import numpy as np
import math as ma
#import numexpr as ne
import scipy as sc

#from scipy.optimize import brentq #fmin

import scipy.integrate as intg
import scipy.signal as sgn

Nx = 100
Ny = 100


Args = [1.3,2.4]

# def f0v(x,y, args):
# 	A, B = args
# 	return A*x*x*y + B*y*y*x#*np.exp(-x*y)
# #f0v = np.vectorize(f0, excluded=[2])


# def fI1v(y, args):
# 	xs = np.linspace(0.0, 10., Nx)
# 	Trap_int = 2.*np.sum(f0v(xs,y,args)) - f0v(xs[0],y,args) - f0v(xs[Nx-1],y,args)
# 	return np.array(10./(2.*Nx) * (Trap_int))
# #fI1v = np.vectorize(fI1, excluded=[1])

# def fI2(args):
# 	ys = np.linspace(0.0, 10., Nx)
# 	Trap_int = 2.*np.sum(fI1v(ys,args)) - fI1v(ys[0],args) - fI1v(ys[Nx-1],args)
# 	return 10./(2.*Nx) * (Trap_int)


# print(fI2(Args))







# def g0(x,y, A,B):
# 	return A*x*x*y + B*y*y*x#*np.exp(-x*y)
# g0 = np.vectorize(g0)


# def gI1(y, A,B):
# 	xs = np.linspace(0.0, 10., Nx)
# 	Trap_int = 2.*np.sum(g0(xs,y,A,B)) - g0(xs[0],y,A,B) - g0(xs[Nx-1],y,A,B)
# 	return 10./(2.*Nx) * (Trap_int)
# gI1 = np.vectorize(gI1)

# def gI2(A,B):
# 	ys = np.linspace(0.0, 10., Nx)
# 	Trap_int = 2.*np.sum(gI1(ys,A,B)) - gI1(ys[0],A,B) - gI1(ys[Nx-1],A,B)
# 	return 10./(2.*Nx) * (Trap_int)


# print(gI2(Args[0], Args[1]))



def h0(x,y):
	return Args[0]*x*x*y + Args[1]*y*y*x#*np.exp(-x*y)
#h0 = np.vectorize(h0)


# def hI1(y):
# 	xs = np.linspace(0.0, 10., Nx)
# 	Trap_int = 2.*np.sum(h0(xs,y), axis=0) - h0(xs[0],y) - h0(xs[Nx-1],y)
# 	return np.array(10./(2.*Nx) * (Trap_int))
# #hI1 = np.vectorize(hI1)

# def hI2():
# 	ys = np.linspace(0.0, 10., Nx)
# 	Trap_int = 2.*np.sum(hI1(ys)) - hI1(ys[0]) - hI1(ys[Nx-1])
# 	return np.array(10./(2.*Nx) * (Trap_int))


# #print(hI2())

# def th1(y):
# 	xs = np.linspace(0.0, 10., Nx)
# 	return np.trapz(h0(xs,y))

xmax = 5.0
ymax = 10.0

def th2():
	xs = np.linspace(0.0, xmax, Nx)
	ys = np.linspace(0.0, ymax, Ny)
	xy = np.meshgrid(xs, ys)
	yx = np.meshgrid(ys, xs)
	xmax/(2.0*len(xs)) * (2.0 * np.sum(h0(xy[0],yx[0]), axis=0) - h0(xy[0][:][0],yx) - h0(xy[0][:][Nx-1],yx))

	#return np.trapz(np.trapz(h0(xy[0],yx[0])))

print(th2())


def h1(y):
	xs = np.linspace(0.0, xmax, Nx)
	return xmax/(2.0*len(xs)) * (2.0 * np.sum([h0(x,y) for x in xs]) - h0(xs[0],y) - h0(xs[Nx-1],y))

def h2():
	ys = np.linspace(0.0, ymax, Ny)
	return ymax/(2.0*len(ys)) * (2.0 * np.sum([h1(y) for y in ys]) - h1(ys[0]) - h1(ys[Ny-1]))

print(h2())

