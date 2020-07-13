from math import *
L = 100.
nu = .2
kx = 2.*pi/L
T = 1./3.
uD = 0.1*sqrt(T) #if I know it
dtscale = 1.#/10. #if I know it

vconv = dtscale

ky = kx
kz = kx
wl = 2*pi/sqrt(kx**2+ky**2+kz**2)
decr = nu*(kx**2+ky**2+kz**2)*T

halftime = - log(0.5)/decr

import numpy as np

#the decay rate of the transverse wave is k**2*T*(tau-0.5)

import sys
#y = np.loadtxt(sys.argv[1],usecols=0)
y = np.loadtxt( sys.stdin,usecols=9)
uD = sqrt((y[0]/vconv)**2-0.05**2) / sqrt(T)
#print(f"#u_d may be {uD}")
x = np.arange(len(y))
A = np.sqrt((y/vconv)**2-uD**2*T)
lkoff = np.polyfit(x,np.log( A ),1)

nu_calc = -lkoff[0]/(kx**2+ky**2+kz**2)/T/vconv

print(f"# sure by {len(y)}/{halftime} = {len(y)/halftime*100} %")
print(f"{uD} {nu:16.9} {nu_calc:16.9}")
