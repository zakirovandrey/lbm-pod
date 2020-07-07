#!/usr/bin/python
from math import *
import sys
import os

Nx=200
Ny=200
Nz=1
FloatPrecision=2

MPIon=0
execfile('compileSource.py')

import lbm
from lbm import G

Nx,Ny,Nz = G.Nx,G.Ny,G.Nz
print "Grid sizes: %d x %d x %d"%(Nx,Ny,Nz)
G.PPhost.setDefault()
G.PPhost.stencilInterpWidth=2;
G.PPhost.stencilFixed=0;
G.PPhost.RegOrder = 4

G.PPhost.dr=1.0;
G.PPhost.dt=1.0;

G.PPhost.visc_atT=.2;

Tinit=1./3.;

Ma_a = 0.

G.PPhost.initial.uDragX = (Ma_a/sqrt(2))*sqrt(Tinit);
G.PPhost.initial.uDragY = (Ma_a/sqrt(2))*sqrt(Tinit);
G.PPhost.initial.u0 = 0.05;
G.PPhost.initial.rho0 = 1;

G.PPhost.initial.T0=Tinit
G.PPhost.fixedTemperature=1

G.PPhost.StepIterPeriod = 1
G.PPhost.MaxSteps = 3000
G.PPhost.set_drop_dir(drop_into);

lbm.run(sys.argv)

############################################################


