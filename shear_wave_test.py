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
G.PPhost.stencilFixed=1;
G.PPhost.RegOrder=-1

G.PPhost.dr=1.0;
G.PPhost.dt=1.0;

G.PPhost.visc_atT=0.1;

Tinit=0.1;

G.PPhost.initial.u0=0.1*sqrt(Tinit);
G.PPhost.initial.rho0=1;

G.PPhost.initial.T0=Tinit

G.PPhost.StepIterPeriod = 10
G.PPhost.MaxSteps = 100
G.PPhost.set_drop_dir(drop_into);

lbm.run(sys.argv)

############################################################


