#!/usr/bin/python
from math import *
import sys
import os

Nx=20
Ny=20
Nz=20
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
G.PPhost.EquilibriumOrder=2;
G.PPhost.IsothermalRelaxation=1;

G.PPhost.dr=1.0;
G.PPhost.dt=1.0;

G.PPhost.visc_atT=0.1;

G.PPhost.initial.planeTVG=0;
G.PPhost.initial.u0=0.2;
G.PPhost.initial.rho0=1;

G.PPhost.initial.uDragX=0.0;
G.PPhost.initial.uDragY=0;
G.PPhost.initial.beta0=10.0#1.0;
G.PPhost.initial.r0=10.0;
G.PPhost.initial.T0=0.1#1/3.#0.8#1./3;

G.PPhost.StepIterPeriod = 1
G.PPhost.set_drop_dir(drop_into);

lbm.run(sys.argv)

############################################################


