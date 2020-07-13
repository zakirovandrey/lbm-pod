#!/usr/bin/python
from math import *
import sys
import os

Nx=100
Ny=100
Nz=100
FloatPrecision=2

MPIon=0
execfile('compileSource.py')

import lbm
from lbm import G

Nx,Ny,Nz = G.Nx,G.Ny,G.Nz
print "Grid sizes: %d x %d x %d"%(Nx,Ny,Nz)
G.PPhost.setDefault()
G.PPhost.stencilInterpWidth = 2;
G.PPhost.stencilFixed = 0;
G.PPhost.RegOrder = 2
G.PPhost.EquilibriumOrder = 4;
G.PPhost.IsothermalRelaxation = 1;

G.PPhost.dr=1.0;
G.PPhost.dt=1.0;

G.PPhost.visc_atT=.2;

Tinit=1./3;

Ma_a = 10.

#G.PPhost.initial.shearWaveDir = 3; ### 2 - XY, 3 - XYZ
#G.PPhost.initial.uDragX = (Ma_a/sqrt(3))*sqrt(Tinit);
#G.PPhost.initial.uDragY = (Ma_a/sqrt(3))*sqrt(Tinit);
#G.PPhost.initial.uDragZ = (Ma_a/sqrt(3))*sqrt(Tinit);
waveD = 3
G.PPhost.initial.shearWaveDir = waveD; ### 2 - XY, 3 - XYZ
#G.PPhost.initial.uDragX = (Ma_a/sqrt(2))*sqrt(Tinit) if G.PPhost.initial.shearWaveDir==2 else Ma_a*sqrt(Tinit);
#G.PPhost.initial.uDragY = (Ma_a/sqrt(2))*sqrt(Tinit) if G.PPhost.initial.shearWaveDir==2 else 0;
#G.PPhost.initial.uDragZ = 0;
G.PPhost.initial.uDragX = (Ma_a/sqrt(waveD))*sqrt(Tinit)
G.PPhost.initial.uDragY = (Ma_a/sqrt(waveD))*sqrt(Tinit) if waveD>1 else 0
G.PPhost.initial.uDragZ = (Ma_a/sqrt(waveD))*sqrt(Tinit) if waveD>2 else 0
G.PPhost.initial.u0 = 0.05;
G.PPhost.initial.rho0 = 1;
G.PPhost.initial.A0 =0.001;

G.PPhost.initial.T0=Tinit
G.PPhost.fixedTemperature=1

G.PPhost.StepIterPeriod = 1
G.PPhost.MaxSteps = 600
G.PPhost.set_drop_dir(drop_into);

lbm.run(sys.argv)

############################################################


