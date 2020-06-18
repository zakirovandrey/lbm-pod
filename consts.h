#ifndef _CONSTS_H
#define _CONSTS_H

#ifndef SET_FROM_PYTHON
#define USE_FLOAT

const int Nx=320;  //Nx%16=0
const int Ny=128;  //Ny%8=0
const int Nz=128;  //Nz%4==0
#else // defined SET_FROM_PYTHON
const int Nx=NX;  //Nx%16=0
const int Ny=NY;  //Ny%8=0
const int Nz=NZ;  //Nz%4==0
#endif // SET_FROM_PYTHON

#ifdef USE_FLOAT
typedef float ftype;
#elif defined USE_DOUBLE
typedef double ftype;
#endif


static int CudaDevs;

#endif//_CONSTS_H
