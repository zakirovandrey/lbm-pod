#ifndef _STRUCTS_H
#define _STRUCTS_H
#include <cuda.h>
#include <curand_kernel.h>
#include "params.h"
#include "data.cuh"

//---------------------------Main Calculation parameters------------------------------------//
struct AllParams{
  unsigned int nFunc;
  int checkFlag;
  int Nt;
  int3 Nroot;
  int iStep;
  double curtime;
  Data_t data;
};

struct AllParamsHost: public AllParams {
  Arr3D_pars arr4im, arr4surf;
  unsigned int MaxFunc;
  void set();
  void reset_im() {
    int ndevs=0; int curdev;
    CHECK_ERROR( cudaGetDevice(&curdev) );
    CHECK_ERROR( cudaGetDeviceCount(&ndevs) );
    #ifdef DRAW_ONLY_LEVEL
    arr4im.reset(Nx, Ny, Nz);
    arr4im.BufSize = sizeof(float)*Nx*Ny*Nz;
    #else
    arr4im.reset(Nx, Ny, Nz);
    arr4im.BufSize = sizeof(float)*Nx*Ny*Nz;
    #endif
    CHECK_ERROR( cudaMallocHost((void**) (&arr4im.Arr3Dbuf), arr4im.BufSize) );
    CHECK_ERROR( cudaMemset(arr4im.Arr3Dbuf, 0, arr4im.BufSize) );
//     CHECK_ERROR( cudaMallocHost((void**) (&arr4surf.Arr3Dbuf), arr4surf.BufSize) );
//     CHECK_ERROR( cudaMemset(arr4surf.Arr3Dbuf, 0, arr4surf.BufSize) );
    arr4im.inGPUmem = true;
    arr4surf.inGPUmem = true;
    CHECK_ERROR(cudaSetDevice(curdev));
  }
  void clear() {
    cudaFree (arr4im.Arr3Dbuf);
    arr4im.clear();
  }
  void checkSizes() { }
};
extern AllParamsHost parsHost;
extern __constant__ AllParams pars;
//----------------------------------------------------------------------------------//

int _main(int argc, char** argv);

#endif //_STRUCTS_H
