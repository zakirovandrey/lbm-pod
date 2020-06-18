#include "structs.cuh"
#include "init.h"
#include "LBMconsts.cuh"
#include "phys.h"
#include <nvfunctional>

#include "materials.cuh"

template<class F> __global__ void fill(F);

void init(){
  parsHost.iStep=0;
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );

  printf("Malloc data\n");
  parsHost.data.malloc_data(Nx,Ny,Nz);
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );

  cuTimer init_timer;
  fill<<<dim3(Nx,Ny),Nz>>>( [] __device__(int ix, int iy,int iz) {return blank_mat(ix,iy,iz);} );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  printf("\n");

  printf("Initialization time: %.2f ms\n", init_timer.gettime());
  
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );

}

template<class F> __global__ void fill(F func){
  /*int ix = threadIdx.x+blockIdx.x*blockDim.x;
  int iy = threadIdx.y+blockIdx.y*blockDim.y;
  int iz = threadIdx.z+blockIdx.z*blockDim.z;*/
  int ix = blockIdx.x;
  int iy = blockIdx.y;
  int iz = threadIdx.x;
  Cell c;
  std::pair<ftype, ftype4> rho_uT = func(ix,iy,iz);

  const ftype rho = rho_uT.first;
  const ftype3 vel = make_ftype3(rho_uT.second.x, rho_uT.second.y, rho_uT.second.z) ;
  const ftype T = rho_uT.second.w;

  assert(rho_uT.second.w>LBMconsts::Tmin);
  assert(rho_uT.second.w<LBMconsts::Tmax);

  ftype feq[LBMconsts::Qn];
  c.calcEq(feq, rho_uT.first, make_ftype3(0,0,0), rho_uT.second.w );
  for(int iq=0; iq<LBMconsts::Qn; iq++) c.f[iq]=feq[iq];
  c.rho = rho;
  c.vel = vel;
  c.T = T;

  pars.data.set_cell(c, 0, ix,iy,iz);
  pars.data.set_cell(c, 1, ix,iy,iz);
}

