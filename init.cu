#include "structs.cuh"
#include "init.h"
#include "LBMconsts.cuh"
#include "phys.h"
#include <nvfunctional>

#include "materials.cuh"

template<class F> __global__ void fill(F);

void init(){
  parsHost.iStep=0;
  parsHost.StepsMax=PPhost.MaxSteps;
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );

  printf("Malloc data\n");
  parsHost.data.malloc_data(Nx,Ny,Nz);
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );

  cuTimer init_timer;
  //fill<<<dim3(Nx,Ny),Nz>>>( [] __device__(int ix, int iy,int iz) {return blank_mat(ix,iy,iz);} );
  fill<<<dim3(Nx,Ny),Nz>>>( [] __device__(int ix, int iy,int iz) {
     //return sinTemperature(ix,iy,iz);
     return shear_wave(ix,iy,iz);
     //return vortex_mat(ix,iy,iz);
     //return TGV_mat(ix,iy,iz);
     //return blank_mat(ix,iy,iz);
  } );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  printf("\n");
  printf("TLat=%16.9g\n\n",LBMconsts::TLat);

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
  #ifdef SEMI_LAGRANGIAN
  c.calcEq(feq, rho_uT.first, make_ftype3(0,0,0), LBMconsts::cs2 );
  #else
  c.calcEq(feq, rho_uT.first, vel, rho_uT.second.w );
  #endif
  for(int iq=0; iq<LBMconsts::Qn; iq++) c.f[iq]=feq[iq];
  c.rho = rho;
  c.vel = vel;
  c.T = T;
  for(int iq=0; iq<LBMconsts::Qn; iq++) {
    if(c.f[iq]<=0) printf("Warning CELL (%4d %4d %4d): f[%d] is negative\n", ix,iy,iz, iq);
  }

  pars.data.set_cell(c, 0, ix,iy,iz);
  pars.data.set_cell(c, 1, ix,iy,iz);

  //using namespace LBMconsts;
  //if(ix==0 && iy==0 && iz==0) printf("TLAT = %g W0123 = %g %g %g %g\n", TLat, W0get(TLat), W1get(TLat), W2get(TLat), W3get(TLat));

}

void PhysPars::setupUnits(){
  RhoUnitConv = 1;
  VelUnitConv = dt/dr;
  ViscUnitConv = dt/(dr*dr);
  const ftype ViscAtTUnitConv = 1./dt;
  tau = 0.5+visc_atT*ViscAtTUnitConv;
  dtau = 1/tau;
  TempUnitConv = dt*dt/(dr*dr);
  printf("Density Units Converstion Coeff = %g\n", RhoUnitConv);
  printf("Velocity Units Converstion Coeff = %g\n", VelUnitConv);
  printf("Temperature Units Converstion Coeff = %g\n", TempUnitConv);
  printf("Viscosity Units Converstion Coeff = %g\n", ViscUnitConv);
  printf("Tau relaxation = %g\n", tau);
}
