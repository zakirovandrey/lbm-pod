#include "structs.cuh"
#include "semaphore.h"
#include "LBMconsts.cuh"
#include "phys.h"

#include "collision.cu"
#include "streaming.cu"

void calcLBM(int it, std::vector<double>& timings);
void simple_drop();
void debug_print();

struct FullIntegrals{
  double mass;
  double3 momentum;
  double Energy, kinEn, Enstropy, Entropy;
};
__managed__ FullIntegrals TotMoments;
__global__ void total_moments( FullIntegrals& totM );

void calcStep(int REV=1){
  cuTimer calct;
  parsHost.iStep++;
  std::vector<double> timings;
  int Ntiles=0;
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );
  calcLBM(parsHost.iStep, timings);
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );
  double phys_time=parsHost.iStep;
  double calc_time = calct.gettime();
  printf("Step %6d (physical time %6.3f ms) | Performance: %.2f ms (%.2f MLU/sec) | timings: ", 
      parsHost.iStep ,phys_time, calc_time,
      Nx*Ny*Nz/calc_time*1e-3     );
  for(auto tmg: timings) printf("%.2f ",tmg);
  printf("\n");
  
  if(parsHost.iStep%1==0) {
  memset( &TotMoments, 0, sizeof(FullIntegrals) );
  total_moments<<<dim3(Nx,Ny),Nz>>>(TotMoments); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  printf("Total Conservations: Mass %.15f Momentum( %.12f %.12f %.12f ), M2 %.12f\n",
                  TotMoments.mass, TotMoments.momentum.x, TotMoments.momentum.y, TotMoments.momentum.z, TotMoments.Energy );
  printf("Total Characteristics: KineticEnergy: %.15f Enstropy: %.15f Entropy: %.15f\n",
                  TotMoments.kinEn, TotMoments.Enstropy, TotMoments.Entropy );
  }
}

template<int n> struct KerRunner {
  static void run() { streaming_collision <n> <<<Nx*Ny*Nz,LBMconsts::Qn>>>(0); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() ); }
};
 
void calcLBM(int it, std::vector<double>& timings){
  cuTimer t0; double new_time=0, prev_time=0;

  cudaDeviceProp prop;
  CHECK_ERROR( cudaGetDeviceProperties( &prop, 0) );
  const int MaxBlocksPerSM=32;
  if(it==1) {
    printf("GPU SM count %3d \n", prop.multiProcessorCount);
    //CHECK_ERROR( cudaDeviceSetLimit(cudaLimitStackSize, 256*1024) );
    CHECK_ERROR( cudaDeviceSetLimit( cudaLimitMallocHeapSize, MaxBlocksPerSM*prop.multiProcessorCount*sizeof(MomentsMatrix) ) );
  }

  TemplateSwitcher<5, KerRunner<5> >::run( PPhost.RegOrder );
  //streaming_collision<<<Nx*Ny*Nz,1>>>(0); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  //using namespace CUDAstep;
  /*for(int ibn=0; ibn<Nx*Ny*Nz; ibn+= MaxBlocksPerSM*prop.multiProcessorCount) {
    printf("step %5d progress = %6d/%6d\r", it, ibn, Nx*Ny*Nz); fflush(stdout);
    streaming_collision<<<MaxBlocksPerSM*prop.multiProcessorCount,1>>>(ibn);
    //streaming_collision<<<1,1>>>(ibn);
    cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  }*/
  parsHost.data.swap_ptrs();
  copy2dev( parsHost, pars );
  debug_print(); timings.push_back( t0.getlaptime() );

  debug_print(); timings.push_back( t0.getlaptime() );

}

__global__ void total_moments( FullIntegrals& totMom ){
  const int ix=blockIdx.x;
  const int iy=blockIdx.y;
  const int iz=threadIdx.x;
  const int index = ix+iy*Nx+iz*Nx*Ny;
  using namespace LBMconsts;

  Cell cell = pars.data.get_cell(0, ix,iy,iz);
  const ftype rho=cell.rho;
  const ftype3 vel = cell.vel;
  const ftype T = cell.T;
  ftype entropy = 0; for(int iq=0;iq<LBMconsts::Qn;iq++) entropy+= cell.f[iq]*log(cell.f[iq]/w[iq]);
  Cell cellMx = pars.data.get_cell(0, (ix-1+Nx)%Nx, iy, iz);
  Cell cellPx = pars.data.get_cell(0, (ix+1   )%Nx, iy, iz);
  Cell cellMy = pars.data.get_cell(0, ix, (iy-1+Ny)%Ny, iz);
  Cell cellPy = pars.data.get_cell(0, ix, (iy+1   )%Ny, iz);
  Cell cellMz = pars.data.get_cell(0, ix, iy, (iz-1+Nz)%Nz);
  Cell cellPz = pars.data.get_cell(0, ix, iy, (iz+1   )%Nz);
  ftype3 vorticity;
  vorticity.x = 0.5*(cellPy.vel.z-cellMy.vel.z) - 0.5*(cellPz.vel.y-cellMz.vel.y);
  vorticity.y = 0.5*(cellPz.vel.x-cellMz.vel.x) - 0.5*(cellPx.vel.z-cellMx.vel.z);
  vorticity.z = 0.5*(cellPx.vel.y-cellMx.vel.y) - 0.5*(cellPy.vel.x-cellMy.vel.x);

  atomicAdd(&totMom.mass      , rho);
  atomicAdd(&totMom.momentum.x, rho*vel.x );
  atomicAdd(&totMom.momentum.y, rho*vel.y );
  atomicAdd(&totMom.momentum.z, rho*vel.z );
  atomicAdd(&totMom.Energy    , rho*T*DIM/2 + rho*dot(vel,vel)/2 );
  atomicAdd(&totMom.kinEn     , rho*dot(vel,vel)/2 );
  atomicAdd(&totMom.Enstropy  , rho*dot(vorticity,vorticity)/2 );
  atomicAdd(&totMom.Entropy   , entropy );
}

inline void debug_print(){
   return;
}
