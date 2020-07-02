//#include <cuda_fp16.h>
#include "cuda_math.h"
//#include "cuda_math_double.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <error.h>
#include "structs.cuh"
#include "init.h"

#include "im2D.h"
#include "im3D.hpp"

#include "phys.h"

int type_diag_flag=0;
extern im3D_pars im3DHost;
AllParamsHost parsHost;
__constant__ AllParams pars;

PhysPars PPhost;
__constant__ PhysPars PPdev;
void PhysPars::MallocData(){};
void PhysPars::setCell(int val, int x,int y){};

void AllParamsHost::set(){
  //PPhost.setUnits();
  sprintf(im3DHost.drop_dir, "%s", PPhost.drop_dir);

  if(PPhost.stencilFixed && PPhost.stencilInterpWidth%2) printf("Warning: stencil Width is Odd and stincil position is Fixed. Is it ok?\n");

  PPhost.setupUnits();

  //if(Nx%(1<<MaxLevel)!=0) error(1,1,"Error: Nx must be dividable by %d\n", 1<<MaxLevel);
}

int print_help();
void launch_im3D(int argc, char** argv);
bool interactive=true, test_only=false;
void reset(im3D_pars* p=0);
void init();
void simple_drop();
void calcStep(int REV=1);
int _main(int argc, char** argv) {
  int Ndevs=0; CHECK_ERROR( cudaGetDeviceCount(&Ndevs) ); CudaDevs=Ndevs;
  CudaDevs=1;
  ::reset();
  argv ++; argc --;
  im3DHost.reset();
  while(argc>0 && strncmp(*argv,"--",2)==0) {
    int pp=1;
    if(strcmp(*argv,"--test")==0) test_only = true;
    else if(strcmp(*argv,"--batch")==0) interactive = false;
    else pp = im3DHost.init_from_command_line(argv);
    if(pp<=0) return print_help();
    else if(pp==1) printf("par: %s; \n", argv[0]);
    else if(pp==2) printf("par: %s; vals: %s\n", argv[0], argv[1]);
    argv += pp; argc -= pp;
  };
  if(test_only) printf("No GL\n");
  else printf("With GL\n");
  im2D.get_device(3,0);
  type_diag_flag = 1;
try {
  if(type_diag_flag>=1) printf("Настройка опций визуализации по умолчанию\n");
  cudaTimer tm; tm.start();
  parsHost.set();
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  init();
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );

  if(test_only) {
    while(parsHost.iStep<parsHost.StepsMax || parsHost.StepsMax<=0 ) {
      tm.start();
      calcStep();
      double tCpu=tm.stop();
//      fprintf(parsHost.fLog,"run time: %.2f msec, %.2f MLU/sec\n", tCpu, 1.e-6*Nx*Ny*Nz/tCpu);
      printf("run time: %.2f msec, %.2f MLU/sec\n", tCpu, 1.e-6*Nx*Ny*Nz/tCpu);
    }
    return 0;
  } else{
    launch_im3D(argc,argv);
  }
} catch(...) {
  printf("Возникла какая-то ошибка.\n");
}
  parsHost.clear();
  return -1;
}
int main(int argc, char** argv) {
  PPhost.setDefault();
  return _main(argc,argv);
}
int run(int argc, char** argv) {
  return _main(argc,argv);
}
