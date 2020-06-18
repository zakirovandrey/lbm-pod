#include "cuda_math.h"
#include <stdio.h>
#include <stdlib.h>
#include "structs.cuh"

#include "im2D.h"
#include "im3D.hpp"

#include "LBMconsts.cuh"

#include "phys.h"

im3D_pars im3DHost;
void calcStep(int REV=1);

const char* FuncStr[] = {
  "rho", "Niter", "Vx", "Vy", "Vz", "T"
};

__global__ void __launch_bounds__(Nz) draw(float* buf) {
  int iz=threadIdx.x;
  int ix=blockIdx.x;
  int iy=blockIdx.y;

  float* pbuf=&buf[ix+gridDim.x*(iy+gridDim.y*iz)];
  Cell cell = pars.data.get_cell(0, ix,iy,iz);
  ftype rho=cell.rho;
  ftype3 vel = cell.vel;
  
  switch(pars.nFunc) {
      case 0 : *pbuf=float(rho); break;
      case 1 : *pbuf=float(cell.Niter); break;
      case 2 : *pbuf=float(vel.x); break;
      case 3 : *pbuf=float(vel.y); break;
      case 4 : *pbuf=float(vel.z); break;
      case 5 : *pbuf=float(cell.T); break;
      default: break;
  }
}
void draw_all(){
  CHECK_ERROR( cudaMemset(parsHost.arr4im.Arr3Dbuf,0,((long long int)parsHost.arr4im.Nx)*parsHost.arr4im.Ny*parsHost.arr4im.Nz*sizeof(float)) );
  draw<<<dim3(parsHost.arr4im.Nx,parsHost.arr4im.Ny),parsHost.arr4im.Nz>>>(parsHost.arr4im.Arr3Dbuf);
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  im3DHost.initCuda(parsHost.arr4im);
}

void idle_func_calc::step() {
  for(int ii=0;ii<PPhost.StepIterPeriod;ii++) {
    calcStep(); 
    t++;
  }
  draw_all();
  im3DHost.save_png(parsHost.iStep/PPhost.StepIterPeriod);
  recalc_at_once=true;
}

static void key_func(unsigned char key, int x, int y) {
  if(type_diag_flag>=2) printf("keyN=%d, coors=(%d,%d)\n", key, x, y);
  if(key == 'h') {
    printf("\
======= Управление:\n\
  <¦>  \tИзменение функции для визуализации: Values_level543210¦isBnd_level543210¦AMR_LEVEL\n\
«Enter»\tПересчёт одного шага\n\
   b   \tвключает пересчёт в динамике (см. «Управление динамикой»)\n\
"); im3DHost.print_help();
    return;
  }
  ftype t0;
  switch(key) {
  //case '>': if(parsHost.nFunc<parsHost.MaxFunc) parsHost.nFunc++; break;
  //case '<': if(parsHost.nFunc>0) parsHost.nFunc--; break;
  case '>': parsHost.nFunc = (parsHost.nFunc+1)%parsHost.MaxFunc; break;
  case '<': parsHost.nFunc = (parsHost.nFunc+parsHost.MaxFunc-1)%parsHost.MaxFunc; break;
  case 13: for(int ii=0;ii<PPhost.StepIterPeriod;ii++) calcStep(+1); /*im3DHost.save_png(parsHost.iStep/PPhost.StepIterPeriod)*/; break;
  case 8 : for(int ii=0;ii<PPhost.StepIterPeriod;ii++) calcStep(-1); /*im3DHost.save_png(parsHost.iStep/PPhost.StepIterPeriod)*/; break;
  default: if(!im3DHost.key_func(key, x, y)) {
  if(type_diag_flag>=0) printf("По клавише %d в позиции (%d,%d) нет никакого действия\n", key, x, y);
  } return;
  }
  copy2dev( parsHost, pars );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  draw_all();
  recalc_at_once=true;
}
static void draw_func() {
  im3DHost.fName = FuncStr[parsHost.nFunc]; 
  im2D.draw(im3DHost.reset_title()); 
}

static void idle_func() { im3DHost.recalc_func(); }
static void mouse_func(int button, int state, int x, int y) { im3DHost.mouse_func(button, state, x, y); }
static void motion_func(int x, int y) { im3DHost.motion_func(x, y); }
static void special_func(int key, int x, int y) { 
  im3DHost.special_func(key, x, y);
  if(key == GLUT_KEY_F2) {
//    parsHost.drawArrows^=1;
    copy2dev( parsHost, pars ); draw_all();
    recalc_at_once=true;
  }
}

int print_help() {
  printf("help | using in test|batch mode:\n ./lbm [--help|--test|--batch]\n");
  printf("using in interactive mode:\n ./lbm %s\n", im3DHost.command_line_help_string());
  im3DHost.print_command_line_help();
  return 0;
}
void read_float3(float* v, char* str);
float read_float(char* str);

void launch_im3D(int argc, char** argv){
  parsHost.nFunc = 1; parsHost.MaxFunc = sizeof(FuncStr)/sizeof(char*);
    
  cudaTimer tm; tm.start();
  parsHost.reset_im();
  im3DHost.reset(parsHost.arr4im);
  copy2dev( parsHost, pars );
  copy2dev( PPhost, PPdev );
  im2D.get_device(3,0);
  im2D.init_image(argc,argv, im3DHost.bNx, im3DHost.bNy, "im3D");
  im3DHost.init3D(parsHost.arr4im); im3DHost.iz0=parsHost.arr4im.Nx-1; im3DHost.key_func('b',0,0);
  im3DHost.initCuda(parsHost.arr4im);
  draw_all();

  if(type_diag_flag>=1) printf("Настройка GLUT и запуск интерфейса\n");
  glutIdleFunc(idle_func);
  glutKeyboardFunc(key_func);
  glutMouseFunc(mouse_func);
  glutMotionFunc(motion_func);
  glutDisplayFunc(draw_func);
  glutSpecialFunc(special_func);
  if(type_diag_flag>=0) printf("Init cuda device: %.1f msec\n", tm.stop());
  glutMainLoop();
}

float get_val_from_arr3D(int ix, int iy, int iz) {
  Arr3D_pars& arr=parsHost.arr4im;
  if(arr.inCPUmem) return arr.Arr3Dbuf[arr.get_ind(ix,iy,iz)];
  float res=0.0;
  if(arr.inGPUmem) CHECK_ERROR(cudaMemcpy(&res, arr.get_ptr(ix,iy,iz), sizeof(float), cudaMemcpyDeviceToHost));
  return res;
}

