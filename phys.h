#ifndef _PHYS_H
#define _PHYS_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <string>

#include "consts.h"

// solution from stackowerflow
static void _mkdir(const char *dir) {
  char tmp[PATH_MAX];
  char *p = NULL;
  size_t len;
  snprintf(tmp, sizeof(tmp),"%s",dir);
  len = strlen(tmp);
  if(tmp[len - 1] == '/') tmp[len - 1] = 0;
  for(p = tmp + 1; *p; p++) if(*p == '/') {
    *p = 0;
    mkdir(tmp, S_IRWXU);
    *p = '/';
  }
  mkdir(tmp, S_IRWXU);
}

struct InitPars{
  ftype uDragX,uDragY,u0,beta0,r0,T0;
  int planeTVG;
};

// Do not forget to copy parameters structure to GPU
struct PhysPars{
  ftype tau,dtau;
  ftype rho0;
  double phys_time;
  char* drop_dir;
  ftype dr,dx,dy,dz,dt;
  ftype RhoUnitConv, VelUnitConv, ViscUnitConv, TempUnitConv;
  //ftype viscosity;
  ftype visc_atT;
 
  int StepIterPeriod, MaxSteps;
  int stencilInterpWidth;
  int stencilFixed;
  int RegOrder;
  InitPars initial;

  void set_drop_dir(std::string dir) {
    //drop_dir=(std::string*)realloc(drop_dir,sizeof(std::string)); *drop_dir=dir;
    drop_dir = new char[dir.length() + 1];
    strcpy(drop_dir, dir.c_str());
    struct stat st = {0};
    if(stat(dir.c_str(), &st)==-1) _mkdir(dir.c_str());
  }
  void setupUnits();
  void setDefault(){
    RhoUnitConv=1;
    VelUnitConv= 1;
    ViscUnitConv= 1;
    TempUnitConv= 1;
    rho0=1;
    //viscosity = 0.01;
    visc_atT = 0.01;
    tau=0.5+visc_atT;
    dtau=1./tau;
    stencilInterpWidth=3;
    stencilFixed=0;
    RegOrder=-1;
  }
  void MallocData();
  void setCell(int val, int x,int y);
};
extern PhysPars PPhost;
#if not defined SWIG and not defined IGNORE_CUDA_EXTENSIONS
extern __constant__ PhysPars PPdev;
#endif
#endif //_PHYS_H
