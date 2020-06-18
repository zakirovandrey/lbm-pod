#ifndef IM3D_HPP
#define IM3D_HPP
#ifndef CUDA_BASIC_PARS
#define CUDA_BASIC_PARS
const int NW=32;//, SW=32, NT=SW*NW;//число тредов cuda в warp-е и block-е
#endif//CUDA_BASIC_PARS

#include "Arr3Dpars.hpp"
#include <string>

//const int Npx=1;//размер ячейки в пикселях
//NT не больше ограничения на число потоков в блоке (1024), также NT*Nregs/Thread для всех cuda-методов должно быть не больше числа регистров в SM (32K/64K) (-Xptxas="-v" как раз и выводит Nregs/Thread, пока что не более 30) 
//unsigned const int M0 = 0x55555555;
//unsigned const int M1 = 0xAAAAAAAA;

struct im3D_pars4save {
  float tstep, density, opacity, wellwidth;
  bool draw_mesh_flag, draw_box_flag, draw_fg_flag, dissect_box_flag, filterMode_flag, draw_sec_xyz_flag, draw_wells_flag;//, grad_mode, surf_mode;
  int mode3D, ix0,iy0,iz0;
  float viewRotation[3];
  float viewTranslation[3];
  float BoxFactor[3], RotPoint[3];
  float MeshBox[3], MeshShift[3], box_shrink[3];
  float Dzoom[3], Dadd[3];
  int Dshrink[3];
  float bkgr_col[3], mesh_col[3], box_col[3];
  float base[3], step[3];
  float Dmesh;
  char drop_dir[60];
  int draw_bmp4backgrownd;
  bool contour_flag;
  float contour_width;
  float cntr_levels[10];
  int cntr_num;
  int Narr[3];
  int2 ld_sz;
  //char ballast[3*128-288];//всего 128 байт
  int init_from_command_line(char** argv);
  void print_command_line_help();
  const char* command_line_help_string();
  void load_from_file(const char* fn);
};

struct im3D_pars: public im3D_pars4save {
  int Nx, Ny, Nz;//размер массива
  int bNx, bNy;//размер картинки
  int pal_sh, secType, secZsh, secYsh, secXsh;//, xyz_sh;
  unsigned int bits_sh;
  float x_zoom, y_zoom, z_zoom;
  float fMin, fMax;
  float2 eyePoint;
  const char* fName,* dfName;
  float randR,* randArr;
  float viewRotationTmp[3];
  float2 wlRange,cameraRange;
  bool checkNcopy(Arr3D_pars& arr) {
    if(Nx != arr.Nx*sizeof(float)/sizeof(floatT4im) || Ny != arr.Ny || Nz != arr.Nz) return false;
    fName = arr.fName;
    dfName = arr.dfName;
    fMin = arr.fMin;
    fMax = arr.fMax;
    return true;
  }
  void reset() {
    for(int i=0; i<3; i++) {
      viewRotation[i] = 0.0;
      viewRotationTmp[i] = 0.0;
      viewTranslation[i] = 0.0;
      BoxFactor[i] = 1.0;
      MeshBox[i] = 100.0;
      MeshShift[i] = 0.0;
      RotPoint[i] = 0.5;
      Dzoom[i] = 1.0;
      Dadd[i] = 0.0;
      Dshrink[i] = 1;
      Narr[i] = 0;
      step[i] = 1.0;
      base[i] = 0.0;
      bkgr_col[i] = 0.1;
      mesh_col[i] = 0.8;
      box_col[i] = 1.0;
      box_shrink[i] = 1.0;
    }
    mesh_col[2] = 0.2;
    viewTranslation[2] =-4.0;
    Dmesh = 2.0;
    tstep = 2.0; density = 0.5; opacity = 0.95; wellwidth=4;
    draw_mesh_flag=true; draw_wells_flag=true; draw_fg_flag=false; draw_box_flag=true; dissect_box_flag=false; filterMode_flag=true; draw_sec_xyz_flag=false;
    mode3D=0;//grad_mode=surf_mode=false;
    cntr_num = 0;
    contour_flag=false;
    contour_width = 0.2;
    ld_sz = make_int2(0,0);
    draw_bmp4backgrownd = 0;
    drop_dir[0] = '.'; drop_dir[1] = 0;
    ix0 = Nx/2; iy0 = Ny/2; iz0 = Nz/2;
    randR = 0.25;
    wlRange=make_float2(500,2500);
    cameraRange=make_float2(0,0);
  }
  void reset(Arr3D_pars& arr) {
    fName = arr.fName;
    dfName = arr.dfName;
    Nx = Narr[0]>0?Narr[0]:arr.Nx*sizeof(float)/sizeof(floatT4im);
    Ny = Narr[1]>0?Narr[1]:arr.Ny;
    Nz = Narr[2]>0?Narr[2]:arr.Nz;
    fMin = arr.fMin; fMax = arr.fMax;
    for(int i=0; i<3; i++) {
      if(Dzoom[i] <= 0.0) Dzoom[i] = 1.0;
      if(Dadd[i] <= 0.0) Dadd[i] = 0.0;
    }
    x_zoom=1.0/Dzoom[0]; y_zoom=1.0/Dzoom[1]; z_zoom=1.0/Dzoom[2];
    int NxZ=(1+Dadd[0])*Nx/x_zoom, NyZ=(1+Dadd[1])*Ny/y_zoom, NzZ=(1+Dadd[2])*Nz/z_zoom;
    int addW=NzZ, addH=NzZ;
    if(NzZ<=NxZ && NzZ<=NyZ) secType=0;
    else if(NyZ<=NxZ) { secType=1; addW = NyZ; }
    else { secType=2; addH = NxZ; }
    int w=NxZ+addW, x_gap=w%NW?(NW-w%NW):0; bNx = w+x_gap;
    int h=NyZ+addH, y_gap=h%NW?(NW-h%NW):0; bNy = h+y_gap;
    pal_sh = bNx*bNy;
    secZsh = bNx*(bNy-NyZ);
    secYsh = (secType<2)?0:bNx-NzZ;
    secXsh = bNx+((secType==1)?-NyZ:secZsh-NzZ);
    bNy+=20;
    ix0 = Nx/2; iy0 = Ny/2; iz0 = Nz/2;
    eyePoint.x = 0.5*bNx; eyePoint.y = 0.5*bNy;
  }
  void reset0(int x, int y);
  void shift0(int x, int y, int x1, int y1);
  void print_help();
  bool key_func(unsigned char key, int x, int y);
  bool special_func(unsigned char key, int x, int y);
  void recalc_func();
  void mouse_func(int button, int state, int x, int y);
  void motion_func(int x, int y);
  char* reset_title();
  void clear4exit();
  void init3D(Arr3D_pars& arr);
  void initCuda(Arr3D_pars& arr);
  void initCuda_surf(Arr3D_pars& arr, size_t sh=0);
  void initTex();
  void initTex_surf();
  void recalc_im3D();
  void recalc_sec_im3D();
  void recalc3D_im3D();
  void save_bmp4backgrownd();
  int init_from_command_line(char** argv);
  std::string getfName();
  std::string getfDropName(const char* ext, int it=-1);
  bool save_png(int it=-1);
  bool save_gp(int it=-1);
  std::string save_section(int it=-1);
  void plot_section();
};
struct any_idle_func_struct {
  virtual void step() {}
};
struct idle_func_calc: public any_idle_func_struct {
  float t;
  idle_func_calc(): t(0.0) {}
  void step();
};
#endif//IM3D_HPP
