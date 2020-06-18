//im3D считывает и визуализирует трёхмерные поля,
//получение исходного кода: bzr checkout bzr+ssh://photon/Save/BZR-for-all/lev/im3D
//автор: Вадим Левченко VadimLevchenko@mail.ru
// запуск: ./im3D <имя-файла-массива> [<имя-файла-массива> ...]
//целевой размер массива от 100 до 1500 элементов по каждой координате
//предполагается, что файлы массивов записаны в формате массивов aivlib-а или drp

#include "cuda_math.h"
#include "fpal.h"
#include "im2D.h"
#include "im3D.hpp"

image2D im2D;
image_pars imHost; __constant__ image_pars im;
__constant__ im3D_pars im3D;

float runTime=0.0, SmoothFPS=0.0;
bool recalc_at_once=true, recalc_always=false, save_anim_flag=false, draw_edges_flag=false;
int anim_acc=0, render_type=3;
texture<floatT4im, cudaTextureType3D> data3D_tex;
cudaArray* data3D_texArray=0;
texture<short2, cudaTextureType3D> data3Dsurf_tex;
cudaArray* data3Dsurf_texArray=0;
const char* optfName="im3DI.opt";//Имя файла для сохранения опций визуализации
FILE* gpPipe=0;
int sec1Daxis=0;

//#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
//#include <time.h>

#include <malloc.h>

char WinTitle[1024], addTitleStr[5]; int TitleStrInd=0;
const char* baseTitleStr="2"; int baseTitleFlag=1;
int optfid=-1; int im3DIopt_shift=0;

void im3D_pars4save::load_from_file(const char* fn) {
  int optfid=open(fn, O_RDWR|O_CREAT, 0644);
  if(optfid<0) { printf("Не могу открыть файл %s, загрузка наборов опций визуализации невозможна\n", fn); return; }
  int sz,rs;
  rs = read(optfid, &sz, sizeof(sz));
  if(sz<=0 || sz>sizeof(fpal_pars)) { printf("Illegal Drop format\n"); return; }
  rs=read(optfid, &imHost, sz); printf("Load %dB fpal of %ldB", rs, sizeof(fpal_pars));
  rs = read(optfid, &sz, sizeof(sz));
  if(sz<=0 || sz>sizeof(im3D_pars4save)) { printf("Illegal Drop format\n"); return; }
  rs=read(optfid, this, sz); printf(" & %dB im3D of %ldB\n", rs, sizeof(im3D_pars4save));
  close(optfid);
}
char* im3D_pars::reset_title() {
  char* pTitle=WinTitle, TitleStr[20];
  if(baseTitleFlag%3>0) strcpy(TitleStr,baseTitleStr);
  strncpy(TitleStr+strlen(TitleStr),addTitleStr,4);
  if(baseTitleFlag%3==1 && fName) { sprintf(pTitle, "%s ", fName); pTitle += strlen(pTitle); }
  if(baseTitleFlag%3==2 && dfName) { sprintf(pTitle, "%s ", dfName); pTitle += strlen(pTitle); }
  if(strpbrk(TitleStr,"23")) { sprintf(pTitle, "(%dx%dx%d)", Nx,Ny,Nz); pTitle += strlen(pTitle); }
  if(strpbrk(TitleStr,"xyzXYZ")) { sprintf(pTitle, "/(%dx%dx%d)", ix0,iy0,iz0); pTitle += strlen(pTitle); }
  if(strpbrk(TitleStr,"aA\001\023=-+_06789")) { sprintf(pTitle, " %g<f<%g", imHost.fmin,imHost.fmax); pTitle += strlen(pTitle); }
  if(strpbrk(TitleStr,"pP[]|?{}")) { sprintf(pTitle, " pal[%d]:(%g)^%g*%g*%g|%g;", imHost.palID, imHost.pscale, imHost.gamma_pal, imHost.brightness_coff, imHost.max_rgb, imHost.base_val); pTitle += strlen(pTitle); }
  if(strpbrk(TitleStr,"dDjJmM")) { sprintf(pTitle, " D/J/M:%g/%g/%g;", density, opacity, tstep); pTitle += strlen(pTitle); }
#ifdef CALC_TIME_DIAG
extern float calcTime, calcPerf; extern int TimeStep;
  if(strpbrk(TitleStr,"bG")) { sprintf(pTitle, " calc: %.2f sec, %.2fG cells/sec; %d steps;", 1e-3*calcTime, calcPerf, TimeStep); pTitle += strlen(pTitle); }
#endif
  if(strpbrk(TitleStr,"tT\20")) { sprintf(pTitle, " transp: %s,%d", imHost.transparency_discrete_flag?"discr":"mode",imHost.transparency_mode); pTitle += strlen(pTitle); }
  //if(strpbrk(TitleStr,"gG")) { sprintf(pTitle, " ", ); pTitle += strlen(pTitle); }
  //sprintf(WinTitle, " %.1f fps", , ,recalc_always?SmoothFPS:1000./runTime);
  //printf(WinTitle, " render: %.1f fps", , recalc_always?SmoothFPS:1000./runTime);
  return WinTitle;
}
struct RotMatr {
  double v[3][3];
  RotMatr(int c, double phi) {
    phi *= M_PI/180;
    int cp=(c+1)%3, cm=(c+2)%3;
    for(int i=0; i<3; i++) v[i][c] = v[c][i] = 0,0;
    v[c][c] = 1.0;
    v[cp][cp] = v[cm][cm] = cos(phi);
    v[cm][cp] = sin(phi); v[cp][cm] =-v[cm][cp];
  }
  void operator *= (RotMatr& M) {
    double vo[3][3];
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) vo[i][j] = v[i][j];
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
      double vn=0.0;
      for(int k=0; k<3; k++) vn += M.v[i][k]*vo[k][j];
      v[i][j] = vn;
    }
  }
};

std::string im3D_pars::getfName() {
  char fN[]="image.__________";
  if(fName) strncpy(fN, fName, sizeof(fN)-1);//[sizeof(fN)-1] = 0;
  if(strrchr(fN,'.')) strrchr(fN,'.')[0] = 0;
  if(strrchr(fN,' ')) strrchr(fN,' ')[0] = 0;
  if(strrchr(fN,'/')) strrchr(fN,'/')[0] = '_';
  return std::string(fN);
}
std::string im3D_pars::getfDropName(const char* ext, int it) {
  char drop_name[1024];
  sprintf(drop_name, "%s/%s_%04d%s", drop_dir, getfName().c_str(), (it>=0)?it:imHost.nFrame, ext);
  return std::string(drop_name);
}

bool im3D_pars::save_png(int it) {
  im2D.out2png(getfDropName(".png", it).c_str());
  imHost.nFpng++; 
  return false;
}
__global__ void save_gp3D();
const int tileSz=16, tilesN=16;

bool im3D_pars::save_gp(int it) {
  std::string png_name=getfDropName(".png", it);
  std::string gp_name=getfDropName(".gp", it);
  im2D.out2png(png_name.c_str());
  //sprintf( gp_name, "a.gp", fName, imHost.nFrame);
  FILE* gp=fopen(gp_name.c_str(), "w"),* old_stdout=stdout;
  fprintf(gp, "unset key\n");
  fprintf(gp, "unset border\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "set x2tics border\n");
  fprintf(gp, "set x2range [%g:%g]\n", imHost.fmin, imHost.fmax);
  fprintf(gp, "unset ytics\n");
  //fprintf(gp, "load \"labels.gp\"\n");
  //printf("viewRotation: %g, %g\n", viewRotation[0], viewRotation[1]);
  //printf("viewTranslation: %g, %g, %g\n", viewTranslation[0], viewTranslation[1], viewTranslation[2]);
  if(render_type==3) {
    const int Sgp=(tilesN-1)*tileSz;
    stdout = gp;
    if(CHECK_ERROR(cudaDeviceSynchronize())) throw(-1);
    save_gp3D <<<dim3((im2D.Nx+Sgp-1)/Sgp,(im2D.Ny+Sgp-1)/Sgp),dim3(tilesN,tilesN)>>>();
    if(CHECK_ERROR(cudaDeviceSynchronize())) throw(-1);
    stdout = old_stdout;
  }
  fprintf(gp, "plot[0:%g][0:%g] \"%s\" binary filetype=png dx=1 dy=1 with rgbimage\n", float(bNx), float(bNy), png_name.c_str());
  fprintf(gp, "pause -1\n");
  fclose(gp);
  if(type_diag_flag>=0) printf("Зарамочное оформление сохранено в %s\n", gp_name.c_str());
  return false;
}

floatT4im get_val_from_arr3D(int ix, int iy, int iz);
void reset(im3D_pars* p=0);
#if DATA_VECTOR_SZ==1
std::string im3D_pars::save_section(int it) {
  printf("f(%d,%d,%d) = %g\n", ix0, iy0, iz0, get_val_from_arr3D(ix0, iy0, iz0));
  std::string dat_name=getfDropName(".dat",it);
  FILE* dat=fopen(dat_name.c_str(), "w");
  for(int i=0; i<Nx; i++) fprintf(dat, "%d %g\n", i, get_val_from_arr3D(i, iy0, iz0));
  fprintf(dat, "\n\n");
  for(int i=0; i<Ny; i++) fprintf(dat, "%d %g\n", i, get_val_from_arr3D(ix0, i, iz0));
  fprintf(dat, "\n\n");
  for(int i=0; i<Nz; i++) fprintf(dat, "%d %g\n", i, get_val_from_arr3D(ix0, iy0, i));
  fclose(dat);
  return dat_name;
}
void im3D_pars::plot_section() {
  const char* re=gpPipe?"re":"";
  if(gpPipe==NULL) gpPipe = popen("gnuplot", "w");
  int sec[]={ix0,iy0,iz0,ix0,iy0};
  if(sec1Daxis<3) fprintf(gpPipe, "set style data l;\n%splot '%s' i %d t '%c:(%d,%d)'\n", re, save_section().c_str(), sec1Daxis, "xyz"[sec1Daxis], sec[sec1Daxis+1], sec[sec1Daxis+2]);
  else fprintf(gpPipe, "set style data l;\n%splot '%s' u ($1-%d):2 i 0 t '(ix-%d)', '' u ($1-%d):2 i 1 t '(iy-%d)', '' u ($1-%d):2 i 2 t '(iz-%d)'\n", re, save_section().c_str(), ix0,ix0,iy0,iy0,iz0,iz0);
  fflush(gpPipe);
}
#elif DATA_VECTOR_SZ==3
std::string im3D_pars::save_section(int it) {
  floatT4im v=get_val_from_arr3D(ix0, iy0, iz0);
  printf("f(%d,%d,%d) = (%g,%g,%g,%g)\n", ix0, iy0, iz0, v.x, v.y, v.z, v.w);
  std::string dat_name=getfDropName(".dat",it);
  FILE* dat=fopen(dat_name.c_str(), "w");
  for(int i=0; i<Nx; i++) {
    floatT4im v=get_val_from_arr3D(i, iy0, iz0);
    fprintf(dat, "%d %g %g %g %g\n", i, v.x, v.y, v.z, v.w);
  }
  fprintf(dat, "\n\n");
  for(int i=0; i<Ny; i++) {
    floatT4im v=get_val_from_arr3D(ix0, i, iz0);
    fprintf(dat, "%d %g %g %g %g\n", i, v.x, v.y, v.z, v.w);
  }
  fprintf(dat, "\n\n");
  for(int i=0; i<Nz; i++) {
    floatT4im v=get_val_from_arr3D(ix0, iy0, i);
    fprintf(dat, "%d %g %g %g %g\n", i, v.x, v.y, v.z, v.w);
  }
  fclose(dat);
  return dat_name;
}
void im3D_pars::plot_section() {
  const char* re=gpPipe?"re":(render_type==3?"s":"");
  if(gpPipe==NULL) gpPipe = popen("gnuplot", "w");
  int sec[]={ix0,iy0,iz0,ix0,iy0};
  if(render_type==3) {
    if(sec1Daxis<3) fprintf(gpPipe, "set ticslevel 0; set style data lp;\n%splot '%s' u 2:3:4 i %d t '%c:(%d,%d)'\n", re, save_section().c_str(), sec1Daxis, "xyz"[sec1Daxis], sec[sec1Daxis+1], sec[sec1Daxis+2]);
    else fprintf(gpPipe, "set ticslevel 0; set style data lp;\n%splot '%s' u 2:3:4 i 0 t '(ix-%d)', '' u 2:3:4 i 1 t '(iy-%d)', '' u 2:3:4 i 2 t '(iz-%d)'\n", re, save_section().c_str(), ix0,iy0,iz0);
  } else if(render_type==2) {
    if(sec1Daxis<3) fprintf(gpPipe, "set style data l;\n%splot '%s' u 1:2 i %d t '[%d,%d].x', '' u 1:3 i %d t '.y', '' u 1:4 i %d t '.z'\n", re, save_section().c_str(), sec1Daxis, sec[sec1Daxis+1], sec[sec1Daxis+2], sec1Daxis, sec1Daxis);
    else fprintf(gpPipe, "set style data l;\n%splot '%s' u ($1-%d):2 i 0 t '(ix-%d).x', '' u ($1-%d):3 i 0 t '.y', '' u ($1-%d):4 i 0 t '.z', '' u ($1-%d):2 i 1 t '(iy-%d).x', '' u ($1-%d):3 i 1 t '.y', '' u ($1-%d):4 i 1 t '.z', '' u ($1-%d):2 i 2 t '(iz-%d).x', '' u ($1-%d):3 i 2 t '.y', '' u ($1-%d):4 i 2 t '.z'\n", re, save_section().c_str(), ix0,ix0,ix0,ix0,iy0,iy0,iy0,iy0,iz0,iz0,iz0,iz0);
  //else fprintf(gpPipe, "set style data l;\n%splot '%s' u ($1-%d):2 i 0 t '(ix-%d)', '' u ($1-%d):2 i 1 t '(iy-%d)', '' u ($1-%d):2 i 2 t '(iz-%d)'\n", re, save_section().c_str(), ix0,ix0,iy0,iy0,iz0,iz0);
  }
  fflush(gpPipe);
}
#endif

void im3D_pars::clear4exit() {
  im2D.clear();
  CHECK_ERROR(cudaFreeArray(data3D_texArray));
  CHECK_ERROR(cudaFreeArray(data3Dsurf_texArray));
  CHECK_ERROR(cudaFree(randArr));
}
void save_bmp4backgrownd();
any_idle_func_struct xyz_void,* xyz=&xyz_void;
struct idle_func_struct3D: public any_idle_func_struct {
  float* par, val;
  void set(float* _par, float _val) { par = _par; val = _val; }
  void step() { *par += val; }
} xyz3D;
struct idle_func_struct2D: public any_idle_func_struct {
  int* i0, N, di;
  void set(int* _i0, int _N, int _di) { i0=_i0; N=_N; di=_di; }
  void step() { *i0 += di; if(*i0<0) *i0=N-1; else if(*i0>=N) *i0=0; }
} xyz2D;
idle_func_calc icalc;
template<class Tflt>
struct idle_func_calcNdrop: public idle_func_calc {
  FILE* sensorsStr;
  int* sensors;
  int Nsensors;
  idle_func_calcNdrop(): Nsensors(0), sensors(0), sensorsStr(0) {}
  ~idle_func_calcNdrop() { delete sensors; }
  void add_sensor(int ix, int iy, int iz) {
    int* pi=sensors;
    for(int i=0; i<Nsensors; i++, pi+=3) if(pi[0] == ix && pi[1] == iy && pi[2] == iz)
      { printf("Сенсор (%d,%d,%d) уже задан. Вы делаете что-то не то!\n", ix, iy, iz); return; }
    Nsensors++;
    printf("Создаю новый сенсор в точке (%d,%d,%d), теперь их %d, файл <sensors.dat> будет очищен.\n", ix, iy, iz, Nsensors);
    if(sensors == 0) sensors = (int*)malloc(Nsensors*3*sizeof(int));
    else sensors = (int*)realloc(sensors, Nsensors*3*sizeof(int));
    pi = sensors+3*(Nsensors-1);
    pi[0] = ix; pi[1] = iy; pi[2] = iz;
    sensorsStr = fopen("sensors.dat", "w");
    pi=sensors;
    fprintf(sensorsStr, "#");
    for(int i=0; i<Nsensors; i++, pi+=3) fprintf(sensorsStr, "\t(%d,%d,%d)", pi[0],pi[1],pi[2]);
    fprintf(sensorsStr, "\n");
    fclose(sensorsStr);
  }
#if DATA_VECTOR_SZ==1
void step() {
  idle_func_calc::step();
  if(Nsensors==0) return;
  sensorsStr = fopen("sensors.dat", "a");
  fprintf(sensorsStr, "%g", t);
  int* pi=sensors;
  for(int i=0; i<Nsensors; i++, pi+=3) fprintf(sensorsStr, "\t%g", get_val_from_arr3D(pi[0], pi[1], pi[2]));
  fprintf(sensorsStr, "\n");
  fclose(sensorsStr);
}
void plot_sensors() {
  if(Nsensors==0) return;
  if(gpPipe==NULL) gpPipe = popen("gnuplot", "w");
  int* pi=sensors;
  fprintf(gpPipe, "set style data l;\nplot 'sensors.dat' u 1:2 t '%d,%d,%d'", pi[0],pi[1],pi[2]); pi+=3;
  for(int i=1; i<Nsensors; i++, pi+=3) fprintf(gpPipe, ", '' u 1:%d t '%d,%d,%d'", i+2, pi[0],pi[1],pi[2]);
  fprintf(gpPipe, "\n");
  fflush(gpPipe);
}
#elif DATA_VECTOR_SZ==3
void step() {
  idle_func_calc::step();
  if(Nsensors==0) return;
  sensorsStr = fopen("sensors.dat", "a");
  fprintf(sensorsStr, "%g", t);
  int* pi=sensors;
  for(int i=0; i<Nsensors; i++, pi+=3) {
    floatT4im v=get_val_from_arr3D(pi[0], pi[1], pi[2]);
    fprintf(sensorsStr, "\t%g\t%g\t%g", v.x, v.y, v.z);
  }
  fprintf(sensorsStr, "\n");
  fclose(sensorsStr);
}
void plot_sensors() {
  if(Nsensors==0) return;
  if(gpPipe==NULL) gpPipe = popen("gnuplot", "w");
  int* pi=sensors;
  if(render_type==3) {
    fprintf(gpPipe, "set style data lp; set ticslevel 0;\nsplot 'sensors.dat' u 2:3:4 t '%d,%d,%d'", pi[0],pi[1],pi[2]); pi+=3;
    for(int i=1; i<Nsensors; i++, pi+=3) fprintf(gpPipe, ", '' u %d:%d:%d t '%d,%d,%d'", 3*i+2,3*i+3,3*i+4, pi[0],pi[1],pi[2]);
  } else if(render_type==2) {
    fprintf(gpPipe, "set style data l;\nplot 'sensors.dat' u 1:2 t '[%d,%d,%d].x', '' u 1:3 t '[%d,%d,%d].y', '' u 1:4 t '[%d,%d,%d].z'", pi[0],pi[1],pi[2], pi[0],pi[1],pi[2], pi[0],pi[1],pi[2]); pi+=3;
    for(int i=1; i<Nsensors; i++, pi+=3) fprintf(gpPipe, ", '' u 1:%d t '[%d,%d,%d].x', '' u 1:%d t '[%d,%d,%d].y', '' u 1:%d t '[%d,%d,%d].z'", 3*i+2, pi[0],pi[1],pi[2],3*i+3, pi[0],pi[1],pi[2],3*i+4, pi[0],pi[1],pi[2]);
  }
  fprintf(gpPipe, "\n");
  fflush(gpPipe);
}
#endif
};
idle_func_calcNdrop<floatT4im> icalcNdrop;
//void add_sensor(int ix, int iy, int iz) { icalcNdrop.add_sensor(ix, iy, iz); }

#include<curand.h>
#include<curand_kernel.h>
__global__ void init_rand(curandState *states, float* randArr) {
  unsigned int tid = threadIdx.x + blockDim.x * blockIdx.x;
  curand_init(1234, tid, 0, &states[tid]);  //  Initialize CURAND
  randArr[tid] = 2.*M_PI*curand_uniform (&states[tid]);     // between 0 and 1
}
__device__ float get_float4lim(float v) { return v; }
__device__ float get_float4lim(float2 v) { return length(v); }
__device__ float get_float4lim(float4 v) { return length(v); }
__global__ void calc_limits3D(uint3 IB, uint3 IE, uint3 blkSz, uint3 Nthr, float2* fLims) {
  float2 fLim;
  IB+=blkSz*blockIdx*make_uint3(blockDim)+make_uint3(threadIdx.x/(Nthr.y*Nthr.z), (threadIdx.x/Nthr.z)%Nthr.y, threadIdx.x%Nthr.z);
  IE=min(IE,IB+blkSz);
  //if(threadIdx.x==0) printf("Blk %d from (%d,%d,%d) to (%d,%d,%d)\n",blockIdx.x+gridDim.x*(blockIdx.y+gridDim.y*blockIdx.z),IB.x,IB.y,IB.z, IE.x,IE.y,IE.z);
  fLim.x = fLim.y = get_float4lim(tex3D(data3D_tex, IB.x,IB.y,IB.z));
  for(int ix=IB.x; ix<IE.x; ix+=Nthr.x) for(int iy=IB.y; iy<IE.y; iy+=Nthr.y) for(int iz=IB.z; iz<IE.z; iz+=Nthr.z) {
    float v=get_float4lim(tex3D(data3D_tex, ix,iy,iz));
    if(v<fLim.x) fLim.x = v;
    if(v>fLim.y) fLim.y = v;
  }
  __shared__ float2 fLim_sh[512];
  fLim_sh[threadIdx.x] = fLim;
  __syncthreads();
  if(threadIdx.x >= warpSize) return;
  for(int i=threadIdx.x; i<blockDim.x; i+=warpSize) {
    float2 v=fLim_sh[i];
    if(v.x<fLim.x) fLim.x = v.x;
    if(v.y>fLim.y) fLim.y = v.y;
  }
  fLim_sh[threadIdx.x] = fLim;
  if(threadIdx.x>0) return;
  for(int i=0; i<warpSize; i++) {
    float2 v=fLim_sh[i];
    if(v.x<fLim.x) fLim.x = v.x;
    if(v.y>fLim.y) fLim.y = v.y;
  }
  fLims[blockIdx.x+gridDim.x*(blockIdx.y+gridDim.y*blockIdx.z)] = fLim;
  //printf("Lim (%d,%d,%d) %d => %g %g\n",blockIdx.x,blockIdx.y,blockIdx.z, blockIdx.x+gridDim.x*(blockIdx.y+gridDim.y*blockIdx.z), fLim.x,fLim.y);
}
float2 set_lim_from_tex(uint3 IB, uint3 N) {
  //if(N.x*N.y*N.z<512) { printf("Too small picture\n"); return make_float2(0.,1.); }
  int ind=0; uint3 Ns=N, Nthr;
  if(Ns.x<Ns.y) { ind += 3; int t=Ns.x; Ns.x=Ns.y; Ns.y=t; }
  if(Ns.y<Ns.z) { ind ++; int t=Ns.y; Ns.y=Ns.z; Ns.z=t; }
  if(Ns.x<Ns.y) { ind ++; int t=Ns.x; Ns.x=Ns.y; Ns.y=t; }
  for(Nthr.z=1; Nthr.z<8&&Nthr.z<Ns.z; Nthr.z*=2);
  for(Nthr.y=1; Nthr.y*Nthr.z<64&&Nthr.y<Ns.y; Nthr.y*=2);
  for(Nthr.x=1; Nthr.x*Nthr.y*Nthr.z<512&&Nthr.x<Ns.x; Nthr.x*=2);
  //printf("set Lim from tex: from (%d,%d,%d) size (%d,%d,%d); ind %d; Nthr: (%d,%d,%d)\n",IB.x,IB.y,IB.z, N.x,N.y,N.z, ind, Nthr.x,Nthr.y,Nthr.z);
  if(ind%3==2) { int t=Nthr.x; Nthr.x=Nthr.y; Nthr.y=t; }
  if(ind%3>=1) { int t=Nthr.y; Nthr.y=Nthr.z; Nthr.z=t; }
  if(ind  >=3) { int t=Nthr.x; Nthr.x=Nthr.y; Nthr.y=t; }

  uint3 Sblk=make_uint3(512), Nblk=(N+(Sblk-1))/Sblk;
  int NNblk=Nblk.x*Nblk.y*Nblk.z;
  float2 fLim,* fLims=0,* fLimsD=0;
  if(CHECK_ERROR(cudaMalloc((void**) &fLimsD, NNblk*sizeof(float2)))) throw(-1);
  //printf("Lim: %d*%d*%d => %d Blks, %d %d %d Thrs\n",Nblk.x,Nblk.y,Nblk.z,NNblk, Nthr.x,Nthr.y,Nthr.z);
  calc_limits3D<<<Nblk,Nthr.x*Nthr.y*Nthr.z>>>(IB, IB+N, Sblk, Nthr, fLimsD);
  fLims=new float2[NNblk];
  if(CHECK_ERROR(cudaMemcpy(fLims, fLimsD, NNblk*sizeof(float2), cudaMemcpyDeviceToHost))) throw(-1);
  CHECK_ERROR(cudaFree(fLimsD));
  fLim = *fLims;
  for(int i=1; i<NNblk; i++) {
    if(fLims[i].x<fLim.x) fLim.x = fLims[i].x;
    if(fLims[i].y>fLim.y) fLim.y = fLims[i].y;
  }
  delete fLims;
  return fLim;
}

int print_help();

void im3D_pars::print_help() {
  ::print_help();
  printf("\
======= Общее управление программой:\n\
 «ESC» \tВыход из программы\n\
  3¦2  \tпереключает рендеринг 3D¦2D в сечениях (%dD)\n\
   4   \tв режиме 3D переключает режим визуализации потенциал/градиентный режим/на поверхности\n\
<Enter¦BackSpace>\tПереход к следующему¦предыдущему массиву\n\
  w¦W  \tСохранение текущего набора опций визуализации в файл «%s»¦то же, но предыдущий набор не переписывается, можно сохранить произвольное число наборов последовательно\n\
  r¦R  \tЗагрузка ранее сохранённых наборов опций последовательно¦загрузка без перехода к следующему набору\n\
«Ctr-r»\tСброс параметров в значения по умолчанию\n\
  f¦F  \tПереход к началу¦концу файла сохранённых наборов опций\n\
  v¦V  \tУвеличение¦уменьшение уровня вывода диагностики (%d)\n\
«Ctr-v»\tПечатает диагностику, особено актуально, если заголовок окна не виден\n\
«Ctr-w»\tПереключает режим вывода в заголовок окна диагностики по умолчанию\n\
  s¦S  \tСохранение картинки в формате png|вместе с зарамочным оформлением в gnuplot\n\
   ~   \tВключает показ положения выделенной точки (x0,y0,z0), xyz при этом работают в режиме 2D\n\
 #¦$¦%% \tПереключение режима зарамочного оформления режима 3D: сетка¦рёбра бокса¦передний план\n\
   @   \tПереключает режим фона: однотонный/сохранённая картинка/2D сечения через выделенную точку/через сетку 3D\n\
   !   \tСохраняет картинку для фона\n\
«Ctr-z»\tУстанавливает координаты точки, относительно которой происходит вращение, в значение выделенной, при этом сдвигается сетка 3D\n\
  k¦K  \tУменьшение¦увеличение ширины линий контура\n\
«Ctr-k»\tВключает режим прорисовки линий контура (в 2D)\n\
  m¦M  \tУменьшение¦увеличение шага вдоль луча для соответствующего изменения точности (%g), ВНИМАНИЕ: при мелком шаге может очень медленно прорисовывать\n\
  e¦E  \tРазмазывание луча по горизонтали для соответствующего изменения муара (%g),\n\
  d¦D  \tУвеличение¦уменьшение плотности цвета при суммировании вдоль луча (%g)\n\
  j¦J  \tУменьшение¦увеличение порога цветовой плотности (%g)\n\
«Ctr-f»\tПереключает режим интерполяции в текстуре с режима по умолчанию (линейный в 3D/point в 2D)\n\
«Ctr-d»\tв 3D режиме отсекает часть массива\n\
«Ctr-L»\tЧитает параметры командной строки из текстового файла <add.opt>, формат: 1 параметр на строку, список значений без кавычек\n\
  a¦A  \tУстановка пределов палитры из пределов текущего массива ¦ из значений fMin..fMax\n\
«Ctr-a»\tУстановка значений fMin..fMax из текущих пределов палитры\n\
«Ctr-s»\tУстановка пределов палитры, используя пределы массива в сечении поперёк выбранной оси\n\
   1   \tпереключает (циклически, по xyz) ось, вдоль которой строится одномерный график в gnuplot  (%c)\n\
o¦«Ctr-o»\tВыводит в окно gnuplot сечение вдоль выбранной оси¦то же с перерисовкой\n\
  O¦Q  \tДля точки (x0,y0,z0): Печатает в терминале значение текущего поля и выводит в файл сечения вдоль лучей, проходящих через неё¦Добавляет сенсор\n\
q¦«Ctr-q»\tсохраняет значения сенсоров в файле sensors.dat¦выводит в окно gnuplot запись сенсоров\n\
======= Управление динамикой:\n\
  g¦G  \tОтключение¦включение постоянной перерисовки в цикле GLUT (%d)\n\
xyz¦XYZ\tВ режиме 2D, а также в 3D в режиме визуализации положения сечений: Увеличение¦уменьшение координат выделенной точки параллелепипеда данных (%d,%d,%d)\n\
xyz¦XYZ\tВ режиме 3D: Вращение вокруг осей x,y,z вперёд¦назад (%g,%g,%g)\n\
  u¦U  \tВ режиме 3D: Приближение¦удаление объекта (%g)\n\
======= Управление мышью (L¦R¦M --- левая¦правая¦средняя кнопки):\n\
   L   \tВ режиме 2D переустанавливает срезы, исходя из координат выбранной точки\n\
 L¦R¦M \tВ режиме 3D: вращение¦изменение масштаба¦сдвиг рисунка\n\
«Ctr-L»\tСдвиг картинки под курсором\n\
 В районе палитры (верхние 20 точек):\n\
 L¦R¦M \tустанавливает нижний¦верхний пределы¦центр палитры, исходя из x-координаты выбранной точки\n\
  L¦R  \tВ режиме «Ctl-t» (бинарной прозрачности) делает цвет прозрачным¦видимым\n\
", render_type, optfName, type_diag_flag, tstep, randR, density, opacity, "xyz"[sec1Daxis], recalc_always, ix0, iy0, iz0, viewRotation[0], viewRotation[1], viewRotation[2], viewTranslation[2]);
  imHost.print_help();
}
// normal          shift           Ctrl
//«DEL»
//`     67            %^&*      `1234567890 
//       i              I         e  yu   []
//         ;'         H   :"      d ghj  ;'\
//     n ,.          BN <>      zx  bnm,./ 
bool im3D_pars::key_func(unsigned char key, int x, int y) {
  recalc_at_once=true;
  size_t rs=0;
  if(key != addTitleStr[TitleStrInd]) addTitleStr[(TitleStrInd++)%4] = key;
  switch(key) {
  case 'A': imHost.set_lim(fMin, fMax); return true;
  case 'a': { float2 fLim=set_lim_from_tex(make_uint3(0,0,0), make_uint3(Nx,Ny,Nz)); imHost.set_lim(fLim.x, fLim.y); } return true;
//  case 'a': { float2 fLim=make_float2(-0.15,+0.15); imHost.set_lim(fLim.x, fLim.y); } return true;
  case 1: { fMin = imHost.fmin; fMax = imHost.fmax; } return true;
  case 19: { float2 fLim=make_float2(-1,1);
    switch(sec1Daxis) {
      case 0: fLim=set_lim_from_tex(make_uint3(ix0,0,0), make_uint3(1,Ny,Nz)); break;
      case 1: fLim=set_lim_from_tex(make_uint3(0,iy0,0), make_uint3(Nx,1,Nz)); break;
      case 2: fLim=set_lim_from_tex(make_uint3(0,0,iz0), make_uint3(Nx,Ny,1)); break;
    }; imHost.set_lim(fLim.x, fLim.y);
  } return true;
  case 18: ::reset(this); return true;
  case 'w': {
    printf("Drop %ldB fpal & %ldB im3D\n", sizeof(fpal_pars), sizeof(im3D_pars4save));
    if(optfid>=0 && im3DIopt_shift) rs=lseek(optfid,-im3DIopt_shift, SEEK_CUR);
  }
  case 'W': if(optfid>=0) {
    int sz=sizeof(fpal_pars); im3DIopt_shift = 0;
    rs=write(optfid, &sz, sizeof(sz)); im3DIopt_shift += rs;
    rs=write(optfid, &imHost, sz); im3DIopt_shift += rs;
    sz = sizeof(im3D_pars4save);
    rs=write(optfid, &sz, sizeof(sz)); im3DIopt_shift += rs;
    rs=write(optfid, this, sz); im3DIopt_shift += rs;
  } recalc_at_once=false; return true;
  case 'R': if(optfid>=0 && im3DIopt_shift) rs=lseek(optfid,-im3DIopt_shift, SEEK_CUR);
  case 'r': if(optfid>=0) {
    int sz=ld_sz.x;
    im3DIopt_shift = 0;
    if(sz==0) {
      rs = read(optfid, &sz, sizeof(sz));
      if(sz<=0 || sz>sizeof(fpal_pars)) { printf("Illegal Drop format\n"); return true; }
      im3DIopt_shift += rs;
    }
    rs=read(optfid, &imHost, sz); printf("Load %ldB fpal of %ldB", rs, sizeof(fpal_pars)); im3DIopt_shift += rs;
    sz=ld_sz.y;
    if(sz==0) {
      rs = read(optfid, &sz, sizeof(sz));
      if(sz<=0 || sz>sizeof(im3D_pars4save)) { printf("Illegal Drop format\n"); return true; }
      im3DIopt_shift += rs;
    }
    rs=read(optfid, this, sz); printf(" & %ldB im3D of %ldB\n", rs, sizeof(im3D_pars4save)); im3DIopt_shift += rs;
    initTex();
  } return true;
  case 'f': if(optfid>=0) lseek(optfid,0, SEEK_SET); recalc_at_once=false; return true;
  case 'F': if(optfid>=0) lseek(optfid,0, SEEK_END); recalc_at_once=false; return true;
  case 23: baseTitleFlag ++; return true; //recalc_at_once=false;
  case 22: recalc_at_once=false;
    printf("%s\nFrame %d (%.2f/%.2f fps), last run Times: %7.2f msec\n", WinTitle, imHost.nFrame, SmoothFPS, 1000./runTime, runTime);
    return true;
  case 'v': recalc_at_once=false; type_diag_flag++; return true;
  case 'V': recalc_at_once=false; type_diag_flag--; return true;
  case 'S': recalc_at_once=save_gp(); return true;
  case 's': recalc_at_once=save_png(imHost.nFpng); return true;
  case 'e': randR *= sqrt(sqrt(2)); return true;
  case 'E': randR /= sqrt(sqrt(2)); return true;
  case 'm': tstep /= sqrt(sqrt(2)); density /= sqrt(sqrt(2)); return true;
  case 'M': tstep *= sqrt(sqrt(2)); density *= sqrt(sqrt(2)); return true;
  case 'd': density *= sqrt(sqrt(2)); return true;
  case 'D': density /= sqrt(sqrt(2)); return true;
  case 'j': opacity = 1.0 - (1.0-opacity)/sqrt(sqrt(2)); return true;
  case 'J': opacity = 1.0 - (1.0-opacity)*sqrt(sqrt(2)); return true;
  case '@': draw_bmp4backgrownd = (draw_bmp4backgrownd+1)%4; return true;
  case '#': draw_mesh_flag ^= true; return true;
  case '$': draw_box_flag ^= true; return true;
  case '%': draw_fg_flag ^= true; return true;
  case  4 : dissect_box_flag ^= true; return true;
  case '~': draw_sec_xyz_flag ^= true; return true;
  case  6 : filterMode_flag ^= true; initTex(); return true;
  case '!': save_bmp4backgrownd(); return true;
  case '2': render_type=2; initTex(); return true;
  case '3': render_type=3; initTex(); return true;
  case '4': mode3D=(mode3D+1)%3; return true;//grad_mode ^= true; imHost.palDim = 1 + 2*grad_mode; return true;
  case '5': imHost.pal3Daxis = (imHost.pal3Daxis+1)%3; return true;
  case 'g': recalc_always=false; return true;
  case 'G': recalc_always=true; return true;
  case 'Q': recalc_at_once=false; icalcNdrop.add_sensor(ix0, iy0, iz0); return true;
  case 'q': icalcNdrop.step(); return true;
  case 11 : contour_flag ^= true; return true;
  case 'k': contour_width *= 1.2; return true;
  case 'K': contour_width /= 1.2; return true;
  case 12 : {
    FILE* cmd=fopen("add.opt", "r"); if(cmd) {
      char str[1024],* argv[2]; argv[0] = str;
      while(fgets(str, 1024, cmd)) {
        char* c=strchr(str, ' ');
        if(c) {
          if(*c==' ') *c = 0;
          argv[1] = c+1;
        } else argv[1] = str;
        init_from_command_line(argv);
      }
      fclose(cmd);
    }}
    //recalc_at_once=false;
    return true;
  case 17: recalc_at_once=false; icalcNdrop.plot_sensors(); return true;
  case 'O': recalc_at_once=false; save_section(); return true;
  case '1': sec1Daxis = (sec1Daxis+1)%4;
    printf("1D section for gnuplot set to %c\n","xyzA"[sec1Daxis]);
  case 15: recalc_at_once=false;
    if(gpPipe) { pclose(gpPipe); gpPipe = 0; }
  case 'o': recalc_at_once=false; plot_section(); return true;
  case 'b': xyz = &icalcNdrop; return true;
  case 26:
    RotPoint[0] = float(ix0)/Nx;
    RotPoint[1] = float(iy0)/Ny;
    RotPoint[2] = float(iz0)/Nz;
    return true;
  case 'x': case 'X': case 'y': case 'Y': case 'z': case 'Z': case 'u': case 'U':
    if(render_type==2 || draw_sec_xyz_flag) { switch(key) {
        case 'x': xyz2D.set(&ix0, Nx, 1); break;
        case 'X': xyz2D.set(&ix0, Nx,-1); break;
        case 'y': xyz2D.set(&iy0, Ny, 1); break;
        case 'Y': xyz2D.set(&iy0, Ny,-1); break;
        case 'z': xyz2D.set(&iz0, Nz, 1); break;
        case 'Z': xyz2D.set(&iz0, Nz,-1); break;
        default: return true;
      } xyz = &xyz2D; xyz->step();
    } else if(render_type==3) { switch(key) {
        case 'x': xyz3D.set(&viewRotation[0], 0.5f); break;
        case 'X': xyz3D.set(&viewRotation[0],-0.5f); break;
        case 'y': xyz3D.set(&viewRotation[1], 0.5f); break;
        case 'Y': xyz3D.set(&viewRotation[1],-0.5f); break;
        case 'z': xyz3D.set(&viewRotation[2], 0.5f); break;
        case 'Z': xyz3D.set(&viewRotation[2],-0.5f); break;
        case 'u': xyz3D.set(&viewTranslation[2], 0.01f); break;
        case 'U': xyz3D.set(&viewTranslation[2],-0.01f); break;
      } xyz = &xyz3D; xyz->step();
    }
    return true;
  case 27: clear4exit(); exit(0);
  default:
    if(imHost.key_func(key, x, y)) return true;
  }
  recalc_at_once=false;
  if(rs==0) return false;
  return false;
}
struct MKstates {
  int ox, oy;
  int buttonState;
  int modState;
  MKstates(): ox(0),oy(0), buttonState(0),modState(0)  {}
  void correct_screen_coor(int& x, int& y) {
    x -= im2D.xPos;
    y += im2D.yPos-(glutGet(GLUT_WINDOW_HEIGHT)-im2D.Ny);
  }
  void grabState(int button, int state, int x, int y) {
    modState = glutGetModifiers();
    if(state == GLUT_DOWN) buttonState  |= 1<<button;
    else if(state == GLUT_UP) buttonState = 0;
    ox = x;
    oy = y;
  }
} mk_state;

bool im3D_pars::special_func(unsigned char key, int x, int y) {
  mk_state.correct_screen_coor(x,y);
  if(type_diag_flag>=2) printf("special_func, keyN=%d, coors=(%d,%d)\n", key, x, y);
  recalc_at_once=true;
  size_t rs=0;
  if(key != addTitleStr[TitleStrInd]) addTitleStr[(TitleStrInd++)%4] = key;
  int modState = glutGetModifiers(), zoom=10;
  if(modState == GLUT_ACTIVE_CTRL) zoom *= 100;
  if(modState == GLUT_ACTIVE_SHIFT) zoom *= 10;
  if(modState == GLUT_ACTIVE_ALT) zoom /= 10;
  switch(key) {
  case GLUT_KEY_PAGE_UP: im2D.yPos = glutGet(GLUT_WINDOW_HEIGHT)-im2D.Ny; return true;
  case GLUT_KEY_PAGE_DOWN: im2D.yPos = 0; return true;
  case GLUT_KEY_DOWN: im2D.yPos += zoom; if(im2D.yPos>0) im2D.yPos=0; return true;
  case GLUT_KEY_UP: im2D.yPos -= zoom; {
    int yPosMax=glutGet(GLUT_WINDOW_HEIGHT)-im2D.Ny;
    if(im2D.yPos<yPosMax) im2D.yPos = yPosMax;
  } return true;
  case GLUT_KEY_HOME: im2D.xPos = 0; return true;
  case GLUT_KEY_END: im2D.xPos = glutGet(GLUT_WINDOW_WIDTH)-im2D.Nx; return true;
  case GLUT_KEY_LEFT: im2D.xPos += zoom; if(im2D.xPos>0) im2D.xPos=0; return true;
  case GLUT_KEY_RIGHT: im2D.xPos -= zoom; {
    int xPosMax=glutGet(GLUT_WINDOW_WIDTH)-im2D.Nx;
    if(im2D.xPos<xPosMax) im2D.xPos = xPosMax;
  } return true;
  }
  recalc_at_once=false;
  if(rs==0) return false;
  return false;
}
void changeCameraRange(float x, float y){};
void draw_scale(){};
void im3D_pars::mouse_func(int button, int state, int x, int y) {
  if(y<20 && state == GLUT_DOWN && !imHost.draw_flag) {
    changeCameraRange(float(x)/float(bNx), -1);
    return;
  }
  if(y<20 && state == GLUT_UP && !imHost.draw_flag) {
    changeCameraRange(-1, float(x)/float(bNx));
    return;
  }
  mk_state.correct_screen_coor(x,y);
  if(y<20 && state == GLUT_DOWN) {
    if(imHost.transparency_discrete_flag) {
      int ic=floor(0.5+(imHost.pscale)*float(x)/float(bNx));
      switch(button) {
        case 0: imHost.transparency_mode |= (1<<ic); break;
        case 1: imHost.transparency_mode ^= (1<<ic); break;
        case 2: imHost.transparency_mode &= ~(1<<ic); break;
      };
    } else {
    float f=imHost.fmin + x/float(bNx)*(imHost.fmax-imHost.fmin);
    switch(button) {
      case 0: imHost.set_lim(f,imHost.fmax); break; 
      case 2: imHost.set_lim(imHost.fmin,f); break; 
      case 1:
      float df=(f-imHost.fmin)>(imHost.fmax-f)?(f-imHost.fmin):(imHost.fmax-f);
      imHost.set_lim(f-df,f+df); break; 
    };
    if(type_diag_flag>=3) printf("mouse pal: %d,%d, button %d, state %d\n", x,y, button, state);
    recalc_at_once=true;
    }
    return;
  }
  mk_state.grabState(button, state, x,y);
  if(render_type==3) {
    if (state == GLUT_UP) {
      RotMatr R=RotMatr(0,viewRotation[0]), Ry=RotMatr(1,viewRotation[1]), Rz=RotMatr(2,viewRotation[2]), RxT=RotMatr(0,viewRotationTmp[0]), RyT=RotMatr(1,viewRotationTmp[1]);
      R *= Ry; R *= Rz; R *= RxT; R *= RyT;
      /*for(int i=0; i<3; i++) {
        printf("(");
        float s2=0.0;
        for(int j=0; j<3; j++) { printf("\t%g", R.v[i][j]); s2 += R.v[i][j]*R.v[i][j]; }
        printf("); %g\n", s2);
      }*/
      //printf("Mouse: (%g,%g,%g)+(%g,%g) -> ", viewRotation[0], viewRotation[1], viewRotation[2], viewRotationTmp[0], viewRotationTmp[1]);
      double Sy=-R.v[2][0], Cy=sqrt(1.-Sy*Sy), phi[3];
      phi[1] = atan2(Sy,Cy);
      if(Cy>0) {
        double Sx=R.v[2][1]/Cy, Cx=R.v[2][2]/Cy; phi[0]=atan2(Sx,Cx);
        double Sz=R.v[1][0]/Cy, Cz=R.v[0][0]/Cy; phi[2]=atan2(Sz,Cz);
      } else {
        double Cxz=R.v[1][1], Sxz=R.v[0][1]*Sy;
        phi[0]=atan2(Sxz, Cxz); phi[2]=0;
      }
      for(int i=0; i<3; i++) viewRotationTmp[i] = 0;
      for(int i=0; i<3; i++) viewRotation[i] = phi[i]*180.0/M_PI;
      //printf(" (%g,%g,%g)\n", viewRotation[0], viewRotation[1], viewRotation[2]);
    }
  } else {
    if (state == GLUT_DOWN  && mk_state.modState != GLUT_ACTIVE_CTRL) { if(0<=x && x<bNx && 0<=y && y<bNy) reset0(x,bNy-1-y); }
  }
  recalc_at_once=true;
  glutPostRedisplay();
}

void im3D_pars::motion_func(int x, int y) {
  mk_state.correct_screen_coor(x,y);
  if(type_diag_flag>=3) printf("motion func: %d,%d -> %d,%d\n",mk_state.ox,mk_state.oy, x,y);
  if(y<20) {
    return;
  }
  float dx, dy;
  dx = (float)(x - mk_state.ox);
  dy = (float)(y - mk_state.oy);

  if(render_type==2) {
  if(mk_state.modState == GLUT_ACTIVE_CTRL) {
    shift0(mk_state.ox,bNy-1-mk_state.oy, x,bNy-1-y);
  }
  } else {
  if(mk_state.modState == GLUT_ACTIVE_CTRL) {
    eyePoint.x = x;
    eyePoint.y = bNy-y;
  } else {
    if (mk_state.buttonState == 4) // right = zoom
      viewTranslation[2] += dy / 100.0f;
    else if (mk_state.buttonState == 2) { // middle = translate
      viewTranslation[0] += dx / 100.0f;
      viewTranslation[1] -= dy / 100.0f;
    }
    else if (mk_state.buttonState == 1) { // left = rotate
      viewRotationTmp[0] += dy / 5.0f; viewRotationTmp[1] += dx / 5.0f;
    }
  }
  }

  mk_state.ox = x;
  mk_state.oy = y;
  recalc_at_once=true;
  glutPostRedisplay();
}
//int cfX=0, cfY=0;

__global__ void im3Dclear(uchar4 bgk_col) {
  int x=blockIdx.x*blockDim.x+threadIdx.x;
  int y=blockIdx.y*blockDim.y+threadIdx.y;
  if(y<im3D.bNy && x<im3D.bNx) im.bmp[x+y*im3D.bNx] = bgk_col;
}
template<int cx, int cy, int cz>
__global__ void im3Ddraw_any(int sh, int i0) {
  int x1=blockIdx.x*blockDim.x+threadIdx.x, x2=blockIdx.y*blockDim.y+threadIdx.y;
  int p1=sh%im3D.bNx+x1, p2=sh/im3D.bNx+x2;
  if(0>p1 || p1>=im3D.bNx || 0>p2 || p2>=im3D.bNy) return;
  int ix = cx==0?i0:((cx==1?x1:x2)*im3D.x_zoom);
  int iy = cy==0?i0:((cy==1?x1:x2)*im3D.y_zoom);
  int iz = cz==0?i0:((cz==1?x1:x2)*im3D.z_zoom);
  if(ix<im3D.Nx && iy<im3D.Ny && iz<im3D.Nz) {
    uchar4 c=im.get_color(tex3D(data3D_tex, ix,iy,iz));
    if(im3D.draw_sec_xyz_flag && (abs(ix-im3D.ix0)<20 && abs(iy-im3D.iy0)<20 && abs(iz-im3D.iz0)<20) && (cx>0 && ix==im3D.ix0 || cy>0 && iy==im3D.iy0 || cz>0 && iz==im3D.iz0)) c = make_uchar4(255-c.x,255-c.y,255-c.z,c.w);
#if DATA_VECTOR_SZ==1
    if(im3D.contour_flag) {
      for(int i=0; i<im3D.cntr_num; i++) {
        float vp=tex3D(data3D_tex, ix+im3D.contour_width,iy,iz);
        float vm=tex3D(data3D_tex, ix-im3D.contour_width,iy,iz);
        float lv=im3D.cntr_levels[i];
        if(vp != 0 && vm != 0 && (vp>0 ^ vm<0) && (vp>lv ^ vm>lv)) { c = make_uchar4(255-c.x,255-c.y,255-c.z,c.w); continue; }
        vp=tex3D(data3D_tex, ix,iy+im3D.contour_width,iz);
        vm=tex3D(data3D_tex, ix,iy-im3D.contour_width,iz);
        if(vp != 0 && vm != 0 && (vp>0 ^ vm<0) && (vp>lv ^ vm>lv)) { c = make_uchar4(255-c.x,255-c.y,255-c.z,c.w); continue; }
        vp=tex3D(data3D_tex, ix,iy,iz+im3D.contour_width);
        vm=tex3D(data3D_tex, ix,iy,iz-im3D.contour_width);
        if(vp != 0 && vm != 0 && (vp>0 ^ vm<0) && (vp>lv ^ vm>lv)) { c = make_uchar4(255-c.x,255-c.y,255-c.z,c.w); continue; }
      }
      //if((1<v && v<1623/1536.) || (-1>v && v>-1623/1536.)) c = make_uchar4(255-c.x,255-c.y,255-c.z,c.w);
    }
#endif
    im.bmp[sh+x1+x2*im3D.bNx] = c;
  }
  //if(x1==128 && x2==128) printf("res(%d,%d,%d)=%g\n", ix,iy,iz, tex3D(data3D_tex, ix,iy,iz));
}
__global__ void draw_pal() {
  int x=blockIdx.x*blockDim.x+threadIdx.x;
  uchar4 col=im.get_color(im.fmin+(float(x)/im3D.bNx)*(im.fmax-im.fmin));
  uchar4* bmp = im.bmp+im3D.pal_sh;
  for(int y=0; y<20; y++, bmp += im3D.bNx) bmp[x] = col;
}
__global__ void draw_wavelength_pal(){
  int x=blockIdx.x*blockDim.x+threadIdx.x;
  uchar4 col=make_uchar4(0,0,0,255);
  float xpos = float(x)/im3D.bNx;
  uchar4* bmp = im.bmp+im3D.pal_sh; 
  for(int y=0; y<20; y++, bmp += im3D.bNx) {
    if(x%10==0 && y>10) col=make_uchar4(255,255,255,255);
    else                col=make_uchar4(0,0,0,255);
    bmp[x] = col;
  }
  float wRange = im3D.wlRange.y-im3D.wlRange.x;
  float2 cameraRange = im3D.cameraRange;
  if(xpos*wRange>cameraRange.x-im3D.wlRange.x && xpos*wRange<cameraRange.y-im3D.wlRange.x) {
    col=make_uchar4(155,0,0,255);
    bmp = im.bmp+im3D.pal_sh;
    for(int y=0; y<20; y++, bmp += im3D.bNx) bmp[x] = col;
  }
}
__global__ void negate() {
  int x=blockIdx.x*blockDim.x+threadIdx.x;
  uchar4 col=make_uchar4(255,255,255,255);
  uchar4* bmp = im.bmp+x;
  for(int y=0; y<im3D.bNy; y++) bmp[y*im3D.bNx] = col-bmp[y*im3D.bNx];
}
float invViewMatrix[12];
typedef struct {
  float4 m[3];
} float3x4;

//Код 3D рендеринга позаимствован из примеров cuda5.5: 2_Graphics/volumeRender/volumeRender_kernel.cu
__constant__ float3x4 c_invViewMatrix;  // inverse view matrix
struct Ray {
  float3 o;   // origin
  float3 d;   // direction
};

__device__ int intersectBox(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar) {
  // compute intersection of ray with all six bbox planes
  float3 invR = make_float3(1.0f) / (r.d+1e-5);
  float3 tbot = invR * (boxmin - r.o);
  float3 ttop = invR * (boxmax - r.o);

  // re-order intersections to find smallest and largest on each axis
  float3 tmin = fminf(ttop, tbot);
  float3 tmax = fmaxf(ttop, tbot);

  // find the largest tmin and the smallest tmax
  float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
  float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

  *tnear = largest_tmin;
  *tfar = smallest_tmax;
  if(im3D.dissect_box_flag) {
    float3 boxmid=boxmin+make_float3(im3D.BoxFactor[0]*(im3D.ix0+1), im3D.BoxFactor[1]*(im3D.iy0+1), im3D.BoxFactor[2]*(im3D.iz0+1));
    float3 ttopC= invR * (boxmid - r.o);
    float3 tminC = fminf(ttopC, tbot);
    float3 tmaxC = fmaxf(ttopC, tbot);
    float largest_tminC = fmaxf(fmaxf(tminC.x, tminC.y), fmaxf(tminC.x, tminC.z));
    float smallest_tmaxC = fminf(fminf(tmaxC.x, tmaxC.y), fminf(tmaxC.x, tmaxC.z));
    if(smallest_tmaxC > largest_tminC && largest_tmin == largest_tminC) *tnear = smallest_tmaxC;
  }

  return smallest_tmax > largest_tmin;
}

// transform vector by matrix (no translation)
__device__
float3 mul(const float3x4 &M, const float3 &v)
{
    float3 r;
    r.x = dot(v, make_float3(M.m[0]));
    r.y = dot(v, make_float3(M.m[1]));
    r.z = dot(v, make_float3(M.m[2]));
    return r;
}

// transform vector by matrix with translation
__device__
float4 mul(const float3x4 &M, const float4 &v)
{
    float4 r;
    r.x = dot(v, M.m[0]);
    r.y = dot(v, M.m[1]);
    r.z = dot(v, M.m[2]);
    r.w = 1.0f;
    return r;
}

__device__ uchar4 rgbaFloatToInt(float4 rgba, uchar4 bk) {
  float a=rgba.w, da=(1.-a)/255.;
  rgba.x = __saturatef(bk.x*da+a*rgba.x);   // clamp to [0.0, 1.0]
  rgba.y = __saturatef(bk.y*da+a*rgba.y);
  rgba.z = __saturatef(bk.z*da+a*rgba.z);
  rgba.w = __saturatef(rgba.w);
  return make_uchar4((rgba.x*255.f), (rgba.y*255.f), (rgba.z*255.f), (rgba.w*255.f));
}

__device__ uchar4 rgbaFloatToInt(float4 rgba) {
  rgba.x = __saturatef(rgba.x);   // clamp to [0.0, 1.0]
  rgba.y = __saturatef(rgba.y);
  rgba.z = __saturatef(rgba.z);
  rgba.w = __saturatef(rgba.w);
  return make_uchar4((rgba.x*255.f), (rgba.y*255.f), (rgba.z*255.f), (rgba.w*255.f));
}
__global__ void draw_pal3D() {
  float x=2.0f*(0.5f+blockIdx.x)/gridDim.x-1.0f, y=2.0f*(0.5f+threadIdx.x)/blockDim.x-1.0f;
  float r2=x*x+y*y;
  if(r2>1.0f) return;
  float r=sqrt(r2), r1=sqrt(1.0f-r2);
  uchar4* bmp = im.bmp+(im3D.pal_sh+im3D.bNx*(int(threadIdx.x)-int(blockDim.x/2))+blockIdx.x);
  bmp[0] = rgbaFloatToInt(im.get_color_for3D(make_float4(x,y,0,1)));
  bmp[blockDim.x] = rgbaFloatToInt(im.get_color_for3D(make_float4(0,y,x,1)));
  bmp[2*blockDim.x] = rgbaFloatToInt(im.get_color_for3D(make_float4(x,0,y,1)));
}

__device__ float smooth(float x) { return __saturatef(1.0f-x*x); } 
//------------------------------
inline __device__ void set_boxMinMax(float3& boxMin, float3& boxMax) {
  float3 boxSize=make_float3(im3D.BoxFactor[0]*im3D.Nx, im3D.BoxFactor[1]*im3D.Ny, im3D.BoxFactor[2]*im3D.Nz);
  //boxMax = 0.5f*boxSize;
  //boxMin =-0.5f*boxSize;
  //boxMax = boxSize;
  float3 cntr=(float3&)im3D.RotPoint*boxSize;
  boxMax = boxSize-cntr;
  boxMin =-cntr;
}
inline __device__ void set_eyeRay(Ray& eyeRay, float x, float y) {
  const float dbNxy=2.0f/(im3D.bNx+im3D.bNy);
  const int Nsum=im3D.Nx+im3D.Ny+im3D.Nz;
  eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 0.32f*Nsum)));
  eyeRay.d = normalize(make_float3((x-im3D.eyePoint.x)*dbNxy, (y-im3D.eyePoint.y)*dbNxy, -2.0f));
  eyeRay.d = mul(c_invViewMatrix, eyeRay.d);
}

__device__ uchar4& get_backgrownd(Ray r, float3 boxmin, float3 boxmax, int bmp_sh) {
  float3 bkgr_col=(float3&)im3D.bkgr_col, box_shrink=(float3&)im3D.box_shrink;
  float3 boxMin=box_shrink*boxmin, boxMax=box_shrink*boxmax;
  float3 fcol=make_float3(0);
  uchar4& vbmp=im.bmp[bmp_sh];
  if(im3D.draw_bmp4backgrownd && im.bmp4backgrownd != 0) vbmp = im.bmp4backgrownd[bmp_sh];
  else { fcol = bkgr_col; vbmp = make_uchar4(0,0,0,0); }
  if(im3D.draw_mesh_flag || im3D.draw_box_flag) {
    float3 invR = make_float3(1.0f) / (r.d+1e-5);
    float3 tB = invR * (boxMin - r.o);
    float3 tT = invR * (boxMax - r.o);

    float tz=r.d.z<0?tB.z:tT.z, xZ=r.o.x+r.d.x*tz, yZ=r.o.y+r.d.y*tz;
    float ty=r.d.y<0?tB.y:tT.y, zY=r.o.z+r.d.z*ty, xY=r.o.x+r.d.x*ty;
    float tx=r.d.x<0?tB.x:tT.x, yX=r.o.y+r.d.y*tx, zX=r.o.z+r.d.z*tx;
    float mval=im3D.Dmesh;
    float3 mb=(float3&)im3D.MeshBox;
    float3 ms=(float3&)im3D.MeshShift;
    if(im3D.draw_box_flag) {
      float xZn=xZ-boxmin.x, yZn=yZ-boxmin.y, xZx=boxmax.x-xZ, yZx=boxmax.y-yZ;
      float zYn=zY-boxmin.z, xYn=xY-boxmin.x, zYx=boxmax.z-zY, xYx=boxmax.x-xY;
      float yXn=yX-boxmin.y, zXn=zX-boxmin.z, yXx=boxmax.y-yX, zXx=boxmax.z-zX;
      float zval=im3D.Dmesh, dm=im3D.Dmesh;
      if(xZn>=-dm && yZn>=-dm && xZx>=-dm && yZx>=-dm) {
        if(im3D.draw_mesh_flag) { mval=fminf(mval,fminf(fabsf(remainderf(xZ-ms.x, mb.x)), fabsf(remainderf(yZ-ms.y, mb.y)))); }
        zval=fminf(zval,fminf(fminf(fabs(xZn), fabs(yZn)), fminf(fabs(xZx), fabs(yZx))));
      }
      if(zYn>=-dm && xYn>=-dm && zYx>=-dm && xYx>=-dm) {
        if(im3D.draw_mesh_flag) { mval=fminf(mval,fminf(fabsf(remainderf(zY-ms.z, mb.z)), fabsf(remainderf(xY-ms.x, mb.x)))); }
        zval=fminf(zval,fminf(fminf(fabs(xYn), fabs(zYn)), fminf(fabs(xYx), fabs(zYx))));
      }
      if(yXn>=-dm && zXn>=-dm && yXx>=-dm && zXx>=-dm) {
        if(im3D.draw_mesh_flag) { mval=fminf(mval,fminf(fabsf(remainderf(yX-ms.y, mb.y)), fabsf(remainderf(zX-ms.z, mb.z)))); }
        zval=fminf(zval,fminf(fminf(fabs(zXn), fabs(yXn)), fminf(fabs(zXx), fabs(yXx))));
      }
      float zdel=smooth(zval/im3D.Dmesh);
      fcol = fcol*(1.0f-zdel)+((float3&)(im3D.box_col))*zdel;
    } else {
           if(xZ>=boxmin.x && yZ>=boxmin.y && xZ<=boxmax.x && yZ<=boxmax.y) mval=fminf(fabsf(remainderf(xZ-ms.x, mb.x)), fabsf(remainderf(yZ-ms.y, mb.y)));
      else if(zY>=boxmin.z && xY>=boxmin.x && zY<=boxmax.z && xY<=boxmax.x) mval=fminf(fabsf(remainderf(zY-ms.z, mb.z)), fabsf(remainderf(xY-ms.x, mb.x)));
      else if(yX>=boxmin.y && zX>=boxmin.z && yX<=boxmax.y && zX<=boxmax.z) mval=fminf(fabsf(remainderf(yX-ms.y, mb.y)), fabsf(remainderf(zX-ms.z, mb.z)));
    }
    if(im3D.draw_mesh_flag) {
      float mdel=smooth(mval/im3D.Dmesh);
      fcol = fcol*(1.0f-mdel)+((float3&)(im3D.mesh_col))*mdel;
    }
  }
  vbmp = vbmp+make_uchar4(__saturatef(fcol.x)*255, __saturatef(fcol.y)*255, __saturatef(fcol.z)*255, 255);
  return vbmp;
}

__device__ uchar4& get_foregrownd(Ray r, float3 boxmin, float3 boxmax, int bmp_sh) {
  float3 box_shrink=(float3&)im3D.box_shrink;
  float3 boxMin=box_shrink*boxmin, boxMax=box_shrink*boxmax;
  float3 fcol=make_float3(0);
  uchar4& vbmp=im.bmp[bmp_sh];
  if(im3D.draw_mesh_flag || im3D.draw_box_flag) {
    float3 invR = make_float3(1.0f) / (r.d+1e-5);
    float3 tB = invR * (boxMin - r.o);
    float3 tT = invR * (boxMax - r.o);

    float tz=r.d.z>0?tB.z:tT.z, xZ=r.o.x+r.d.x*tz, yZ=r.o.y+r.d.y*tz;
    float ty=r.d.y>0?tB.y:tT.y, zY=r.o.z+r.d.z*ty, xY=r.o.x+r.d.x*ty;
    float tx=r.d.x>0?tB.x:tT.x, yX=r.o.y+r.d.y*tx, zX=r.o.z+r.d.z*tx;
    float mval=im3D.Dmesh;
    float3 mb=(float3&)im3D.MeshBox;
    float3 ms=(float3&)im3D.MeshShift;
    if(im3D.draw_box_flag) {
      float xZn=xZ-boxmin.x, yZn=yZ-boxmin.y, xZx=boxmax.x-xZ, yZx=boxmax.y-yZ;
      float zYn=zY-boxmin.z, xYn=xY-boxmin.x, zYx=boxmax.z-zY, xYx=boxmax.x-xY;
      float yXn=yX-boxmin.y, zXn=zX-boxmin.z, yXx=boxmax.y-yX, zXx=boxmax.z-zX;
      float zval=im3D.Dmesh, dm=im3D.Dmesh;
      if(xZn>=-dm && yZn>=-dm && xZx>=-dm && yZx>=-dm) {
        if(im3D.draw_mesh_flag) { mval=fminf(mval,fminf(fabsf(remainderf(xZ-ms.x, mb.x)), fabsf(remainderf(yZ-ms.y, mb.y)))); }
        zval=fminf(zval,fminf(fminf(fabs(xZn), fabs(yZn)), fminf(fabs(xZx), fabs(yZx))));
      }
      if(zYn>=-dm && xYn>=-dm && zYx>=-dm && xYx>=-dm) {
        if(im3D.draw_mesh_flag) { mval=fminf(mval,fminf(fabsf(remainderf(zY-ms.z, mb.z)), fabsf(remainderf(xY-ms.x, mb.x)))); }
        zval=fminf(zval,fminf(fminf(fabs(xYn), fabs(zYn)), fminf(fabs(xYx), fabs(zYx))));
      }
      if(yXn>=-dm && zXn>=-dm && yXx>=-dm && zXx>=-dm) {
        if(im3D.draw_mesh_flag) { mval=fminf(mval,fminf(fabsf(remainderf(yX-ms.y, mb.y)), fabsf(remainderf(zX-ms.z, mb.z)))); }
        zval=fminf(zval,fminf(fminf(fabs(zXn), fabs(yXn)), fminf(fabs(zXx), fabs(yXx))));
      }
      float zdel=smooth(zval/im3D.Dmesh);
      fcol = fcol*(1.0f-zdel)+((float3&)(im3D.box_col))*zdel;
    } else {
           if(xZ>=boxmin.x && yZ>=boxmin.y && xZ<=boxmax.x && yZ<=boxmax.y) mval=fminf(fabsf(remainderf(xZ-ms.x, mb.x)), fabsf(remainderf(yZ-ms.y, mb.y)));
      else if(zY>=boxmin.z && xY>=boxmin.x && zY<=boxmax.z && xY<=boxmax.x) mval=fminf(fabsf(remainderf(zY-ms.z, mb.z)), fabsf(remainderf(xY-ms.x, mb.x)));
      else if(yX>=boxmin.y && zX>=boxmin.z && yX<=boxmax.y && zX<=boxmax.z) mval=fminf(fabsf(remainderf(yX-ms.y, mb.y)), fabsf(remainderf(zX-ms.z, mb.z)));
    }
    if(im3D.draw_mesh_flag) {
      float mdel=smooth(mval/im3D.Dmesh);
      fcol = fcol*(1.0f-mdel)+((float3&)(im3D.mesh_col))*mdel;
    }
  }
  vbmp = vbmp+make_uchar4(__saturatef(fcol.x)*255, __saturatef(fcol.y)*255, __saturatef(fcol.z)*255, 255);
  return vbmp;
}

__device__ void mk_pts(int x, int y, uchar4 col) {
  const int ps=2;
  if(x+1<ps || x+ps>=im3D.bNx || y+1<ps || y+ps>=im3D.bNy) return;
  for(int ix=1-ps; ix<ps; ix++) for(int iy=1-ps; iy<ps; iy++)
    im.bmp[(iy+y)*im3D.bNx + x+ix] = col;
}
__device__ void mk_box(int x, int y, uchar4 col) {
  if(x<0 || x+tileSz>=im3D.bNx || y<0 || y+tileSz>=im3D.bNy) return;
  for(int ix=0; ix<tileSz; ix++) im.bmp[y*im3D.bNx + x+ix] = im.bmp[(tileSz+y)*im3D.bNx + x+ix] = col;
  for(int iy=0; iy<tileSz; iy++) im.bmp[(iy+y)*im3D.bNx + x] = im.bmp[(iy+y)*im3D.bNx + x+tileSz] = col;
}
inline bool __device__ is_inside(float2 pt, float2 p0, float2 px, float2 py) {
  float v1=(p0.x - pt.x) * (px.y - p0.y) - (px.x - p0.x) * (p0.y - pt.y);
  float v2=(px.x - pt.x) * (py.y - px.y) - (py.x - px.x) * (px.y - pt.y);
  float v3=(py.x - pt.x) * (p0.y - py.y) - (p0.x - py.x) * (py.y - pt.y);
  return (v1*v2>=0.0 && v1*v3>=0.0 && v2*v3>=0.0);
}
inline float2 __device__ pt_inside(float2 pt, float2 p0, float2 px, float2 py) {
  float2 res;
  res.x = ((pt.x-p0.x)*(py.y-p0.y)-(pt.y-p0.y)*(py.x-p0.x))/((px.x-p0.x)*(py.y-p0.y)-(px.y-p0.y)*(py.x-p0.x));
  res.y = ((pt.x-p0.x)*(px.y-p0.y)-(pt.y-p0.y)*(px.x-p0.x))/((py.x-p0.x)*(px.y-p0.y)-(py.y-p0.y)*(px.x-p0.x));
  return res;
}
__global__ void save_gp3D() {
  __shared__ float2 fm[3][tilesN][tilesN];//координаты точки в области с сеткой
  __shared__ int hit[tilesN][tilesN];//индекс области попадания луча: 1-z 2-y 4-x 0-молоко
  const int Sgp=(tilesN-1)*tileSz;
  int x=blockIdx.x*Sgp+threadIdx.x*tileSz, y=blockIdx.y*Sgp+threadIdx.y*tileSz;
  float3 boxMin, boxMax; set_boxMinMax(boxMin, boxMax);
  boxMax=((float3&)im3D.box_shrink)*boxMax;
  boxMin=((float3&)im3D.box_shrink)*boxMin;
  Ray r; set_eyeRay(r, x,y);
  float3 invR = make_float3(1.0f) / (r.d+1e-5);
  float3 tB = invR * (boxMin - r.o);
  float3 tT = invR * (boxMax - r.o);
  float tz=r.d.z<0?tB.z:tT.z, xZ=r.o.x+r.d.x*tz, yZ=r.o.y+r.d.y*tz;
  float ty=r.d.y<0?tB.y:tT.y, zY=r.o.z+r.d.z*ty, xY=r.o.x+r.d.x*ty;
  float tx=r.d.x<0?tB.x:tT.x, yX=r.o.y+r.d.y*tx, zX=r.o.z+r.d.z*tx;
  fm[2][threadIdx.x][threadIdx.y] = make_float2(xZ, yZ);
  fm[1][threadIdx.x][threadIdx.y] = make_float2(zY, xY);
  fm[0][threadIdx.x][threadIdx.y] = make_float2(yX, zX);
  if(xZ>=boxMin.x && yZ>=boxMin.y && xZ<=boxMax.x && yZ<=boxMax.y) hit[threadIdx.x][threadIdx.y] = 1; //mk_pts(x,y, red);}
  else if(zY>=boxMin.z && xY>=boxMin.x && zY<=boxMax.z && xY<=boxMax.x) hit[threadIdx.x][threadIdx.y] = 2; //mk_pts(x,y, green);}
  else if(yX>=boxMin.y && zX>=boxMin.z && yX<=boxMax.y && zX<=boxMax.z) hit[threadIdx.x][threadIdx.y] = 4; //mk_pts(x,y, blue);}
  else hit[threadIdx.x][threadIdx.y] = 0;
  __syncthreads();

  int hitA=0, hitM=0;
  if(threadIdx.x<tilesN-1 && threadIdx.y<tilesN-1) {
    for(int i=0;i<2;i++) for(int j=0;j<2;j++) {
      int h=hit[threadIdx.x+i][threadIdx.y+j];
      if(h>0) { hitA++; hitM |= h; }
    }
  }
  int cs=abs(2*hitM-7)/2;
  if(hitA==0 || hitA==4 || cs>=3) return;
  bool is4tick=false, is4bnd=false, is4axis=false;
  is4bnd = hitM==1 || hitM==2 || hitM==4;
  is4axis= hitM==3 || hitM==5 || hitM==6;
  int cp=(cs+1)%3, cm=(cs+2)%3;
  float2 tick_sh={0.0,0.0}, tick2sh={0.0,0.0}; float tick_val;
  const float axis_gap=60., tick_gap=20.;
  float2 pt, spt={0.,0.}; float bMax[]={boxMax.x,boxMax.y,boxMax.z}, bMin[]={boxMin.x,boxMin.y,boxMin.z};
  int labN=(blockIdx.x*(tilesN-1)+threadIdx.x)+gridDim.x*(tilesN-1)*(blockIdx.y*(tilesN-1)+threadIdx.y);
  if(is4axis) {
    float2 p0=fm[cm][threadIdx.x][threadIdx.y], px=fm[cm][threadIdx.x+1][threadIdx.y], py=fm[cm][threadIdx.x][threadIdx.y+1];
    if(fabs(p0.x-bMax[cs])<fabs(p0.x-bMin[cs])) { pt.x = bMax[cs]; spt.x = axis_gap; }
    else { pt.x = bMin[cs]; spt.x = -axis_gap; }
    pt.y = fabs(p0.y-bMax[cp])<fabs(p0.y-bMin[cp])?bMax[cp]:bMin[cp];
    tick_sh = pt_inside(pt, p0,px,py);
    tick2sh = pt_inside(pt+spt, p0,px,py);
    printf("set arrow %d from %g,%g to %g,%g front nohead\n", labN, x+tick_sh.x*tileSz,y+tick_sh.y*tileSz, x+tick2sh.x*tileSz,y+tick2sh.y*tileSz);
    printf("set label %d \"%c\" at %g,%g front center\n", labN, "xyz?"[cs], x+tick2sh.x*tileSz,y+tick2sh.y*tileSz+tick_gap*((tick2sh.y<tick_sh.y)?-1.:1.));
  } else if(is4bnd) {
    float2 fmin,fmax; fmin = fmax = fm[cs][threadIdx.x][threadIdx.y];
    for(int i=0;i<2;i++) for(int j=0;j<2;j++) {
      float2 f = fm[cs][threadIdx.x+i][threadIdx.y+j];
      if(f.x<fmin.x) fmin.x = f.x;
      if(f.y<fmin.y) fmin.y = f.y;
      if(f.x>fmax.x) fmax.x = f.x;
      if(f.y>fmax.y) fmax.y = f.y;
    }
    if(fmin.x<bMin[cp] || fmax.x>bMax[cp]) {// cM = cm;
      int mmin=floorf(fmin.y/im3D.MeshBox[cm]), mmax=floorf(fmax.y/im3D.MeshBox[cm]);
      if(mmin != mmax) is4tick = true;
      pt.x = fmin.x<bMin[cp]?bMin[cp]:bMax[cp]; spt.x = fmin.x<bMin[cp]?-tick_gap:tick_gap;
      pt.y = mmax*im3D.MeshBox[cm];
      tick_val = im3D.base[cm] + pt.y*im3D.step[cm];
    } else if(fmin.y<bMin[cm] || fmax.y>bMax[cm]) {// cM = cp;
      int mmin=floorf(fmin.x/im3D.MeshBox[cp]), mmax=floorf(fmax.x/im3D.MeshBox[cp]);
      if(mmin != mmax) is4tick = true;
      pt.x = mmax*im3D.MeshBox[cp];
      pt.y = fmin.y<bMin[cm]?bMin[cm]:bMax[cm]; spt.y = fmin.y<bMin[cm]?-tick_gap:tick_gap;
      tick_val = im3D.base[cp] + pt.x*im3D.step[cp];
    }
    if(is4tick) {
      float2 p0=fm[cs][threadIdx.x][threadIdx.y], px=fm[cs][threadIdx.x+1][threadIdx.y], py=fm[cs][threadIdx.x][threadIdx.y+1], p1=fm[cs][threadIdx.x+1][threadIdx.y+1];
      if(is_inside(pt, p0,px,py)) {
        tick_sh = pt_inside(pt, p0,px,py);
        tick2sh = pt_inside(pt+spt, p0,px,py);
      } else if(is_inside(pt, p1,py,px)) {
        tick_sh = 1.0-pt_inside(pt, p1,py,px);
        tick2sh = 1.0-pt_inside(pt+spt, p1,py,px);
      } else is4tick = false;
      if(is4tick) printf("set label %d \"%g\" at %g,%g front %s\n", labN, tick_val, x+tick2sh.x*tileSz,y+tick2sh.y*tileSz, (tick2sh.x<tick_sh.x)?"right":"left");
    }
  }
  uchar4 red=make_uchar4(255,0,0,0), green=make_uchar4(0,255,0,0), blue=make_uchar4(0,0,255,0);
  uchar4 ltred=make_uchar4(128,0,0,0), ltgreen=make_uchar4(0,128,0,0), ltblue=make_uchar4(0,0,128,0);
  if(is4axis) {
    mk_box(x,y, red);
    mk_pts(x+tick2sh.x*tileSz,y+tick2sh.y*tileSz, red);
  } else if(is4tick) {
    mk_box(x,y, blue);
    mk_pts(x+tick2sh.x*tileSz,y+tick2sh.y*tileSz, blue);
  } else if(is4bnd) mk_box(x,y, green);
  else mk_box(x,y, ltblue);
}
__global__ void __launch_bounds__(1024,1) grad_render3D() {
#if DATA_VECTOR_SZ==1
  const float opacityThreshold = im3D.opacity;
  const float density=im3D.density, brightness=im.max_rgb;
  float3 boxMin, boxMax; set_boxMinMax(boxMin, boxMax);

  int x = blockIdx.x*blockDim.x + threadIdx.x;
  int y = blockIdx.y*blockDim.y + threadIdx.y;
  Ray eyeRay; set_eyeRay(eyeRay, x,y);

  uchar4& vbmp=get_backgrownd(eyeRay, boxMin, boxMax, y*im3D.bNx + x);
  float phi=im3D.randArr[threadIdx.x+threadIdx.y*blockDim.x];
  set_eyeRay(eyeRay, x+im3D.randR*cos(phi),y+im3D.randR*sin(phi));
  float tnear, tfar;
  int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

  if (!hit) return;

  if(tnear < 0.0f) tnear = 0.0f;     // clamp to near plane
  float4 sum = make_float4(0.0f);
  const float3 SzfdBox=make_float3(im3D.Nx,im3D.Ny,im3D.Nz)/(boxMax-boxMin);
  float3 pos_sc = (eyeRay.o + eyeRay.d*tnear-boxMin)*SzfdBox-0.5f;
  const float3 step_sc = (eyeRay.d*im3D.tstep)*SzfdBox;
  for(float t=tnear; t<tfar; t+=im3D.tstep, pos_sc += step_sc) {
    // cross stencil:
    float d=im3D.tstep, dd=im.max_rgb_step*0.5/d;
    float dfdx=dd*(tex3D(data3D_tex, pos_sc.x+d, pos_sc.y, pos_sc.z)-tex3D(data3D_tex, pos_sc.x-d, pos_sc.y, pos_sc.z));
    float dfdy=dd*(tex3D(data3D_tex, pos_sc.x, pos_sc.y+d, pos_sc.z)-tex3D(data3D_tex, pos_sc.x, pos_sc.y-d, pos_sc.z));
    float dfdz=dd*(tex3D(data3D_tex, pos_sc.x, pos_sc.y, pos_sc.z+d)-tex3D(data3D_tex, pos_sc.x, pos_sc.y, pos_sc.z-d));
    float4 col = im.get_color_for3D(make_float4(dfdx,dfdy,dfdz,tex3D(data3D_tex, pos_sc.x, pos_sc.y, pos_sc.z)));
    float w=col.w*density*(1.0f - sum.w); col.w = 1;
    sum += col * w;
    if(sum.w >= opacityThreshold) {
      sum -= col*(sum.w - opacityThreshold);
      break;
    }/*
    col.w *= density;
    col.x *= col.w;
    col.y *= col.w;
    col.z *= col.w;
    sum = sum + col*(1.0f - sum.w);

    if (sum.w > opacityThreshold) break;*/
  }
  sum.x *= brightness; sum.y *= brightness; sum.z *= brightness;
  vbmp = rgbaFloatToInt(sum, vbmp);
  if(im3D.draw_fg_flag) vbmp=get_foregrownd(eyeRay, boxMin, boxMax, y*im3D.bNx + x);
  if(im3D.draw_sec_xyz_flag) {
    if(fabs(pos_sc.x-im3D.ix0)<=0.5|| fabs(pos_sc.y-im3D.iy0)<=0.5|| fabs(pos_sc.z-im3D.iz0)<=0.5) vbmp = make_uchar4(255-vbmp.x,255-vbmp.y,255-vbmp.z,vbmp.w);
  }
#endif
}
__global__ void __launch_bounds__(1024,1) surf_render3D() {
#if DATA_VECTOR_SZ==1
  const float opacityThreshold = im3D.opacity;
  const float density=im3D.density, brightness=im.max_rgb;
  float3 boxMin, boxMax; set_boxMinMax(boxMin, boxMax);

  int x = blockIdx.x*blockDim.x + threadIdx.x;
  int y = blockIdx.y*blockDim.y + threadIdx.y;
  Ray eyeRay; set_eyeRay(eyeRay, x,y);

  uchar4& vbmp=get_backgrownd(eyeRay, boxMin, boxMax, y*im3D.bNx + x);
  float phi=im3D.randArr[threadIdx.x+threadIdx.y*blockDim.x];
  set_eyeRay(eyeRay, x+im3D.randR*cos(phi),y+im3D.randR*sin(phi));
  float tnear, tfar;
  int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

  if (!hit) return;

  if(tnear < 0.0f) tnear = 0.0f;     // clamp to near plane
  float4 sum = make_float4(0.0f);
  const float3 SzfdBox=make_float3(im3D.Nx,im3D.Ny,im3D.Nz)/(boxMax-boxMin);
  float3 pos_sc = (eyeRay.o + eyeRay.d*tnear-boxMin)*SzfdBox-0.5f;
  const float3 step_sc = (eyeRay.d*im3D.tstep)*SzfdBox;
  for(float t=tnear; t<tfar; t+=im3D.tstep, pos_sc += step_sc) {
    // cross stencil:
    short2 s2=tex3D(data3Dsurf_tex, pos_sc.x, pos_sc.y, pos_sc.z);
    const short MAX_SHORT=(1<<15)-1; const float dMS=1.0f/MAX_SHORT;
    float3 f={0,0,0};
    if(s2.x!=-MAX_SHORT-1 && s2.y!=-MAX_SHORT-1) {
      f.z = s2.x*dMS; float fxy=sqrt(1-f.z*f.z), phi=s2.y*dMS*M_PI;
      f.y = fxy*sin(phi);
      f.x = fxy*cos(phi);
    }
    float4 col = im.get_color_for3D(make_float4(f.x,f.y,f.z,tex3D(data3D_tex, pos_sc.x, pos_sc.y, pos_sc.z)));
    float w=col.w*density*(1.0f - sum.w); col.w = 1;
    sum += col * w;
    if(sum.w >= opacityThreshold) {
      sum -= col*(sum.w - opacityThreshold);
      break;
    }
  }
  sum.x *= brightness; sum.y *= brightness; sum.z *= brightness;
  vbmp = rgbaFloatToInt(sum, vbmp);
  if(im3D.draw_fg_flag) vbmp=get_foregrownd(eyeRay, boxMin, boxMax, y*im3D.bNx + x);
  if(im3D.draw_sec_xyz_flag) {
    if(fabs(pos_sc.x-im3D.ix0)<=0.5|| fabs(pos_sc.y-im3D.iy0)<=0.5|| fabs(pos_sc.z-im3D.iz0)<=0.5) vbmp = make_uchar4(255-vbmp.x,255-vbmp.y,255-vbmp.z,-vbmp.w);
  }
  //if(x==im3D.bNx/2 && y==im3D.bNy/2) printf("Surf: %f,%f,%f,%f*%f/%f*%f => %d,%d,%d\n", sum.x,sum.y,sum.z,sum.w, last_mul,opacityThreshold, brightness, vbmp.x,vbmp.y,vbmp.z);
#endif
}
__global__ void __launch_bounds__(1024,1) render3D() {
  const float opacityThreshold = im3D.opacity;//0.95f;
  const float density=im3D.density, brightness=im.max_rgb;
  float3 boxMin, boxMax; set_boxMinMax(boxMin, boxMax);

  int x = blockIdx.x*blockDim.x + threadIdx.x;
  int y = blockIdx.y*blockDim.y + threadIdx.y;
  //bool isCnt=blockIdx.x==gridDim.x/2 && blockIdx.y==gridDim.y/2 && threadIdx.x == blockDim.x/2 && threadIdx.y == blockDim.y/2;
  //if ((x >= im3D.bNx) || (y >= im3D.bNy)) return;
  //if(x==0 && y==0) printf("block: %gx%gx%g\n", boxMax.x, boxMax.y, boxMax.z);

  // calculate eye ray in world space
  Ray eyeRay; set_eyeRay(eyeRay, x,y);
  //const int Nsum=im3D.Nx+im3D.Ny+im3D.Nz;
  //const float dbNxy=2.0f/(im3D.bNx+im3D.bNy);
  //eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 0.32f*Nsum)));
  //eyeRay.d = normalize(make_float3((x-im3D.bNx/2)*dbNxy, (y-im3D.bNy/2)*dbNxy, -2.0f));
  //eyeRay.d = mul(c_invViewMatrix, eyeRay.d);

  uchar4& vbmp=get_backgrownd(eyeRay, boxMin, boxMax, y*im3D.bNx + x);
  float phi=im3D.randArr[threadIdx.x+threadIdx.y*blockDim.x];
  set_eyeRay(eyeRay, x+im3D.randR*cos(phi),y+im3D.randR*sin(phi));
  float tnear, tfar;
  int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

  if (!hit) return;

  if(tnear < 0.0f) tnear = 0.0f;     // clamp to near plane
  //if(tnear+im3D.tstep*Nsum<tfar) tfar = tnear+im3D.tstep*Nsum;
  // march along ray from front to back, accumulating color
  float4 sum = make_float4(0.0f);
  //float3 pos = eyeRay.o + eyeRay.d*tnear;
  //float3 step = eyeRay.d*im3D.tstep;
  const float3 SzfdBox=make_float3(im3D.Nx,im3D.Ny,im3D.Nz)/(boxMax-boxMin);
  float3 pos_sc = (eyeRay.o + eyeRay.d*tnear-boxMin)*SzfdBox-0.5f;
  const float3 step_sc = (eyeRay.d*im3D.tstep)*SzfdBox;
  //const float pscale=im.pscale*0.01f, fscale=100.0f*im.fscale, fmin=0.5f-im.fmin*fscale;
//if(isCnt) printf("I am ray: %f(%f)%f step: %f,%f,%f; pos: %f,%f,%f of %d,%d,%d\n", tnear,im3D.tstep,tfar, step_sc.x,step_sc.y,step_sc.z,  pos_sc.x, pos_sc.y, pos_sc.z, im3D.Nx,im3D.Ny,im3D.Nz);
  for(float t=tnear; t<tfar; t+=im3D.tstep, pos_sc += step_sc) {
    // read from 3D texture
    float4 col = im.get_color_for3D(tex3D(data3D_tex, pos_sc.x, pos_sc.y, pos_sc.z));
    float w=col.w*density*(1.0f - sum.w); col.w = 1;
    sum += col * w;
    if(sum.w >= opacityThreshold) {
      sum -= col*(sum.w - opacityThreshold);
      break;
    }/*
    //float f = tex3D(data3D_tex, pos_sc.x, pos_sc.y, pos_sc.z);
    //float4 col = tex1D(fpal_col_tex, 0.5f + pscale*tex1D(fpal_scale_tex, fmin+f*fscale));
    col.w *= density;

    // "under" operator for back-to-front blending
    //sum = lerp(sum, col, col.w);

    // pre-multiply alpha
    col.x *= col.w;
    col.y *= col.w;
    col.z *= col.w;
    // "over" operator for front-to-back blending
    sum = sum + col*(1.0f - sum.w);

    // exit early if opaque
    if (sum.w > opacityThreshold) break;
   // pos_sc += step_sc;
   */
  }
//if(isCnt) printf("I am ray: %f\n",sum.w);
  sum.x *= brightness; sum.y *= brightness; sum.z *= brightness;
  //sum *= brightness;

  // write output color
  vbmp = rgbaFloatToInt(sum, vbmp);
  if(im3D.draw_fg_flag) vbmp=get_foregrownd(eyeRay, boxMin, boxMax, y*im3D.bNx + x);
  //if(threadIdx.x==0 && threadIdx.y==0) vbmp = make_uchar4(255,255,255,255);
  if(im3D.draw_sec_xyz_flag) {
    if(fabs(pos_sc.x-im3D.ix0)<=0.5|| fabs(pos_sc.y-im3D.iy0)<=0.5|| fabs(pos_sc.z-im3D.iz0)<=0.5) vbmp = make_uchar4(255-vbmp.x,255-vbmp.y,255-vbmp.z,-vbmp.w);
  }
}

void im3D_pars::save_bmp4backgrownd() {
try {
  uchar4* devPtr; size_t size;
  if(CHECK_ERROR(cudaGraphicsMapResources(1, &im2D.resource, NULL))) throw(-1);
  if(imHost.negate_flag) negate <<<bNx/NW,NW>>>();
  if(CHECK_ERROR(cudaGraphicsResourceGetMappedPointer((void**) &devPtr, &size, im2D.resource))) throw(-1);
  if(imHost.bmp4backgrownd != 0) CHECK_ERROR(cudaFree(imHost.bmp4backgrownd));
  if(CHECK_ERROR(cudaMalloc((void**) &imHost.bmp4backgrownd, size))) throw(-1);
  if(CHECK_ERROR(cudaMemcpy(imHost.bmp4backgrownd, devPtr, size, cudaMemcpyDeviceToDevice))) throw(-1);
  im2D.unmapAfterDraw();
} catch(...) {
  printf("save_bmp4backgrownd: Возникла какая-то ошибка.\n");
}
}
void im3D_pars::recalc_sec_im3D() {
try {
  imHost.bmp = im2D.map4draw();
  imHost.bind2draw();
  if(CHECK_ERROR(cudaMemcpyToSymbol(im, &imHost, sizeof(imHost)))) throw(-1);
  if(CHECK_ERROR(cudaMemcpyToSymbol(im3D, this, sizeof(im3D_pars)))) throw(-1);
  int NxZ=Nx/x_zoom, NyZ=Ny/y_zoom, NzZ=Nz/z_zoom;
  int NxB=(NxZ+NW-1)/NW, NyB=(NyZ+NW-1)/NW, NzB=(NzZ+NW-1)/NW;
  unsigned char ub[3];
  for(int i=0; i<3; i++) { ub[i] = bkgr_col[i]<0?0:(bkgr_col[i]>1?255:255.*bkgr_col[i]); }
  im3Dclear <<<dim3(bNx/NW,bNy/NW),dim3(NW,NW)>>>(make_uchar4(ub[0], ub[1], ub[2], 255));
  int shX=0,shY=0;
  for(int ix=int(Nx*RotPoint[0])%int(MeshBox[0]); ix<Nx; ix+=MeshBox[0]) {
    if(shX+NyZ>bNx) { shX=0; shY += NzZ+2; } if(shY+NzZ>bNy) break;
    im3Ddraw_any<0,1,2> <<<dim3(NyB,NzB),dim3(NW,NW)>>>(shX+shY*bNx,ix);
    shX += NyZ+2;
  }// if(shX>0) { shX=0; shY += NzZ+2; }
  for(int iy=int(Ny*RotPoint[1])%int(MeshBox[1]); iy<Ny; iy+=MeshBox[1]) {
    if(shX+NxZ>bNx) { shX=0; shY += NzZ+2; } if(shY+NzZ>bNy) break;
    im3Ddraw_any<1,0,2> <<<dim3(NxB,NzB),dim3(NW,NW)>>>(shX+shY*bNx,iy);
    shX += NxZ+2;
  } if(shX>0) { shX=0; shY += NzZ+2; }
  for(int iz=int(Nz*RotPoint[2])%int(MeshBox[2]); iz<Nz; iz+=MeshBox[2]) {
    //printf("draw xy at iz=%d; (%d,%d) -> (%d,%d)..\n", iz, shX,shY, 0,shY +(NyZ+2), );
    if(shX+NxZ>bNx) { shX=0; shY += NyZ+2; } if(shY+NyZ>bNy) break;
    im3Ddraw_any<1,2,0> <<<dim3(NxB,NyB),dim3(NW,NW)>>>(shX+shY*bNx,iz);
    shX += NxZ+2;
  }// if(shX>0) { shX=0; shY += NyZ+2; }
  if(imHost.draw_flag) draw_pal <<<bNx/NW,NW>>>(); else draw_wavelength_pal <<<bNx/NW,NW>>>();
  if(imHost.negate_flag) negate <<<bNx/NW,NW>>>();
  imHost.nFrame++;
  imHost.unbindAfterDraw();
  im2D.unmapAfterDraw();
} catch(...) {
  printf("recalc_im3D: Возникла какая-то ошибка.\n");
}
}
void im3D_pars::shift0(int x, int y, int x1, int y1) {
  int ix,iy, dx=x1-x, dy=y1-y, sh=dx+dy*bNx;
  if(secType!=1) {
    ix=(x-secXsh%bNx)*z_zoom; iy=(y-secXsh/bNx)*y_zoom;
    if(0<=ix && ix<Nz && 0<=iy && iy<Ny) { if(secXsh%bNx+dx>=0 && secXsh/bNx+dy>=0) secXsh += sh; return; }
  } else {
    ix=(x-secXsh%bNx)*y_zoom; iy=(y-secXsh/bNx)*z_zoom;
    if(0<=ix && ix<Ny && 0<=iy && iy<Nz) { if(secXsh%bNx+dx>=0 && secXsh/bNx+dy>=0) secXsh += sh; return; }
  }
  if(secType<2) {
    ix=(x-secYsh%bNx)*x_zoom; iy=(y-secYsh/bNx)*z_zoom;
    if(0<=ix && ix<Nx && 0<=iy && iy<Nz) { if(secYsh%bNx+dx>=0 && secYsh/bNx+dy>=0) secYsh += sh; return; }
  } else {
    ix=(x-secYsh%bNx)*z_zoom; iy=(y-secYsh/bNx)*x_zoom;
    if(0<=ix && ix<Nz && 0<=iy && iy<Nx) { if(secYsh%bNx+dx>=0 && secYsh/bNx+dy>=0) secYsh += sh; return; }
  }
  ix=(x-secZsh%bNx)*x_zoom; iy=(y-secZsh/bNx)*y_zoom;
  if(0<=ix && ix<Nx && 0<=iy && iy<Ny) { if(secZsh%bNx+dx>=0 && secZsh/bNx+dy>=0) secZsh += sh; return; }
}
void im3D_pars::reset0(int x, int y) {
  int ix,iy;
  ix=(x-secZsh%bNx)*x_zoom; iy=(y-secZsh/bNx)*y_zoom;
  if(0<=ix && ix<Nx && 0<=iy && iy<Ny) { ix0 = ix; iy0 = iy; return; }
  if(secType<2) {
    ix=(x-secYsh%bNx)*x_zoom; iy=(y-secYsh/bNx)*z_zoom;
    if(0<=ix && ix<Nx && 0<=iy && iy<Nz) { ix0 = ix; iz0 = iy; return; }
  } else {
    ix=(x-secYsh%bNx)*z_zoom; iy=(y-secYsh/bNx)*x_zoom;
    if(0<=ix && ix<Nz && 0<=iy && iy<Nx) { iz0 = ix; ix0 = iy; return; }
  }
  if(secType!=1) {
    ix=(x-secXsh%bNx)*z_zoom; iy=(y-secXsh/bNx)*y_zoom;
    if(0<=ix && ix<Nz && 0<=iy && iy<Ny) { iz0 = ix; iy0 = iy; return; }
  } else {
    ix=(x-secXsh%bNx)*y_zoom; iy=(y-secXsh/bNx)*z_zoom;
    if(0<=ix && ix<Ny && 0<=iy && iy<Nz) { iy0 = ix; iz0 = iy; return; }
  }
}
void im3D_pars::recalc_im3D() {
try {
  imHost.bmp = im2D.map4draw();
  imHost.bind2draw();
  if(CHECK_ERROR(cudaMemcpyToSymbol(im, &imHost, sizeof(imHost)))) throw(-1);
  if(CHECK_ERROR(cudaMemcpyToSymbol(im3D, this, sizeof(im3D_pars)))) throw(-1);
  //if(CHECK_ERROR(cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeof(float4)*3))) throw(-1);
  //if(CHECK_ERROR(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared))) throw(-1);
  //Pal via Tex
  int NxB=(Nx/x_zoom+NW-1)/NW, NyB=(Ny/y_zoom+NW-1)/NW, NzB=(Nz/z_zoom+NW-1)/NW;
  unsigned char ub[3];
  for(int i=0; i<3; i++) { ub[i] = bkgr_col[i]<0?0:(bkgr_col[i]>1?255:255.*bkgr_col[i]); }
  im3Dclear <<<dim3(bNx/NW,bNy/NW),dim3(NW,NW)>>>(make_uchar4(ub[0], ub[1], ub[2], 255));
  im3Ddraw_any<1,2,0> <<<dim3(NxB,NyB),dim3(NW,NW)>>>(secZsh,iz0);
  if(secType<2) im3Ddraw_any<1,0,2> <<<dim3(NxB,NzB),dim3(NW,NW)>>>(secYsh,iy0);
  else          im3Ddraw_any<2,0,1> <<<dim3(NzB,NxB),dim3(NW,NW)>>>(secYsh,iy0);
  if(secType!=1) im3Ddraw_any<0,2,1> <<<dim3(NzB,NyB),dim3(NW,NW)>>>(secXsh,ix0);
  else           im3Ddraw_any<0,1,2> <<<dim3(NyB,NzB),dim3(NW,NW)>>>(secXsh,ix0);
  if(imHost.draw_flag) draw_pal <<<bNx/NW,NW>>>(); else draw_wavelength_pal <<<bNx/NW,NW>>>();
  if(imHost.negate_flag) negate <<<bNx/NW,NW>>>();
  imHost.nFrame++;
  imHost.unbindAfterDraw();
  im2D.unmapAfterDraw();
} catch(...) {
  printf("recalc_im3D: Возникла какая-то ошибка.\n");
}
}
void im3D_pars::recalc3D_im3D() {
try {
  // use OpenGL to build view matrix
  GLfloat modelView[16];
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  switch(mk_state.modState) {
    case GLUT_ACTIVE_SHIFT:
  glRotatef(-viewRotation[1], 0.0, 1.0, 0.0);
  glRotatef(-viewRotation[0], 1.0, 0.0, 0.0);
    break;
    case GLUT_ACTIVE_CTRL:
    default:
  glRotatef(-viewRotation[0], 1.0, 0.0, 0.0);
  glRotatef(-viewRotation[1], 0.0, 1.0, 0.0);
  glRotatef(-viewRotation[2], 0.0, 0.0, 1.0);
  glRotatef(-viewRotationTmp[0], 1.0, 0.0, 0.0);
  glRotatef(-viewRotationTmp[1], 0.0, 1.0, 0.0);
  }
  glTranslatef(-viewTranslation[0], -viewTranslation[1], -viewTranslation[2]);
  glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
  glPopMatrix();
  for(int i=0; i<12; i++) invViewMatrix[i] = modelView[4*(i&3)+i/4];
  if(CHECK_ERROR(cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeof(float4)*3))) throw(-1);
  //copyInvViewMatrix(invViewMatrix, sizeof(float4)*3);
  imHost.bmp = im2D.map4draw();
  imHost.bind2draw();
  if(CHECK_ERROR(cudaMemcpyToSymbol(im, &imHost, sizeof(imHost)))) throw(-1);
  if(CHECK_ERROR(cudaMemcpyToSymbol(im3D, this, sizeof(im3D_pars)))) throw(-1);
  //if(CHECK_ERROR(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared))) throw(-1);
  switch(mode3D) {
    case 0: render3D <<<dim3(bNx/NW,bNy/NW),dim3(NW,NW)>>>(); break;
    case 1:
#ifdef SURF
      surf_render3D <<<dim3(bNx/NW,bNy/NW),dim3(NW,NW)>>>();
#else//SURF
      printf("Для задействования визуализации на поверхности скомпилируйте im3D.cu с опцией -DSURF или используйте im3Dsurf\n");
#endif//SURF
      break;
    case 2: grad_render3D <<<dim3(bNx/NW,bNy/NW),dim3(NW,NW)>>>(); break;
  }
  if(imHost.draw_flag) {
    //if(mode3D<=1) 
      draw_pal <<<bNx/NW,NW>>>();
    if(mode3D>0) draw_pal3D <<<NW,NW>>>();
    //if(imHost.palDim <= 2) draw_pal <<<bNx/NW,NW>>>();
    //else if(imHost.palDim > 1) draw_pal3D <<<NW,NW>>>();
  } else draw_wavelength_pal <<<bNx/NW,NW>>>();
  if(imHost.negate_flag) negate <<<bNx/NW,NW>>>();
  imHost.nFrame++;
  imHost.unbindAfterDraw();
  im2D.unmapAfterDraw();
} catch(...) {
  printf("recalc3D_im3D: Возникла какая-то ошибка.\n");
}
}

#include <cufft.h>

//inline __device__ float my_fabsC(float2& v) { return v.x;}//hypotf(v.x, v.y); }
inline __device__ float my_fabsC(float2& v) { return hypotf(v.x, v.y); }
inline __device__ int my_abs(int v) { return v>=0?v:-v; }
//inline __device__ int my_abs(int v) { return v==0?1:v>=0?v:-v; }

__global__ void cmplx2abs(cufftComplex *dataC, cufftReal *dataR) {
  //float* pC=(float*)(dataC+blockIdx.x*(blockDim.x/2+1));
  //dataR[blockIdx.x*blockDim.x+threadIdx.x] = pC[threadIdx.x];
  dataR[blockIdx.x*blockDim.x+threadIdx.x] = my_fabsC(dataC[blockIdx.x*(blockDim.x/2+1)+my_abs(blockDim.x/2-threadIdx.x)]);
}
#define CHECK_ERROR_FFT(err) CheckErrorFFT( err, __FILE__,__LINE__)
bool CheckErrorFFT(cufftResult rs, const char *file, int line) {
  if(rs == CUFFT_SUCCESS) return false;
  const char* err="Непонятная ошибка в cuFFT";
  switch(rs) {
  case CUFFT_SUCCESS: err = "0, // The cuFFT operation was successful";
  case CUFFT_INVALID_PLAN: err = "1, // cuFFT was passed an invalid plan handle";
  case CUFFT_ALLOC_FAILED: err = "2, // cuFFT failed to allocate GPU or CPU memory";
  case CUFFT_INVALID_TYPE: err = "3, // No longer used";
  case CUFFT_INVALID_VALUE: err = "4, // User specified an invalid pointer or parameter";
  case CUFFT_INTERNAL_ERROR: err = "5, // Driver or internal cuFFT library error";
  case CUFFT_EXEC_FAILED: err = "6, // Failed to execute an FFT on the GPU";
  case CUFFT_SETUP_FAILED: err = "7, // The cuFFT library failed to initialize";
  case CUFFT_INVALID_SIZE: err = "8, // User specified an invalid transform size";
  case CUFFT_UNALIGNED_DATA: err = "9, // No longer used";
  case CUFFT_INCOMPLETE_PARAMETER_LIST: err = "10, // Missing parameters in call";
  case CUFFT_INVALID_DEVICE: err = "11, // Execution of a plan was on different GPU than plan creation";
  case CUFFT_PARSE_ERROR: err = "12, // Internal plan database error";
  case CUFFT_NO_WORKSPACE: err = "13 // No workspace has been provided prior to plan execution";
  };
  fprintf(stderr, "%s in %s at line %d\n", err, file, line);
  return true;
}
void makeFFTz(float* buf, int Nx, int Ny, int Nz) {
try {
  cufftHandle plan;
  cufftComplex *dataC; cufftReal *dataR;
  if(CHECK_ERROR(cudaMalloc((void**)&dataC, sizeof(cufftComplex)*(Nz/2+1)*Nx*Ny))) throw(-1);
  if(CHECK_ERROR(cudaMalloc((void**)&dataR, sizeof(cufftReal)*Nz*Nx*Ny))) throw(-1);
  if(CHECK_ERROR(cudaMemcpy(dataR, buf, 4*Nz*Nx*Ny, cudaMemcpyHostToDevice))) throw(-1);
  if(CHECK_ERROR_FFT(cufftPlan1d(&plan, Nz, CUFFT_R2C, Nx*Ny))) throw(-1);
  if(CHECK_ERROR_FFT(cufftExecR2C(plan, dataR, dataC))) throw(-1);
  if(CHECK_ERROR(cudaDeviceSynchronize())) throw(-1);
  cmplx2abs <<<Nx*Ny,Nz>>>(dataC, dataR);
  if(CHECK_ERROR(cudaDeviceSynchronize())) throw(-1);
  if(CHECK_ERROR(cudaMemcpy(buf, dataR, 4*Nz*Nx*Ny, cudaMemcpyDeviceToHost))) throw(-1);
  if(CHECK_ERROR_FFT(cufftDestroy(plan))) throw(-1);
  if(CHECK_ERROR(cudaFree(dataC))) throw(-1);
  if(CHECK_ERROR(cudaFree(dataR))) throw(-1);
} catch(...) {
  printf("Ошибка в makeFFTz.\n");
}
}
void im3D_pars::initCuda(Arr3D_pars& arr) {
    //printf("==============\n");
    //for(int ix=0; ix<Nx; ix++) for(int iy=0; iy<Ny; iy++) for(int iz=0; iz<Nz; iz++) arr.Arr3Dbuf[iz*Ny*Nx+iy*Nx+ix]=exp(-0.01*ix);
  // create transfer function texture
  //cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  //if(CHECK_ERROR(cudaMalloc3DArray(&data3D_texArray, &channelDesc, make_cudaExtent(Nx,Ny,Nz)))) throw(-1);
  cudaMemcpy3DParms myparms={0};
  myparms.srcPos = make_cudaPos(0,0,0);
  myparms.dstPos = make_cudaPos(0,0,0);
  myparms.srcPtr = make_cudaPitchedPtr(arr.Arr3Dbuf, Nx*sizeof(floatT4im), Nx, Ny);
  myparms.dstArray = data3D_texArray;
  myparms.extent = make_cudaExtent(Nx,Ny,Nz);
  myparms.kind = arr.inGPUmem?cudaMemcpyDeviceToDevice:cudaMemcpyHostToDevice;
  if(CHECK_ERROR(cudaMemcpy3D(&myparms))) throw(-1);
  //if(draw_edges_flag) draw_edges(imHost.fmax);
  initTex();
}
void im3D_pars::initTex() {
  data3D_tex.normalized = false;//true;
  data3D_tex.filterMode = ((render_type==3)==filterMode_flag)?cudaFilterModeLinear:cudaFilterModePoint; //Point;//filter_pal?cudaFilterModePoint:cudaFilterModeLinear;
  data3D_tex.addressMode[0] = cudaAddressModeClamp;//cyclic_pal?cudaAddressModeWrap:cudaAddressModeClamp;
  data3D_tex.addressMode[1] = cudaAddressModeClamp;//cyclic_pal?cudaAddressModeWrap:cudaAddressModeClamp;
  data3D_tex.addressMode[2] = cudaAddressModeClamp;//cyclic_pal?cudaAddressModeWrap:cudaAddressModeClamp;
  if(CHECK_ERROR(cudaBindTextureToArray(data3D_tex, data3D_texArray))) throw(-1);
}
void im3D_pars::initCuda_surf(Arr3D_pars& arr, size_t sh) {
#ifdef SURF
  cudaMemcpy3DParms myparms={0};
  myparms.srcPos = make_cudaPos(0,0,0);
  myparms.dstPos = make_cudaPos(0,0,0);
  size_t N=Nx; N*=Ny; N*=Nz;
  myparms.srcPtr = make_cudaPitchedPtr(arr.Arr3Dbuf+sh, Nx*sizeof(short2), Nx, Ny);
  myparms.dstArray = data3Dsurf_texArray;
  myparms.extent = make_cudaExtent(Nx,Ny,Nz);
  myparms.kind = arr.inGPUmem?cudaMemcpyDeviceToDevice:cudaMemcpyHostToDevice;
  if(CHECK_ERROR(cudaMemcpy3D(&myparms))) throw(-1);
  initTex_surf();
#endif//SURF
}
void im3D_pars::initTex_surf() {
  data3Dsurf_tex.normalized = false;//true;
  data3Dsurf_tex.filterMode = cudaFilterModePoint; //Point;//filter_pal?cudaFilterModePoint:cudaFilterModeLinear;
  data3Dsurf_tex.addressMode[0] = cudaAddressModeClamp;//cyclic_pal?cudaAddressModeWrap:cudaAddressModeClamp;
  data3Dsurf_tex.addressMode[1] = cudaAddressModeClamp;//cyclic_pal?cudaAddressModeWrap:cudaAddressModeClamp;
  data3Dsurf_tex.addressMode[2] = cudaAddressModeClamp;//cyclic_pal?cudaAddressModeWrap:cudaAddressModeClamp;
  if(CHECK_ERROR(cudaBindTextureToArray(data3Dsurf_tex, data3Dsurf_texArray))) throw(-1);
}
void reset(im3D_pars* p) {
  imHost.reset();
  imHost.set_lim(-1.f,1.f);
  imHost.draw_flag = imHost.negate_flag = imHost.centric_pal = true;
  imHost.cyclic_pal = false;
  if(p) p->reset();
}
void im3D_pars::init3D(Arr3D_pars& arr) {
  //::reset();
  optfid = open(optfName, O_RDWR|O_CREAT, 0644);
  if(optfid<0) printf("Не могу открыть файл %s, сохранение/загрузка наборов опций визуализации невозможна\n", optfName);
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<floatT4im>();
  printf("im3D_pars::init3D: Nx,Ny,Nz=%d,%d,%d\n", Nx,Ny,Nz);
  if(CHECK_ERROR(cudaMalloc3DArray(&data3D_texArray, &channelDesc, make_cudaExtent(Nx,Ny,Nz)))) throw(-1);
  if(CHECK_ERROR(cudaMalloc(&randArr, NW*NW*sizeof(float)))) throw(-1);
  curandState *devStates;
  cudaMalloc( (void **)&devStates, NW*NW*sizeof(curandState) );
  init_rand<<<NW,NW>>>(devStates,randArr);
  if(CHECK_ERROR(cudaDeviceSynchronize())) throw(-1);
  cudaFree(devStates);

  //initCuda(arr); ---- !!!!!!!!!!!!!!!!!!!
#ifdef SURF
  if(CHECK_ERROR(cudaDeviceSynchronize())) throw(-1);
  channelDesc = cudaCreateChannelDesc<short2>();
  if(CHECK_ERROR(cudaMalloc3DArray(&data3Dsurf_texArray, &channelDesc, make_cudaExtent(Nx,Ny,Nz)))) throw(-1);
  if(CHECK_ERROR(cudaDeviceSynchronize())) throw(-1);
  //initCuda_surf(arr); ----- !!!!!!!!!!!!!!!!!!1
#endif//SURF
}
void im3D_pars::recalc_func() {
  if(recalc_always || recalc_at_once) {
    if(recalc_at_once) recalc_at_once=false;
    else xyz->step();
    cudaTimer tm; tm.start();
    if(draw_bmp4backgrownd>=2 && render_type==3) {
      switch(draw_bmp4backgrownd) {
      case 2: recalc_im3D(); break;
      case 3: recalc_sec_im3D(); break;
      }
      save_bmp4backgrownd();
    }
    switch(render_type) {
    case 2: recalc_im3D(); break;
    case 3: recalc3D_im3D(); break;
    }
    runTime=tm.stop(); SmoothFPS = 0.9*SmoothFPS+100./runTime;
    if(type_diag_flag>=2) printf("Frame %d (%.2f/%.2f fps), last run Times: %7.2f msec\n", imHost.nFrame, SmoothFPS, 1000./runTime, runTime);
  }
}
int im3D_pars::init_from_command_line(char** argv) {
  if(strcmp(*argv,"--sensor")==0) { float v[3]; read_float3(v, argv[1]); icalcNdrop.add_sensor(v[0], v[1], v[2]); return 2; }
  return im3D_pars4save::init_from_command_line(argv);
}
floatT4im Arr3D_pars::get_val_from_arr3D(int ix, int iy, int iz) {
  if(inCPUmem) return ((floatT4im*)Arr3Dbuf)[get_ind(ix,iy,iz)];
  floatT4im res;
  if(inGPUmem) CHECK_ERROR(cudaMemcpy(&res, get_ptr((sizeof(floatT4im)/sizeof(float))*ix,iy,iz), sizeof(floatT4im), cudaMemcpyDeviceToHost));
  return res;
}
/*
__global__ void calc_limits(float* buf, float* fLims, int Nxv, int Nxa, int Nxs) {
  float2 fLim;
  float* pf=buf+blockIdx.x*Nxv+threadIdx.x;
  fLim.x = fLim.y = *pf;

  for(int i=0; i<Nxs; i++,pf+=Nxa*Nxv) {
    float v=*pf;
    if(v<fLim.x) fLim.x = v;
    if(v>fLim.y) fLim.y = v;
  }
  __shared__ float2 fLim_sh[Nxv];
  fLim_sh[threadIdx.x] = fLim;
  __syncthreads();
  if(threadIdx.x>warpSize) return;
  for(int i=threadIdx.x; i<Nxv; i+=warpSize) {
    float2 v=fLim_sh[i];
    if(v.x<fLim.x) fLim.x = v.x;
    if(v.y>fLim.y) fLim.y = v.y;
  }
  fLim_sh[threadIdx.x] = fLim;
  if(threadIdx.x>0) return;
  for(int i=0; i<warpSize; i++) {
    float2 v=fLim_sh[i];
    if(v.x<fLim.x) fLim.x = v.x;
    if(v.y>fLim.y) fLim.y = v.y;
  }
  fLims[2*blockIdx.x  ] = fLim.x;
  fLims[2*blockIdx.x+1] = fLim.y;
}

void Arr3D_pars::set_lim_from_arr3D() {
  if(inCPUmem) reset_min_max();
  if(inGPUmem) {
    float* fLims=0,* fLimsD=0;
    CHECK_ERROR(cudaMalloc((void**) &fLimsD, 2*Ny*sizeof(float)));
    calc_limits<<<Ny,Nx>>>(Arr3Dbuf, fLimsD, Nx, Ny, Nz);
    fLims=new float[2*Ny];
    CHECK_ERROR(cudaMemcpy(fLims, fLimsD, 2*Ny*sizeof(float), cudaMemcpyDeviceToHost));
    CHECK_ERROR(cudaFree(fLimsD));
    fMin = fLims[0]; fMax = fLims[1];
    for(int i=0; i<Ny; i++) {
      if(fLims[2*i  ]<fMin) fMin = fLims[2*i  ];
      if(fLims[2*i+1]>fMax) fMax = fLims[2*i+1];
    }
    delete fLims;
  }
}*/
