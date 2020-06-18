#include <cuda.h>
#include <stdio.h>
#include "err.h"
namespace errors{
__managed__ cudaError_t last_err=cudaSuccess;
};

void PrintLastError(const char *file, int line) {
  cudaError_t err=cudaGetLastError();
  if(err!=cudaSuccess) fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err), file, line);
}
bool CheckError(cudaError_t err, const char *file, int line) {
  cudaError_t dev_err; cudaMemcpy(&dev_err, &errors::last_err, sizeof(cudaError_t), cudaMemcpyDefault);
  if(err==cudaSuccess && dev_err==cudaSuccess) return false;
  if(err==cudaSuccess) err = errors::last_err;
  fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err), file, line);
  return true;
}
bool __device__ CheckErrorDevice(cudaError_t err, const char *file, int line) {
  if(err==cudaSuccess) return false;
  atomicCAS((int*)(&errors::last_err), cudaSuccess, err);
  return true;
}

void deviceDiagnostics(){
  int deviceCount;
  CHECK_ERROR( cudaGetDeviceCount(&deviceCount) );  
  printf("GPU devices :: %d \n", deviceCount);
  cudaDeviceProp devProp[deviceCount];
  for(int i = 0; i < deviceCount; ++i) {
    printf("*** CUDA Device #%d ***", i);
    CHECK_ERROR( cudaGetDeviceProperties(&devProp[i], i) );
    printf("%s ***\n", devProp[i].name);
    printf("\t%d.%d compute capability\n", devProp[i].major, devProp[i].minor);
    printf("\t%d multiprocessors\n", devProp[i].multiProcessorCount);
    printf("\t%.2fGB max mem pitch of %.2fGB global memory\n", devProp[i].memPitch/(1024.*1024.*1024), devProp[i].totalGlobalMem/(1024.*1024.*1024));
    printf("\t%.2fKB total shared memory per block\n", devProp[i].sharedMemPerBlock/1024.);
    printf("\t%.2fKB total constant memory\n", devProp[i].totalConstMem/1024.);
    printf("\t%.2fK registers per block\n", devProp[i].regsPerBlock/1024.);
    printf("\t%d/%d threads per Warp/block\n", devProp[i].warpSize, devProp[i].maxThreadsPerBlock);
    printf("\tClock rate: %.2fGHz\n", devProp[i].clockRate*1e-6);
    printf("\tTexture alignment: %luB\n", devProp[i].textureAlignment);
    printf("\tConcurrent copy and execution: %s\n", (devProp[i].deviceOverlap ? "Yes" : "No"));
    printf("\tKernel execution timeout: %s\n", (devProp[i].kernelExecTimeoutEnabled ? "Yes" : "No"));
  }
}
#include "im3D.hpp"
extern bool recalc_always; 
int im3D_pars4save::init_from_command_line(char** argv) {
  if(strncmp(*argv,"--help",6)==0) return -1;
  if(strncmp(*argv,"--devQ",6)==0) { deviceDiagnostics(); return 1; }
  if(strcmp(*argv,"--box")==0) read_float3(BoxFactor, argv[1]);
  else if(strcmp(*argv,"--load")==0) load_from_file(argv[1]);
  else if(strcmp(*argv,"--mesh")==0) read_float3(MeshBox, argv[1]);
  else if(strcmp(*argv,"--sh_mesh")==0) read_float3(MeshShift, argv[1]);
  else if(strcmp(*argv,"--Dmesh")==0) Dmesh=read_float(argv[1]);
  else if(strcmp(*argv,"--zoom")==0) read_float3(Dzoom, argv[1]);
  else if(strcmp(*argv,"--add")==0) read_float3(Dadd, argv[1]);
  else if(strcmp(*argv,"--shrink")==0) read_int3(Dshrink, argv[1]);
  else if(strcmp(*argv,"--Narr")==0) read_int3(Narr, argv[1]);
  else if(strcmp(*argv,"--step")==0) read_float3(step, argv[1]);
  else if(strcmp(*argv,"--base")==0) read_float3(base, argv[1]);
  else if(strcmp(*argv,"--bkgr_col")==0) read_float3(bkgr_col, argv[1]);
  else if(strcmp(*argv,"--mesh_col")==0) read_float3(mesh_col, argv[1]);
  else if(strcmp(*argv,"--box_col")==0) read_float3(box_col, argv[1]);
  else if(strcmp(*argv,"--rot_point")==0) read_float3(RotPoint, argv[1]);
  else if(strcmp(*argv,"--box_shrink")==0) read_float3(box_shrink, argv[1]);
  else if(strcmp(*argv,"--drop_dir")==0) strcpy(drop_dir,argv[1]);
  else if(strcmp(*argv,"--cntr")==0) cntr_levels[cntr_num++]=read_float(argv[1]);
  else if(strcmp(*argv,"--ld_sz")==0) read_int2(ld_sz, argv[1]);
  else if(strcmp(*argv,"--recalc_always")==0) { recalc_always=true; return 1; }
  else if(strcmp(*argv,"--cntr_clear")==0) { cntr_num=0; return 1; }
  else if(strcmp(*argv,"--nocomp")==0) return 1;
  else if(strcmp(*argv,"--norun")==0) return 1;
  else if(strcmp(*argv,"--redefine")==0) return 2;
  else { printf("Illegal parameters' syntax notation\n<%s>", *argv); return 0; }
  //else if(strcmp(*argv,"--")==0) read_float3(, argv[1]);
  //printf("par: %s; vals: %s\n", argv[0], argv[1]);
  return 2;
}
const char* im3D_pars4save::command_line_help_string() {
  return "[--devQ] [--load <opt-file>] [--zoom \"1. 1. 1.\"] [--shrink \"1 1 1\"] [--step \"1. 1. 1.\"] [--base \"1. 1. 1.\"] [--box \"1. 1. 1.\"] [--mesh \"200. 200. 200.\"] [--Dmesh 5.] [--drop_dir \".\"] [--bkgr_col \"0.1 0.1 0.1\"] [--mesh_col \"0.8 0.8 0.2\"] [--box_col \"1. 1. 1.\"] [--box_shrink \"1. 1. 1.\"] [--sensor \"1 1 1\"]";
}
void im3D_pars4save::print_command_line_help() {
  printf("  --devQ\tВыдаёт информацию о видеокартах на компьютере;\n");
  printf("  --load\tВводит параметры из файла <opt-file>, сохранённые ранее клавишей <w/W>\n");
  printf("  --zoom\tмасштабный фактор, действует на 2D режим и размер окна, [1. 1. 1.];\n");
  printf("  --add \tдобавляет пространство к размеру окна. Требуется для вывода 3D на фоне 2D, [0. 0. 0.];\n");
  printf("  --shrink\tмасштабный фактор, действует везде, сокращает требования к памяти, [1 1 1];\n");
  printf("  --Narr\tявно заданный размер массива (если =0, берётся из первого файла) [0 0 0];\n");
  printf("  --box \tкоррекция пропорций размера бокса в 3D режиме, [1. 1. 1.];\n");
  printf("  --step \tшаги между точками, действует только на тики, [1. 1. 1.];\n");
  printf("  --base \tшаги между точками, действует только на тики, [0. 0. 0.];\n");
  printf("  --mesh\tрасстояние между линиями сетки в боксе по координатам в ячейках (до коррекции), [100. 100. 100.];\n");
  printf("  --sh_mesh\tсдвиг линий сетки в боксе по координатам в ячейках (до коррекции), [0. 0. 0.];\n");
  printf("  --Dmesh\tширина линии сетки в пикселях (со сглаживанием выглядит несколько уже), [5.];\n");
  printf("  --drop_dir\tимя директории, в которую будут сохраняться различные файлы, [.];\n");
  printf("  --bkgr_col\tцвет фона, [0.1 0.1 0.1];\n");
  printf("  --mesh_col\tцвет линий сетки, [0.8 0.8 0.2];\n");
  printf("  --box_col\tцвет линий бокса, [1.0 1.0 1.0];\n");
  printf("  --box_shrink\t коэффициент растяжения размеров бокса, [1.0 1.0 1.0];\n");
  printf("  --rot_point\t точка в боксе, относительно которой производится вращение, [0.5 0.5 0.5];\n");
  printf("  --sensor\tкоординаты сенсора, можно задавать несколько сенсоров;\n");
  printf("  --cntr\tзначение уровня контура, можно задавать несколько уровней;\n");
  printf("  --cntr_clear\tочищает все ранее заданные значения уровней контура;\n");
  printf("  --ld_sz\tчтение сохраненных ранее параметров в режиме совместимости, [80 288];\n");
}
