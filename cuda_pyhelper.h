#ifndef CUDA_PYHELPER
#define CUDA_PYHELPER
#include <stdio.h>
#include <cuda.h>

#define CHECK_ERROR(err) CheckError( err, __FILE__,__LINE__)
bool CheckError(cudaError_t err, const char *file, int line) {
  if(err==cudaSuccess) return false;
  fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err), file, line);
  return true;
}


int getBusID(int device){
  cudaDeviceProp devProp;
  CHECK_ERROR( cudaGetDeviceProperties(&devProp, device) );
  return devProp.pciBusID;
}
//struct cudafuncs{
//  void cudaGetDeviceProperties(cudaDeviceProp* devProp, int device){
//    CHECK_ERROR( cudaGetDeviceProperties(devProp, device) );
//  }
//} cf;

#endif
