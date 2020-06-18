#ifndef CHECK_ERROR_H
#define CHECK_ERROR_H
#define PRINT_LAST_ERROR() PrintLastError(__FILE__,__LINE__)
void PrintLastError(const char *file, int line);
#define CHECK_ERROR(err) CheckError( err, __FILE__,__LINE__)
bool CheckError( cudaError_t err, const char *file, int line);
#define CHECK_ERROR_DEVICE(err) CheckErrorDevice( err, __FILE__,__LINE__)
bool __device__ CheckErrorDevice( cudaError_t err, const char *file, int line);
void deviceDiagnostics();

struct cudaTimer {
  cudaEvent_t start_event, stop_event;
  cudaTimer() { cudaEventCreate(&start_event); cudaEventCreate(&stop_event); }
  ~cudaTimer() { cudaEventDestroy(start_event); cudaEventDestroy(stop_event); }
  void start() { cudaEventRecord (start_event, 0);/* cudaEventSynchronize (start_event);*/ }
  float stop() {
    float res;
    cudaEventRecord (stop_event, 0); cudaEventSynchronize (stop_event);
    cudaEventElapsedTime (&res, start_event, stop_event);
    return res;
  }
  float restart() { float res=stop(); start(); return res; }
};
inline void read_int2(int2& v, char* str) { v.x = strtol(str, &str, 10); str++; v.y = strtol(str, &str, 10); }
inline void read_int3(int* v, char* str) { for(int i=0; i<3; i++) { v[i] = strtol(str, &str, 10); str++; } }
inline void read_float2(float2& v, char* str) { v.x = strtof(str, &str); str++; v.y = strtof(str, &str); }
inline void read_float3(float* v, char* str) { for(int i=0; i<3; i++) { v[i] = strtof(str, &str); str++; } }
inline float read_float(char* str) { return atof(str); }
#endif//CHECK_ERROR_H
