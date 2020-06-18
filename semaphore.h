#ifndef SEMAPHORE_H
#define SEMAPHORE_H
struct Semaphore{
  int sem;
  __device__ inline int wait(const int lock_val=0){
      int* semval=(int*)&sem;
      volatile int sem_val=0;
      long long int start = clock64();
      int niters=0;
      volatile int* semaddr=&sem;
      while(*semaddr<0) {niters++;}
      return 0;
  }
  __device__ void post(const int val=1){
      __threadfence();
      atomicAdd(&sem,val);
  }
};
#endif
