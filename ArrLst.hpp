#include <semaphore.h>
#include <pthread.h>
#include <stdio.h>
#include "Arr3Dpars.hpp"

struct Arr3DList {
  sem_t sem_load, sem_drop;
  char** fNames;
  int Narr, Iarr, arr4load, shrink[3], size[3];
  size_t BufSumSize;
  Arr3D_pars* pars_arr;
  Arr3DList(int _Narr, char** _fNames, int _size[3], int _shrink[3]): Narr(_Narr), fNames(_fNames), arr4load(0), BufSumSize(0) {
    for(int i=0; i<3; i++) { shrink[i] = _shrink[i]; size[i] = _size[i]; }
    sem_init(&sem_load, 0,0);
    sem_init(&sem_drop, 0,Narr);
    pars_arr = new Arr3D_pars[Narr];
  }
  ~Arr3DList() {
    for(int s=0; s<Narr; s++) pars_arr[s].clear();
    delete [] pars_arr;
  }
  void arr_clear_gpu(int iarr);
  void arr2gpu(int iarr);
  bool first2gpu();
  bool next_arr2gpu() {
    if(Iarr>=Narr-1) { printf("Больше массивов нет, %s последний из %d\n", fNames[Narr-1], Narr); return false; }
    sem_wait(&sem_load);
    arr_clear_gpu(Iarr);
    arr2gpu(++Iarr);
    return true;
  }
  bool prev_arr2gpu() {
    if(Iarr<=0) { printf("Больше массивов нет, %s первый из %d\n", fNames[0], Narr); return false; }
    arr_clear_gpu(Iarr);
    arr2gpu(--Iarr);
    sem_post(&sem_load);
    return true;
  }
  bool diff_arr2gpu();
  int read_next_file();
};
extern Arr3DList* arr3D_list;

#include "cubeLR.hpp"
