#ifndef ARR3D_PARS_HPP
#define ARR3D_PARS_HPP

#ifndef floatT4im
#define floatT4im float
#define DATA_VECTOR_SZ 1
#endif

//------------------------------------
struct Arr3D_pars {
  float* Arr3Dbuf;
  size_t BufSize, ArrSize;
  float fMin, fMax;
  char* fName,* dfName;
  int Nx, Ny, Nz;//размер массива
  int Sx, Sy, Sz;//shrink
  int Nx0, Ny0, Nz0;//размер массива в файле
  bool inGPUmem, inCPUmem;
  Arr3D_pars(): inGPUmem(false), inCPUmem(false), Arr3Dbuf(0), BufSize(0), ArrSize(0) {}
  //int read_from_file(char* fName);
  void clear() { delete Arr3Dbuf; Arr3Dbuf = 0; BufSize = 0; ArrSize = 0; inGPUmem = inCPUmem = false; fMin=0.0; fMax=1.0; }
  inline unsigned int get_ind(int ix, int iy, int iz) { return ix+(iy+iz*Ny)*Nx; }
  inline float* get_ptr(int ix, int iy, int iz) { return Arr3Dbuf+get_ind(ix,iy,iz); }
  void reset_min_max() {
    if(Arr3Dbuf==NULL || Nx*Ny*Nz==0) return;
    fMin = fMax = Arr3Dbuf[0];
    size_t N=Nx; N*=Ny; N*=Nz;
    for(size_t i=0; i<N; i++) {
      float v=Arr3Dbuf[i];
      if(v<fMin) fMin = v;
      if(v>fMax) fMax = v;
    }
    //printf("reset_min_max: fMin=%g, fMax=%g\n", fMin, fMax);
  }
  void reset(int _Nx, int _Ny, int _Nz, int _Sx=1, int _Sy=1, int _Sz=1) {
    Nx0 = _Nx; Ny0 = _Ny; Nz0 = _Nz;
    Sx = _Sx; Sy = _Sy; Sz = _Sz;
    Nx = _Nx/Sx; Ny = _Ny/Sy; Nz = _Nz/Sz;
    fName = 0;
    fMin = -1.0; fMax = 1.0;
    ArrSize = Nx; ArrSize *= Ny; ArrSize *= Nz;
  }
  floatT4im get_val_from_arr3D(int ix, int iy, int iz);
  //void set_lim_from_arr3D();
};

#endif//ARR3D_PARS_HPP
