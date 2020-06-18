#include "LBMconsts.cuh"
struct Cell{
  static const int Qn=LBMconsts::Qn;
  static const int EqOrder=2;
  ftype f[Qn];
  ftype rho;
  ftype3 vel;
  ftype T;
  int Niter;
  __host__ __device__ static void calcEq(ftype feq[Qn], const ftype Rho, const ftype3 Velocity, const ftype Tempr);
  __host__ __device__ void operator=(const ftype v){
    for(int iq=0;iq<Qn;iq++) f[iq]=v;
  }
};
inline __device__ bool isConv(Cell& c1, Cell& c2) {
  ftype vals1[] = {c1.rho, c1.vel.x, c1.vel.y, c1.vel.z, c1.T};
  ftype vals2[] = {c2.rho, c2.vel.x, c2.vel.y, c2.vel.z, c2.T};
  #ifdef USE_DOUBLE
  const ftype err_abs=1e-12;
  const ftype err_rel=1e-10;
  #else
  const ftype err_abs=1e-6;
  const ftype err_rel=1e-5;
  #endif
  for(int i=0; i<sizeof(vals1)/sizeof(vals1[0]); i++) {
    if(fabs(vals1[i]-vals2[i])>= err_abs+err_rel*fabs(vals1[i])) return false;
  }
  for(int i=0; i<LBMconsts::Qn; i++) {
    if(fabs(c1.f[i]-c2.f[i])>= err_abs+err_rel*fabs(c1.f[i])) return false;
  }
  return true;
}

template<int _Ns> struct Tile_t{
  static const int Ns=_Ns;
  ftype f[Cell::Qn*Ns*Ns*Ns];
  ftype4 uT[Ns*Ns*Ns];
  int Niter[Ns*Ns*Ns];
  __host__ __device__ Cell construct_cell(const int3 loc_crd) {
    const int Ns3 = Ns*Ns*Ns; 
    Cell c;
    for(int iq=0; iq<Cell::Qn; iq++) c.f[iq]=f[iq + (loc_crd.x+loc_crd.y*Ns+loc_crd.z*Ns*Ns)*Ns3 ];
    c.rho=0; for(auto _f: c.f)  c.rho+= _f;
    const int loc_ind = (loc_crd.x+loc_crd.y*Ns+loc_crd.z*Ns*Ns);
    c.vel = make_ftype3( uT[loc_ind].x, uT[loc_ind].y, uT[loc_ind].z );
    c.T = uT[loc_ind].w;
    c.Niter = Niter[loc_ind];

    return c;
  }
};

//typedef Tile_t<2> Tile;
typedef Tile_t<1> Tile;

struct Data_t{
  Tile* tiles[2];
  Tile* tilesHost;
  #ifdef __CUDA_ARCH__
  __device__
  #else
  __host__
  #endif
  Cell get_cell(const int ipar, const int ix, const int iy, const int iz) {
    Tile* ctile = &tiles[ipar][ ix/Tile::Ns + iy/Tile::Ns*(Nx/Tile::Ns) + iz/Tile::Ns*(Nx/Tile::Ns)*(Ny/Tile::Ns) ];
    return ctile->construct_cell( make_int3(ix,iy,iz)%Tile::Ns ); 
  }
  inline __host__ __device__ void set_cell(const Cell& c, const int ipar, const int ix, const int iy, const int iz);
  void malloc_data(const int Nx, const int Ny, const int Nz);
  void swap_ptrs();
  void copyHost2Dev();
  void copyDev2Host(const int);
};

namespace CUDAstep{
//  static constexpr uint Nbsz = 8;
//  static constexpr uint3 Nb = (const uint3){Nbsz,Nbsz,Nbsz};
  static constexpr uint3 Nb = (const uint3){4,4,4};
  static constexpr uint Nblk = Nb.x*Nb.y*Nb.z;
//  static constexpr uint Nblk = 512;//Nb.x*Nb.y*Nb.z;
//  static_assert(Nblk<=Nbsz*Nbsz*Nbsz);
};


#include "data-inl.cu"
