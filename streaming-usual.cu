
__global__ __launch_bounds__(LBMconsts::Qn) void  streaming_collision(int ibn) {
  if(threadIdx.x!=0) return;

  /*const int ix = blockIdx.x*CUDAstep::Nb.x + threadIdx.x;
  const int iy = blockIdx.y*CUDAstep::Nb.y + threadIdx.y;
  const int iz = blockIdx.z*CUDAstep::Nb.z + threadIdx.z;*/

  ibn+= blockIdx.x;
  if(ibn>=Nx*Ny*Nz) return;

  const int ix = ibn%Nx;
  const int iy = ibn/Nx%Ny;
  const int iz = ibn/(Nx*Ny);

  const int ild = 0;
  const int ist = 1;

  Cell cell = pars.data.get_cell(ild, ix,iy,iz);
  ftype rho=cell.rho;
  ftype3 vel = cell.vel;

  const int3 ic = make_int3(ix, iy, iz);

  const int3 Nxyz = make_int3(Nx,Ny,Nz);
  using namespace LBMconsts;

  ftype fnew[Qn], feq[Qn];
  ftype4 Vrho = make_ftype4(0,0,0,0);
  ftype T=0;

  for(int iq=0; iq<Qn; iq++) {
    const int3 icn = (ic-e[iq]+Nxyz)%Nxyz;
    const int nind = icn.x + icn.y*Nx + icn.z*Nx*Ny;
    fnew[iq] = pars.data.tiles[ild][nind].f[iq];
    Vrho+= make_ftype4(e[iq].x,e[iq].y,e[iq].z,1)*fnew[iq];
    T+= dot(e[iq],e[iq])*fnew[iq];
  }

  const ftype3 Vel = make_ftype3(Vrho.x,Vrho.y,Vrho.z)/Vrho.w;
  T = T/Vrho.w - dot(Vel,Vel);
  T/= DIM;

  Cell::calcEq(feq, Vrho.w, Vel, T);
  collision(fnew,feq);

  cell.rho = 0;
  cell.vel = make_ftype3(0,0,0);
  cell.T = 0;
  for(int iq=0; iq<Qn; iq++) {
    cell.f[iq] = fnew[iq];
    cell.rho+= fnew[iq];
    ftype3 ef = make_ftype3(e[iq]);
    cell.vel+= ef*fnew[iq];
    cell.T+= dot(ef,ef)*fnew[iq];
  }
  cell.vel/= cell.rho;
  cell.T = cell.T/cell.rho - dot(cell.vel,cell.vel);
  cell.T/= DIM;

  pars.data.set_cell(cell, ist, ix,iy,iz);
}
