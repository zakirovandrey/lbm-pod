
template<int RO=-1> __global__ __launch_bounds__(LBMconsts::Qn) void  streaming_collision(int ibn) {
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

  ftype3 Qm=make_ftype3(0,0,0);
  ftype3 Qp=make_ftype3(0,0,0);
  const ftype R=1;
  const ftype lx=1,ly=1,lz=1;
  Cell ncell;
  ncell = pars.data.get_cell(ild, (ix-1+Nx)%Nx,iy,iz); Qm.x = ncell.rho*ncell.vel.x*(lx*lx - 3*R*ncell.T-ncell.vel.x*ncell.vel.x);
  ncell = pars.data.get_cell(ild, (ix+1   )%Nx,iy,iz); Qp.x = ncell.rho*ncell.vel.x*(lx*lx - 3*R*ncell.T-ncell.vel.x*ncell.vel.x);
  ncell = pars.data.get_cell(ild, ix,(iy-1+Ny)%Ny,iz); Qm.y = ncell.rho*ncell.vel.y*(ly*ly - 3*R*ncell.T-ncell.vel.y*ncell.vel.y);
  ncell = pars.data.get_cell(ild, ix,(iy+1   )%Ny,iz); Qp.y = ncell.rho*ncell.vel.y*(ly*ly - 3*R*ncell.T-ncell.vel.y*ncell.vel.y);
  ncell = pars.data.get_cell(ild, ix,iy,(iz-1+Nz)%Nz); Qm.z = ncell.rho*ncell.vel.z*(lz*lz - 3*R*ncell.T-ncell.vel.z*ncell.vel.z);
  ncell = pars.data.get_cell(ild, ix,iy,(iz+1   )%Nz); Qp.z = ncell.rho*ncell.vel.z*(lz*lz - 3*R*ncell.T-ncell.vel.z*ncell.vel.z);

  const ftype3 difQ = 0.5*(Qp-Qm);

  Cell::calcEq(feq, Vrho.w, Vel, T, difQ);
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
