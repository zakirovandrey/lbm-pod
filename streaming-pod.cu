#include "momentsMatrix.cuh"
struct InterpolateStruct{
  ftype3 shifts;
  int3 stencilMinPos;
  __device__ __forceinline__ InterpolateStruct(ftype3 xf, int3 _sminpos): stencilMinPos(_sminpos) {
    shifts = xf - make_ftype3(stencilMinPos);
  }
  template<class T,class F> __device__ inline T calc(F);
};

template<int RegOrder=-1> __global__ __launch_bounds__(LBMconsts::Qn) void  streaming_collision(int ibn) {
  ibn+= blockIdx.x;
  if(ibn>=Nx*Ny*Nz) return;

  const int ix = ibn%Nx;
  const int iy = ibn/Nx%Ny;
  const int iz = ibn/(Nx*Ny);

  const int ild = 0;
  const int ist = 1;

  Cell cell = pars.data.get_cell(ild, ix,iy,iz);

  const int3 ic = make_int3(ix, iy, iz);

  const int3 Nxyz = make_int3(Nx,Ny,Nz);
  using namespace LBMconsts;

  int useSHmem4momMatrix = 0;
  if(sizeof(MomentsMatrix)<48*1024) useSHmem4momMatrix = 1;

  /*if(blockIdx.x==0) {
  MomentsMatrix mtest;
  for(int i=0;i<Qn;i++) for(int j=0;j<Qn;j++) mtest.m[i][j] = 2.+sqrt(abs(sin(100*i)*cos(i*j)));//1.23456+1.23*i/(j+1)+j/(i+1)*0.21;
  MomentsMatrix mtestorig=mtest;
  mtest.inverse();
  ftype mmul[Qn][Qn];
  for(int i=0;i<Qn;i++) for(int j=0;j<Qn;j++) {
    ftype sum=0;
    for(int k=0;k<Qn;k++) sum+= mtestorig.m[i][k]*mtest.m[k][j+Qn];
    if(abs(sum)<1e-8) sum=0;
    printf("%g ",sum);
    if(j==Qn-1) printf("\n");
  }
  assert(0);
  }*/

  MomentsMatrix* Mm,*Mm0;
  __shared__ MomentsMatrix _mmsh;
  __shared__ MomentsMatrix* shMm;
  //__shared__ MomentsMatrix _mmsh0;
  //__shared__ MomentsMatrix* shMm0;

  if(useSHmem4momMatrix) {
    Mm = &_mmsh;
    //Mm0 = &_mmsh0;
  } else {
    if(threadIdx.x==0) {
      shMm = (MomentsMatrix*)malloc(sizeof(MomentsMatrix));
      //shMm0 = (MomentsMatrix*)malloc(sizeof(MomentsMatrix));
      assert(shMm);
      //assert(shMm0);
    }
    __syncthreads();
    Mm=shMm;
    //Mm0=shMm0;
  }

  //cell.vel*=0; cell.T=PPdev.initial.T0/10.0;
  __shared__ Cell cell_new;
  int Niter=0;
  while(Niter<100) {
    ftype T = cell.T;
    ftype rho=cell.rho;
    ftype3 vel = cell.vel;

    ftype4 gauge = make_ftype4(vel.x, vel.y, vel.z, sqrt(T/TLat));

    __syncthreads();
    Mm->init(gauge);
    __syncthreads();
    Mm->inverse();
    __syncthreads();
    /* Mm0->init(make_ftype4(0,0,0,1));
      __syncthreads();
      Mm0->inverse();
      __syncthreads(); */

    const int iq = threadIdx.x;
    ftype3 v = ef[iq]*gauge.w + make_ftype3(gauge.x,gauge.y,gauge.z);

    int3 interpStencilMinPos = ic - make_int3(PPdev.stencilInterpWidth/2);
    const ftype3 xf = make_ftype3(ic)-v;
    if(PPdev.stencilFixed==0) {
      ftype3 pos = xf-make_ftype3(0.5*PPdev.stencilInterpWidth);
      interpStencilMinPos = make_int3( round(pos.x), round(pos.y), round(pos.z) );
    }
    if(DIM<2) interpStencilMinPos.y = ic.y;
    if(DIM<3) interpStencilMinPos.z = ic.z;

    InterpolateStruct interpolate(xf, interpStencilMinPos);
    //cell_new.f[iq] = pars.data.tiles[ild][interpStencilMinPos.x+interpStencilMinPos.y+interpStencilMinPos.z].f[iq];

    if(RegOrder<0) {
      cell_new.f[iq] = interpolate.calc<ftype>( [&] __device__ (int index) {
        ftype mVec[Qn];

        ftype4 igauge = pars.data.tiles[ild][index].uT[0];
        igauge.w = sqrt(igauge.w/TLat);
        calc_moments_vec( igauge, pars.data.tiles[ild][index].f, mVec );

        /*ftype all_fi[Qn]; for(int ii=0;ii<Qn;ii++) all_fi[ii] = Mm->get_inv(iq,mVec);
         calc_moments_vec( make_ftype4(0,0,0,1), all_fi, mVec );
         for(int ii=0;ii<Qn;ii++) if(MomentsPower[ii].x+MomentsPower[ii].y+MomentsPower[ii].z>4) mVec[ii]=0;
         const ftype fi_reg = Mm0->get_inv(iq,mVec);
         return fi_reg;*/

        const ftype fi = Mm->get_inv(iq, mVec);
        return fi;
        //return pars.data.tiles[ild][index].f[iq];
      } );
    } else {
      TensorCoeffs<RegOrder> an = interpolate.calc< TensorCoeffs<RegOrder> >( [&] __device__ (int index) {
        TensorCoeffs<RegOrder> an_p;
      
        ftype4 igauge = pars.data.tiles[ild][index].uT[0];
        igauge.w = sqrt(igauge.w/TLat);
      
        calc_moments_tensors( igauge, pars.data.tiles[ild][index].f, an_p);
        return an_p;
      } );

      TensorCoeffs<RegOrder> dn =  convertAtoD(an, gauge);
      cell_new.f[iq] = eval_fi_Hermit(dn, iq);

      /*ftype4 tmpgauge = pars.data.tiles[ild][174].uT[0];
      tmpgauge.w = sqrt(tmpgauge.w/TLat);
      ftype mVec[Qn];
      calc_moments_vec( tmpgauge, pars.data.tiles[ild][174].f, mVec );
      if(ix==174 && iq==0) printf("Niter=%d moments=(%g %g %g %g %g %g)\n     An=(%g %g %g %g %g %g)\n      Dn=(%g %g %g %g %g %g)\n", 
                                    Niter, mVec[0],mVec[1],mVec[2],mVec[3],mVec[4],mVec[5],
                                    an.k[0],an.k[1],an.k[2],an.k[3],an.k[4],an.k[5],
                                    dn.k[0],dn.k[1],dn.k[2],dn.k[3],dn.k[4],dn.k[5]
                             );*/



    }
    __syncthreads();
    if(threadIdx.x==0) {
      ftype4 Vrho = make_ftype4(0,0,0,0);
      ftype M2 = 0;
      for(int ik=0; ik<Qn; ik++) {
        ftype3 v_k = ef[ik]*gauge.w + make_ftype3(gauge.x,gauge.y,gauge.z);
        Vrho+= make_ftype4(v_k.x,v_k.y,v_k.z,1)*cell_new.f[ik];
        M2+= dot(v_k,v_k)*cell_new.f[ik];
      }
      cell_new.rho = Vrho.w;
      cell_new.vel = make_ftype3(Vrho.x,Vrho.y,Vrho.z)/cell_new.rho;
      cell_new.T = M2/cell_new.rho-dot(cell_new.vel,cell_new.vel); cell_new.T/=DIM;
      if(PPdev.fixedTemperature) cell_new.T=cell.T;
      if(cell_new.T<0) {
        printf("Convergency problem: cell %d %d %d (iteration %d) got negative T=%g, reset to positive\n",
                                          ix,iy,iz,Niter, cell_new.T );
        cell_new.T=-cell_new.T;
      }
    }
    __syncthreads();

    Niter++;
    if( isConv(cell,cell_new) ) { cell=cell_new; break; }
    cell=cell_new;
    __syncthreads();
  }
  //printf("ixyz=%d %d %d Niter=%d\n",ix,iy,iz, Niter);
  __syncthreads();
  if(threadIdx.x==0) {
    if(!useSHmem4momMatrix) free(Mm);

    ftype feq[Qn];
    Cell::calcEq(feq, cell.rho, make_ftype3(0,0,0), TLat);
    collision(cell.f,feq);

    pars.data.set_cell(cell, ist, ix,iy,iz);

    pars.data.tiles[ist][ix+iy*Nx+iz*Nx*Ny].Niter[0] = Niter;
  }
}

inline __device__ ftype LagrPol(int ix,int iy,int iz, const ftype3 shifts, const int N);

template<class Interp_t, class F> __device__ inline Interp_t InterpolateStruct::calc(F func) {
  const int3 Nxyz = make_int3(Nx,Ny,Nz);
  Interp_t val(0);
  const int Npoints = PPdev.stencilInterpWidth+1;
  for(int xs=0; xs<Npoints; xs++) {
    for(int ys=0; ys<((DIM<2)?1:Npoints); ys++) {
      for(int zs=0; zs<((DIM<3)?1:Npoints); zs++) {
        const int3 crd = ( stencilMinPos+make_int3(xs,ys,zs)+Nxyz )%Nxyz;
        const int index = crd.x + crd.y*Nx + crd.z*Nx*Ny;
        const ftype coeff = LagrPol(xs,ys,zs, shifts, Npoints);
        auto Tcoffsp = func(index);
        Tcoffsp*= coeff;
        val+= Tcoffsp;
      }
    }
  }
  return val;
}
inline __device__ ftype LagrPol(int ix,int iy,int iz, const ftype3 shifts, const int N){
  ftype a=1;
  if(DIM>0) for(int ixp=0; ixp<N; ixp++) if(ixp!=ix) a*= (shifts.x-ixp)/(ix-ixp);
  if(DIM>1) for(int iyp=0; iyp<N; iyp++) if(iyp!=iy) a*= (shifts.y-iyp)/(iy-iyp);
  if(DIM>2) for(int izp=0; izp<N; izp++) if(izp!=iz) a*= (shifts.z-izp)/(iz-izp);
  return a;
}
