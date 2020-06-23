__host__ __device__ inline void Cell::calcEq(ftype feq[Qn], const ftype Rho, const ftype3 Velocity, const ftype Tempr){
  using namespace LBMconsts;
  const ftype rho = Rho;
  ftype3 u = Velocity;
  if(rho==0) u = make_ftype3(0,0,0);
  ftype dT = dcs2;
  ftype Tcur = cs2;
  const ftype T0=cs2;
  #ifdef NON_ISOTHERMAL_RELAXATION
  Tcur=Tempr;
  #endif

  const ftype dT2 = dT*dT;
  const ftype dT4 = dT2*dT2;
  
  const ftype u2 = dot(u,u);
  const ftype u4 = u2*u2;
  const int TERM1 = (EqOrder>=1);
  const int TERM2 = (EqOrder>=2);
  #ifdef NON_ISOTHERMAL_RELAXATION
  const int TERM3 = 1;//(EqOrder>=3);
  const int TERM4 = 1;//(EqOrder>=4);
  #else
  const int TERM3 = 0;
  const int TERM4 = 0;
  #endif
  const ftype mxwU = 1 - TERM2*u2*0.5*dT;
  for(int i=0; i<Qn; i++) {
    ftype3 eidx = make_ftype3(e[i]);
    const ftype ei2 = dot(eidx,eidx);
    const ftype ei4 = ei2*ei2;
    const ftype eu =  dot(eidx,u);
    const ftype eu2 = eu*eu;
    const ftype eu4 = eu2*eu2;
    ftype mxw  = mxwU +
                 TERM1*eu*dT +
                 TERM2*eu2*0.5*dT2 +
                 TERM2*(Tcur-T0)*0.5*dT*(ei2*dT-DIM)+
                 TERM3*1./6.*eu*dT*( eu2*dT2-3*u2*dT + 3*(Tcur-T0)*dT*(ei2*dT-DIM-2) )+
                 TERM4*1./24.*dT4*( eu4 + 3*T0*T0*u4 - 6*T0*eu2*u2 + 6*(Tcur-T0)*eu2*ei2 + 3*(Tcur-T0)*(Tcur-T0)*ei4
                                    - 6*T0*(Tcur-T0)*(Tcur-T0)*(DIM+2)*ei2 + 3*T0*T0*(Tcur-T0)*(Tcur-T0)*DIM*(DIM+2)
                                    - 6*T0*(Tcur-T0)*u2*ei2 - 6*T0*(Tcur-T0)*(DIM+4)*eu2 + 6*T0*T0*(Tcur-T0)*(DIM+2)*u2
                                  );
    feq[i] = w[i]*rho*mxw;
  }
}
/*__host__ __device__ inline void Cell::calcEq(ftype feq[Qn], const ftype Rho, const ftype3 Velocity, const ftype Tempr){
  using namespace LBMconsts;
  const ftype rho = Rho;
  ftype3 u = Velocity;
  if(rho==0) u = make_ftype3(0,0,0);
  ftype dT = dcs2;
  ftype Tcur = cs2;
  #ifdef NON_ISOTHERMAL_RELAXATION
  Tcur=Tempr
  #endif

  const ftype dT2 = dT*dT;
  const ftype dT3 = dT*dT*dT;
  const ftype dT4 = dT*dT*dT*dT;
  const ftype dT5 = dT*dT*dT*dT*dT;
  
  const ftype u2 = dot(u,u);
  const ftype u4 = u2*u2;
  const int TERM1 = (EqOrder>=1);
  const int TERM2 = (EqOrder>=2);
  const int TERM3 = (EqOrder>=3);
  const int TERM4 = (EqOrder>=4);
  const int TERM5 = (EqOrder>=5);
  const ftype mxwU = 1 - TERM2*u2*0.5*dT + TERM4*u4*0.125*dT2;
  for(int i=0; i<Qn; i++) {
    ftype3 eidx = make_ftype3(e[i]);
    const ftype eu =  dot(eidx,u);
    const ftype eu2 = eu*eu;
    const ftype eu3 = eu*eu*eu;
    const ftype eu4 = eu*eu*eu*eu;
    const ftype eu5 = eu*eu*eu*eu*eu;
    ftype mxw  = mxwU +
                 TERM1*eu*dT +
                 TERM2*eu2*0.5*dT2 +
                 TERM2*(Tempr-cs2)*0.5*dT*(dot(eidx,eidx)*dT-DIM) +
                 TERM3*eu3*ftype(1./6.)*dT3  - TERM3*eu*u2*0.5*dT2 +
                 TERM4*eu4*ftype(1./24)*dT4  - TERM4*eu2*u2*0.25*dT3 +
                 TERM5*eu5*ftype(1./120)*dT5 - TERM5*eu3*u2*ftype(1./12)*dT4 + TERM5*eu*u4*0.125*dT3;
    feq[i] = w[i]*rho*mxw;
    #ifdef NON_ISOTHERMAL_RELAXATION
    feq[i] = w_get(i,Tempr)*rho*mxw;
    #endif
  }
}*/

inline __host__ __device__ void Data_t::set_cell(const Cell& c, const int ipar, const int ix, const int iy, const int iz){
  static_assert(Tile::Ns==1);
  for(int iq=0; iq<Cell::Qn; iq++) {
    const int3 gCrd = make_int3(ix, iy, iz);
    const int gInd =  gCrd.x + gCrd.y*Nx+ gCrd.z*Nx*Ny;

    if(Tile::Ns==1) {
      tiles[ipar][gInd].f[iq] = c.f[iq];
      tiles[ipar][gInd].uT[0] = make_ftype4(c.vel.x, c.vel.y, c.vel.z, c.T);
    } else {
      Tile* ctile = &tiles[ipar][ gCrd.x/Tile::Ns + gCrd.y/Tile::Ns*(Nx/Tile::Ns) + gCrd.z/Tile::Ns*(Nx/Tile::Ns)*(Ny/Tile::Ns) ];
      const int3 intileCrd = gCrd%Tile::Ns;
      const int Ns3 = Tile::Ns*Tile::Ns*Tile::Ns; 
      ctile->f[iq + (intileCrd.x+intileCrd.y*Tile::Ns+intileCrd.z*Tile::Ns*Tile::Ns)*Ns3 ] = c.f[iq];
    }
  }
}
