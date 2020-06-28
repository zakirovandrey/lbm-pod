__device__ inline std::pair<ftype, ftype4> blank_mat(int ix, int iy, int iz){
  const ftype3 r = make_ftype3( (ix-Nx/2), (iy-Ny/2), (iz-Nz/2) );
  ftype rho=1;
  if(length(r)<6.5) rho = 1.25;
  else rho=1;

  //if(ix<Nx/2) rho=1; else rho=1.2;

  //rho=1;

  ftype vx=0.0,vy=0,vz=0,T=PPdev.initial.T0;

  rho*= PPdev.RhoUnitConv;
  vx*= PPdev.VelUnitConv;
  vy*= PPdev.VelUnitConv;
  vz*= PPdev.VelUnitConv;
  T*= PPdev.TempUnitConv;

  return std::make_pair(rho, make_ftype4(vx,vy,vz,T));
};

__device__ inline std::pair<ftype, ftype4> TGV_mat(int ix, int iy, int iz){
  const ftype3 r = make_ftype3( (ix-Nx/2), (iy-Ny/2), (iz-Nz/2) );
  
  const ftype u0 = PPdev.initial.u0;
  const ftype beta0 = PPdev.initial.beta0;
  const ftype r0 = PPdev.initial.r0;
  const ftype T0 = PPdev.initial.T0;

  const ftype CS2 = T0;
  
  const ftype rho0 = 1.0;

  const ftype gamma = double(DIM+2)/DIM;

  ftype vx=0.0,vy=0,vz=0,T=T0,rho=1;
  ftype pressure;

  int plane_TGV=PPdev.initial.planeTVG;

  if(plane_TGV) {
    vx =  u0*sin(2*M_PI*ix/Nx)*cos(2*M_PI*iy/Ny)*cos(2*M_PI*iz/Nz);
    vy = -u0*cos(2*M_PI*ix/Nx)*sin(2*M_PI*iy/Ny)*cos(2*M_PI*iz/Nz);
    vz =0;
    pressure = rho0*(T0 + u0*u0/4*( cos(4.0*M_PI*ix/Nx)+cos(4.0*M_PI*iy/Ny) ) );
  } else {
    vx =  u0*sin(2*M_PI*ix/Nx)*cos(2*M_PI*iy/Ny)*cos(2*M_PI*iz/Nz);
    vy = -u0*cos(2*M_PI*ix/Nx)*sin(2*M_PI*iy/Ny)*cos(2*M_PI*iz/Nz);
    vz =0;
    pressure = rho0*(T0 + u0*u0/16*( cos(4.0*M_PI*ix/Nx)+cos(4.0*M_PI*iy/Ny) )*(cos(4.0*M_PI*iz/Nz)+2.0) );
  }
  rho = pressure/T0;

  rho*= PPdev.RhoUnitConv;
  vx*= PPdev.VelUnitConv;
  vy*= PPdev.VelUnitConv;
  vz*= PPdev.VelUnitConv;
  T*= PPdev.TempUnitConv;

  return std::make_pair(rho, make_ftype4(vx,vy,vz,T));
};


__device__ inline std::pair<ftype, ftype4> vortex_mat(int ix, int iy, int iz){
  const ftype3 r = make_ftype3( (ix-Nx/2), (iy-Ny/2), (iz-Nz/2) );
  
  const ftype udragX = PPdev.initial.uDragX;
  const ftype udragY = PPdev.initial.uDragY;
  const ftype u0 = PPdev.initial.u0;
  const ftype beta0 = PPdev.initial.beta0;
  const ftype r0 = PPdev.initial.r0;
  const ftype T0 = PPdev.initial.T0;

  //Initial conditions from https://www.cfd-online.com/Wiki/2-D_vortex_in_isentropic_flow

  const ftype gamma = double(DIM+2)/DIM;

  ftype T = T0 - (gamma-1)*beta0*beta0/(8*gamma*M_PI*M_PI)*exp(1-dot(r/r0,r/r0) );
  ftype rho = pow( T, 1.0/(gamma-1) );
  ftype pressure = rho*T;

  ftype vx=0.0,vy=0,vz=0;
  vx = udragX - beta0/(2*M_PI)*r.y/r0*exp( (1-dot(r/r0,r/r0))/2 );
  vy = udragY + beta0/(2*M_PI)*r.x/r0*exp( (1-dot(r/r0,r/r0))/2 );

  rho*= PPdev.RhoUnitConv;
  vx*= PPdev.VelUnitConv;
  vy*= PPdev.VelUnitConv;
  vz*= PPdev.VelUnitConv;
  T*= PPdev.TempUnitConv;

  return std::make_pair(rho, make_ftype4(vx,vy,vz,T));
};





