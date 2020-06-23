__device__ inline std::pair<ftype, ftype4> blank_mat(int ix, int iy, int iz){
  const ftype3 r = make_ftype3( (ix-Nx/2), (iy-Ny/2), (iz-Nz/2) );
  ftype rho=1;
  if(length(r)<25.5) rho = 1.25;
  else rho=1;

  //if(ix<Nx/2) rho=1; else rho=1.2;

  //rho=1;

  ftype vx=0.0,vy=0,vz=0,T=0.1;

  rho*= PPdev.RhoUnitConv;
  vx*= PPdev.VelUnitConv;
  vy*= PPdev.VelUnitConv;
  vz*= PPdev.VelUnitConv;
  T*= PPdev.TempUnitConv;

  return std::make_pair(rho, make_ftype4(vx,vy,vz,T));
};



