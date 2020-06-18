__device__ inline std::pair<ftype, ftype4> blank_mat(int ix, int iy, int iz){
  const ftype3 r = make_ftype3(ix-Nx/2, (iy-Ny/2), (iz-Nz/2));
  ftype rho=1;
  if(length(r)<5) rho = 1.2;
  else rho=1;

  //rho=1;

  return std::make_pair(rho, make_ftype4(0,0,0,0.1));
};



