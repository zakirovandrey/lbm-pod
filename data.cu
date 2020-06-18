#include "data.cuh"

void Data_t::malloc_data(const int Nx, const int Ny, const int Nz){
  const int Ns3 = Tile::Ns*Tile::Ns*Tile::Ns; 
  const size_t sz = long(Nx)*Ny*Nz/Ns3*sizeof(Tile);
  printf("Total data size = %g GB\n",double(sz)/1024/1024/1024); 
  for(int ipar:{0,1} ) CHECK_ERROR( cudaMalloc((void**)&tiles[ipar], sz ) );
  CHECK_ERROR( cudaMallocHost((void**)&tilesHost, sz ) );
  for(auto itile: tiles ) CHECK_ERROR( cudaMemset(itile, 0, sz ) );
  CHECK_ERROR( cudaMemset(tilesHost, 0, sz ) );
};
void Data_t::copyHost2Dev(){
  const int Ns3 = Tile::Ns*Tile::Ns*Tile::Ns; 
  const size_t sz = long(Nx)*Ny*Nz/Ns3*sizeof(Tile);
  for(auto itiles: tiles ) CHECK_ERROR( cudaMemcpy(itiles, tilesHost, sz, cudaMemcpyHostToDevice ) );
}
void Data_t::copyDev2Host(const int ipar){
  const int Ns3 = Tile::Ns*Tile::Ns*Tile::Ns; 
  const size_t sz = long(Nx)*Ny*Nz/Ns3*sizeof(Tile);
  CHECK_ERROR( cudaMemcpy(tilesHost, tiles[ipar], sz, cudaMemcpyDeviceToHost ) );
}

void Data_t::swap_ptrs(){
  Tile* tmpptr = tiles[0];
  tiles[0] = tiles[1];
  tiles[1] = tmpptr;
}

