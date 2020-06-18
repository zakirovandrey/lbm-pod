__forceinline__ __device__ void collision(ftype f[LBMconsts::Qn], ftype feq[LBMconsts::Qn]){
  using namespace LBMconsts;
  register ftype _f[Qn];
  const ftype dtau = PPdev.dtau;
  for(int i=0; i<Qn; i++) _f[i] = f[i]-dtau*(f[i]-feq[i]);
  //for (auto& [a, b] : zip(f, _f)) a=b;
  for(int i=0; i<Qn; i++) f[i] = _f[i];
}
