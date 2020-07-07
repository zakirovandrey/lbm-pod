#define SEMI_LAGRANGIAN
#define DIM 2
#define QN 7
#define QN_IN_DIM 28  // QN*2^DIM
constexpr const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=2;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

constexpr const ftype3 _ef[QN] = {
  {0   , 0   , 0 },
  { 1, +sprout::math::sqrt(3.0) , 0 },
  {-1, +sprout::math::sqrt(3.0) , 0 },
  {-2, 0, 0 },
  {-1, -sprout::math::sqrt(3.0) , 0 },
  { 1, -sprout::math::sqrt(3.0) , 0 },
  { 2, 0 , 0 },
};
constexpr const ftype cs2 = TLat;

const ftype W0 = 1./2.;
const ftype W1 = 1./12;
const ftype W2 = 0;

constexpr ftype _w[QN] = {
  W0,
  W1,W1,W1,W1,W1,W1
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0},{1,0,0},{2,0,0},
  {0,1,0},{1,1,0},{0,2,0},
  {2,2,0}
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


