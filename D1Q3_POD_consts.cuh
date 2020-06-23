#define SEMI_LAGRANGIAN
#define DIM 1
#define QN 3
#define QN_IN_DIM 6   // QN*2^DIM
constexpr const ftype TLat=1;
const ftype Dim=1;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

const ftype ec1=1.7320508075688772;
const ftype ec1d=ec1*ec1;
constexpr const ftype3 _ef[QN] = {
  {0   , 0   , 0 }, {+ec1, 0   , 0 }, {-ec1, 0   , 0 }
};
constexpr const ftype cs2 = TLat;
const ftype W0=2./3.;
const ftype W1=1./6.;
constexpr ftype _w[QN] = {
  W0, W1, W1
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0},{1,0,0},{2,0,0}
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


