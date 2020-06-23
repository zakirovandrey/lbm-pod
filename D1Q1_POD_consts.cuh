#define SEMI_LAGRANGIAN
#define DIM 1
#define QN 1
#define QN_IN_DIM 2   // QN*2^DIM
constexpr const ftype TLat=1;
const ftype Dim=1;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

constexpr const ftype3 _ef[QN] = {
  {0   , 0   , 0 }
};
constexpr const ftype cs2 = TLat;
const ftype W0=1;
constexpr ftype _w[QN] = {
  W0
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0}
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


