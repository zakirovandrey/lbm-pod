#define SEMI_LAGRANGIAN
#define DIM 2
#define QN 17
#define QN_IN_DIM 68   // QN*2^DIM
constexpr const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=2;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

const ftype ec1 = sprout::math::sqrt((125+5*sprout::math::sqrt(193))/72.0);
const ftype ec1d=ec1*ec1;
constexpr const ftype3 _ef[QN] = {
  {0,0,0},
  {  +ec1, 0     , 0 }, {  -ec1, 0     , 0 }, {0     ,   +ec1, 0 }, {0     ,   -ec1, 0 },
  {  +ec1,   +ec1, 0 }, {  -ec1,   +ec1, 0 }, {  +ec1,   -ec1, 0 }, {  -ec1,   -ec1, 0 },
  {+2*ec1, +2*ec1, 0 }, {-2*ec1, +2*ec1, 0 }, {+2*ec1, -2*ec1, 0 }, {-2*ec1, -2*ec1, 0 },
  {+3*ec1, 0     , 0 }, {-3*ec1, 0     , 0 }, {0     , +3*ec1, 0 }, {0     , -3*ec1, 0 },
};
constexpr const ftype cs2 = TLat;

const ftype W0 = (573 +193*sprout::math::sqrt(193))/8100.0;
const ftype W1 = (3355- 91*sprout::math::sqrt(193))/18000.0;
const ftype W2 = (655+  17*sprout::math::sqrt(193))/27000.0;
const ftype W3 = (685  -49*sprout::math::sqrt(193))/54000.0;
const ftype W4 = (1445-101*sprout::math::sqrt(193))/162000.0;

constexpr ftype _w[QN] = {
  W0,
  W1,W1,W1,W1,
  W2,W2,W2,W2,
  W3,W3,W3,W3,
  W4,W4,W4,W4,
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0},
  {1,0,0},{0,1,0},{2,0,0},{0,2,0},
  {1,1,0},
  {2,1,0},{1,2,0},{3,0,0},{0,3,0},
  {1,3,0},{3,1,0},{2,2,0},{4,0,0},{0,4,0},
  {3,2,0},{2,3,0}
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


