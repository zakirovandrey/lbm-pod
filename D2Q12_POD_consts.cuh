#define SEMI_LAGRANGIAN
#define DIM 2
#define QN 12
#define QN_IN_DIM 48   // QN*2^DIM
constexpr const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=2;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

const ftype ec1 = sprout::math::sqrt(6);
const ftype ec2 = sprout::math::sqrt((9-3*sprout::math::sqrt(5))/4.0);
const ftype ec3 = sprout::math::sqrt((9+3*sprout::math::sqrt(5))/4.0);
const ftype ec1d=ec1*ec1;
const ftype ec2d=ec2*ec2;
const ftype ec3d=ec3*ec3;
constexpr const ftype3 _ef[QN] = {
  {+ec1, 0   , 0 }, {-ec1, 0   , 0 }, {0   , +ec1, 0 }, {0   , -ec1, 0 },
  {+ec2, +ec2, 0 }, {-ec2, +ec2, 0 }, {+ec2, -ec2, 0 }, {-ec2, -ec2, 0 },
  {+ec3, +ec3, 0 }, {-ec3, +ec3, 0 }, {+ec3, -ec3, 0 }, {-ec3, -ec3, 0 }
};
constexpr const ftype cs2 = TLat;

const ftype W1 = 1./36;
const ftype W2 = (5+2*sprout::math::sqrt(5))/45.0;
const ftype W3 = (5-2*sprout::math::sqrt(5))/45.0;

constexpr ftype _w[QN] = {
  W1,W1,W1,W1,
  W2,W2,W2,W2,
  W3,W3,W3,W3,
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0},
  {1,0,0},{0,1,0},
  {2,0,0},{0,2,0},
  {1,1,0},
  {2,1,0},{1,2,0},
  {3,0,0},{0,3,0},
  {1,3,0},{3,1,0}
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


