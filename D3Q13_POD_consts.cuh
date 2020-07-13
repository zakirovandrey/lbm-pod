#define SEMI_LAGRANGIAN
#define DIM 3
#define QN 13
#define QN_IN_DIM 152   // QN*2^DIM
constexpr const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=3;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

const ftype ec1 = sprout::math::sqrt((5+sprout::math::sqrt(5))/2.0);
const ftype ec2 = sprout::math::sqrt((5-sprout::math::sqrt(5))/2.0);
const ftype ec1d=ec1*ec1;
const ftype ec2d=ec2*ec2;
constexpr const ftype3 _ef[QN] = {
  {0   , 0   , 0 },
  {+ec1, +ec2, 0 }, {+ec1, -ec2, 0 }, {-ec1, +ec2, 0 }, {-ec1, -ec2, 0 }, 
  {0, +ec1, +ec2 }, {0, +ec1, -ec2 }, {0, -ec1, +ec2 }, {0, -ec1, -ec2 }, 
  {+ec2, 0, +ec1 }, {-ec2, 0, +ec1 }, {+ec2, 0, -ec1 }, {-ec2, 0, -ec1 },
};
constexpr const ftype cs2 = TLat;

const ftype W0 = 2./5.;
const ftype W1 = 1./20.;

constexpr ftype _w[QN] = {
  W0,
  W1,W1,W1,W1,W1,W1,W1,W1,W1,W1,W1,W1,
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0},
  {1,0,0},{0,1,0},{0,0,1}, {1,1,0},{0,1,1},{1,0,1},
  {2,0,0},{0,2,0},{0,0,2},
  {2,1,1},{1,2,1},{1,1,2}
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


