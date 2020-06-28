#define SEMI_LAGRANGIAN
#define DIM 3
#define QN 15
#define QN_IN_DIM 120   // QN*2^DIM
constexpr const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=3;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

//strange_coeff = 0.8883588741731
const ftype ec1=0.8883588741731*sprout::math::sqrt(5.0+sprout::math::sqrt(10.));
const ftype ec2=0.8883588741731*sprout::math::sqrt(5.0-sprout::math::sqrt(10.));
const ftype ec1d=ec1*ec1;
const ftype ec2d=ec2*ec2;
constexpr const ftype3 _ef[QN] = {
  {0   , 0   , 0 },
  {+ec1, 0   , 0    }, {-ec1, 0   , 0   }, {0   , +ec1, 0   }, {0   , -ec1, 0   }, {0   , 0, +ec1   }, {0   , 0, -ec1   },
  {+ec2, +ec2, +ec2 }, {-ec2, +ec2, +ec2 }, {+ec2, -ec2, +ec2 }, {+ec2, +ec2, -ec2 }, {-ec2, -ec2, +ec2 }, {-ec2, +ec2, -ec2 }, {+ec2, -ec2, -ec2 }, {-ec2, -ec2, -ec2 },
};
constexpr const ftype cs2 = TLat;

const ftype W2 = 1.0/(8*ec2d*ec2d);
const ftype W1 = (1-8*ec2d*W2)/(2*ec1d);
const ftype W0 = 1-6*W1-8*W2;
constexpr ftype _w[QN] = {
  W0,
  W1,W1,W1,W1,W1,W1,
  W2,W2,W2,W2,W2,W2,W2,W2
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0},
  {1,0,0},{0,1,0},{0,0,1},
  {2,0,0},{0,2,0},{0,0,2},
  {1,1,1},{2,2,2},
  {2,2,1},{2,1,2},{1,2,2},
  {2,1,1},{1,2,1},{1,1,2}
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


