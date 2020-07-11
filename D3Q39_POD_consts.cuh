#define SEMI_LAGRANGIAN
#define DIM 3
#define QN 39
#define QN_IN_DIM 312   // QN*2^DIM
constexpr const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=3;
const ftype Tmin=0;
const ftype Tmax=TLat*10;

const ftype ec1 = sprout::math::sqrt(3.0/2.0);
const ftype ec2 = 2*ec1;
const ftype ec3 = 3*ec1;
const ftype ec1d=ec1*ec1;
constexpr const ftype3 _ef[QN] = {
  {0   , 0   , 0 },
  {+ec1, 0   , 0 }, {-ec1, 0, 0   }, {0   , +ec1, 0   }, {0   , -ec1, 0   }, {0   , 0, +ec1   }, {0   , 0, -ec1   },
  {+ec1, +ec1, +ec1}, {-ec1, +ec1, +ec1}, {+ec1, -ec1, +ec1}, {+ec1, +ec1, -ec1}, {-ec1, -ec1, +ec1}, {-ec1, +ec1, -ec1}, {+ec1, -ec1, -ec1}, {-ec1, -ec1, -ec1},
  {+ec2, 0   , 0 }, {-ec2, 0, 0   }, {0   , +ec2, 0   }, {0   , -ec2, 0   }, {0   , 0, +ec2   }, {0   , 0, -ec2   },
  {+ec2, +ec2, 0 }, {+ec2, 0, +ec2}, {0, +ec2, +ec2 },
  {+ec2, -ec2, 0 }, {+ec2, 0, -ec2}, {0, +ec2, -ec2 },
  {-ec2, +ec2, 0 }, {-ec2, 0, +ec2}, {0, -ec2, +ec2 },
  {-ec2, -ec2, 0 }, {-ec2, 0, -ec2}, {0, -ec2, -ec2 },
  {+ec3, 0   , 0 }, {-ec3, 0, 0   }, {0   , +ec3, 0   }, {0   , -ec3, 0   }, {0   , 0, +ec3   }, {0   , 0, -ec3   }
};
constexpr const ftype cs2 = TLat;

const ftype W0 = 1.0/12.0;
const ftype W1 = 1.0/12.0;
const ftype W2 = 1.0/27.0;
const ftype W3 = 2.0/135.0;
const ftype W4 = 1.0/432.0;
const ftype W5 = 1.0/1620.0;

constexpr ftype _w[QN] = {
  W0,
  W1,W1,W1,W1,W1,W1,
  W2,W2,W2,W2,W2,W2,W2,W2,
  W3,W3,W3,W3,W3,W3,
  W4,W4,W4,W4,W4,W4,W4,W4,W4,W4,W4,W4,
  W5,W5,W5,W5,W5,W5
};

constexpr const int3 _MomentsPower[QN] = {
  {0,0,0},
  {1,0,0},{0,1,0},{0,0,1},
  {2,0,0},{0,2,0},{0,0,2},
  {1,1,0},{1,0,1},{0,1,1},
  {1,1,1},
  {2,1,1},{1,2,1},{1,1,2},
  {2,2,1},{1,2,2},{2,1,2},
  {2,2,2},
  {2,2,0},{2,0,2},{0,2,2},
  {2,1,0},{2,0,1},{0,2,1},{1,2,0},{0,1,2},{1,0,2},
  {3,0,0},{0,3,0},{0,0,3},
  {4,0,0},{0,4,0},{0,0,4},
  {3,1,0},{1,3,0},{3,0,1},{1,0,3},{0,1,3},{0,3,1}
/*
  {0,0,0},{1,0,0},{2,0,0}, {0,1,0},{1,1,0},{2,1,0}, {0,2,0},{1,2,0},{2,2,0},
  {0,0,1},{1,0,1},{2,0,1}, {0,1,1},{1,1,1},{2,1,1}, {0,2,1},{1,2,1},{2,2,1},
  {0,0,2},{1,0,2},{2,0,2}, {0,1,2},{1,1,2},{2,1,2}, {0,2,2},{1,2,2},{2,2,2},*/
};

HOST_DEV_CHOOSE const int3 MomentsPower[QN] = { SEQ_LIST_COLLECTION(_MomentsPower,QN) };
HOST_DEV_CHOOSE const ftype3 ef[QN] = { SEQ_LIST_COLLECTION(_ef,QN) };
HOST_DEV_CHOOSE const int3 e[QN] = {0};


