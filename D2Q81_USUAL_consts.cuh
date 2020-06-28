#define DIM 2
#define QN 81
#define QN_IN_DM 324   // QN*2^DIM

const int ec1=1;
const int ec2=2;
const int ec3=3;
const int ec4=5;
const int ec1d=ec1*ec1;
const int ec2d=ec2*ec2;
const int ec3d=ec3*ec3;
const int ec4d=ec4*ec4;

#include "REFERENCE_TEMP_AND_WEIGHTS-SET-9.h"

const ftype Tmin = 1./100.;
const ftype Tmax = 10.0;

#define NON_ISOTHERMAL_RELAXATION

constexpr const int3 _e[QN] = {
  {0   , 0   , 0 }, {0   , +ec1, 0 }, {0   , -ec1, 0 }, {0   , +ec2, 0 }, {0   , -ec2, 0 }, {0   , +ec3, 0 }, {0   , -ec3, 0 }, {0   , +ec4, 0 }, {0   , -ec4, 0 },
  {+ec1, 0   , 0 }, {+ec1, +ec1, 0 }, {+ec1, -ec1, 0 }, {+ec1, +ec2, 0 }, {+ec1, -ec2, 0 }, {+ec1, +ec3, 0 }, {+ec1, -ec3, 0 }, {+ec1, +ec4, 0 }, {+ec1, -ec4, 0 },
  {-ec1, 0   , 0 }, {-ec1, +ec1, 0 }, {-ec1, -ec1, 0 }, {-ec1, +ec2, 0 }, {-ec1, -ec2, 0 }, {-ec1, +ec3, 0 }, {-ec1, -ec3, 0 }, {-ec1, +ec4, 0 }, {-ec1, -ec4, 0 },
  {+ec2, 0   , 0 }, {+ec2, +ec1, 0 }, {+ec2, -ec1, 0 }, {+ec2, +ec2, 0 }, {+ec2, -ec2, 0 }, {+ec2, +ec3, 0 }, {+ec2, -ec3, 0 }, {+ec2, +ec4, 0 }, {+ec2, -ec4, 0 },
  {-ec2, 0   , 0 }, {-ec2, +ec1, 0 }, {-ec2, -ec1, 0 }, {-ec2, +ec2, 0 }, {-ec2, -ec2, 0 }, {-ec2, +ec3, 0 }, {-ec2, -ec3, 0 }, {-ec2, +ec4, 0 }, {-ec2, -ec4, 0 },
  {+ec3, 0   , 0 }, {+ec3, +ec1, 0 }, {+ec3, -ec1, 0 }, {+ec3, +ec2, 0 }, {+ec3, -ec2, 0 }, {+ec3, +ec3, 0 }, {+ec3, -ec3, 0 }, {+ec3, +ec4, 0 }, {+ec3, -ec4, 0 },
  {-ec3, 0   , 0 }, {-ec3, +ec1, 0 }, {-ec3, -ec1, 0 }, {-ec3, +ec2, 0 }, {-ec3, -ec2, 0 }, {-ec3, +ec3, 0 }, {-ec3, -ec3, 0 }, {-ec3, +ec4, 0 }, {-ec3, -ec4, 0 },
  {+ec4, 0   , 0 }, {+ec4, +ec1, 0 }, {+ec4, -ec1, 0 }, {+ec4, +ec2, 0 }, {+ec4, -ec2, 0 }, {+ec4, +ec3, 0 }, {+ec4, -ec3, 0 }, {+ec4, +ec4, 0 }, {+ec4, -ec4, 0 },
  {-ec4, 0   , 0 }, {-ec4, +ec1, 0 }, {-ec4, -ec1, 0 }, {-ec4, +ec2, 0 }, {-ec4, -ec2, 0 }, {-ec4, +ec3, 0 }, {-ec4, -ec3, 0 }, {-ec4, +ec4, 0 }, {-ec4, -ec4, 0 },
};

#define WEIGHTS_MANUALLY

#define W0 W0get(Tc)
#define W1 W1get(Tc)
#define W2 W2get(Tc)
#define W3 W3get(Tc)
#define W4 W4get(Tc)
#define W5 W5get(Tc)

constexpr ftype w_get(const int i, const ftype Tc) {
  const ftype arr[QN] = {
  W0*W0, W0*W1, W0*W1, W0*W2, W0*W2, W0*W3, W0*W3, W0*W4, W0*W4,
  W1*W0, W1*W1, W1*W1, W1*W2, W1*W2, W1*W3, W1*W3, W1*W4, W1*W4,
  W1*W0, W1*W1, W1*W1, W1*W2, W1*W2, W1*W3, W1*W3, W1*W4, W1*W4,
  W2*W0, W2*W1, W2*W1, W2*W2, W2*W2, W2*W3, W2*W3, W2*W4, W2*W4,
  W2*W0, W2*W1, W2*W1, W2*W2, W2*W2, W2*W3, W2*W3, W2*W4, W2*W4,
  W3*W0, W3*W1, W3*W1, W3*W2, W3*W2, W3*W3, W3*W3, W3*W4, W3*W4,
  W3*W0, W3*W1, W3*W1, W3*W2, W3*W2, W3*W3, W3*W3, W3*W4, W3*W4,
  W4*W0, W4*W1, W4*W1, W4*W2, W4*W2, W4*W3, W4*W3, W4*W4, W4*W4,
  W4*W0, W4*W1, W4*W1, W4*W2, W4*W2, W4*W3, W4*W3, W4*W4, W4*W4,
  };
  return arr[i];
}
#define TO_SEQ_ELEM_WGET(z, n, text) w_get(n,TLat),
constexpr ftype _w[QN] = {
  BOOST_PP_REPEAT( QN, TO_SEQ_ELEM_WGET, w_get )
};

/*constexpr ftype _w[QN] = {
  W0*W0, W0*W1, W0*W1, W0*W2, W0*W2,
  W1*W0, W1*W1, W1*W1, W1*W2, W1*W2,
  W1*W0, W1*W1, W1*W1, W1*W2, W1*W2,
  W2*W0, W2*W1, W2*W1, W2*W2, W2*W2,
  W2*W0, W2*W1, W2*W1, W2*W2, W2*W2,
};*/


