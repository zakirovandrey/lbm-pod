#define DIM 1
#define QN 5
#define QN_IN_DIM 10   // QN*2^DIM

// 1 and 3 is ZOT lattice
const int ec1=1;
const int ec2=3;
const int ec1d=ec1*ec1;
const int ec2d=ec2*ec2;

//#include "REFERENCE_TEMP_AND_WEIGHTS-SET-5.h"

//---- Reference Temperature
constexpr const ftype _cs2 = ( -sprout::math::sqrt(double(ec1d*ec1d+ec2d*ec2d-14./3.*ec1d*ec2d)) + ec1d + ec2d )/10.;
constexpr const ftype cs2 = _cs2;
const ftype TLat=cs2;

constexpr ftype W0get(const ftype Tv) { return (3*Tv*Tv-(ec1d+ec2d)*Tv+ec1d*ec2d)/(ec1d*ec2d); }
constexpr ftype W1get(const ftype Tv) { return Tv*(ec2d-3*Tv)/(2*ec1d*(ec2d-ec1d)); }
constexpr ftype W2get(const ftype Tv) { return Tv*(ec1d-3*Tv)/(2*ec2d*(ec1d-ec2d)); }
/*constexpr ftype W0get(const ftype Tv) { return (-3.*ec1d*ec1d-3*ec2d*ec2d+54*ec1d*ec2d+(ec1d+ec2d)*sprout::math::sqrt(9.0*ec1d*ec1d+9*ec2d*ec2d-42*ec1d*ec2d))/(75.0*ec1d*ec2d); }
constexpr ftype W1get(const ftype Tv) { return (9.*ec1d*ec1d-6*ec2d*ec2d-27*ec1d*ec2d-(3*ec1d-2*ec2d)*sprout::math::sqrt(9.0*ec1d*ec1d+9*ec2d*ec2d-42*ec1d*ec2d))/(300.*ec1d*(ec1d-ec2d)); }
constexpr ftype W2get(const ftype Tv) { return (9.*ec2d*ec2d-6*ec1d*ec1d-27*ec1d*ec2d-(3*ec2d-2*ec1d)*sprout::math::sqrt(9.0*ec1d*ec1d+9*ec2d*ec2d-42*ec1d*ec2d))/(300.*ec2d*(ec2d-ec1d)); }*/

const ftype Tmin = 1./3.;
const ftype Tmax = 3.0;

#define NON_ISOTHERMAL_RELAXATION

constexpr const int3 _e[QN] = {
  {0   , 0   , 0 }, {+ec1, 0   , 0 }, {-ec1, 0   , 0 }, {+ec2, 0   , 0 }, {-ec2, 0   , 0 },
};

#define WEIGHTS_MANUALLY
/*constexpr ftype _w[QN] = {
  W0, W1, W1, W2, W2
};*/

#define W0 W0get(Tc)
#define W1 W1get(Tc)
#define W2 W2get(Tc)
#define W3 W3get(Tc)
#define W4 W4get(Tc)
#define W5 W5get(Tc)

constexpr ftype w_get(const int i, const ftype Tc) {
  const ftype arr[QN] = {W0, W1, W1, W2, W2};
  return arr[i];
}
#define TO_SEQ_ELEM_WGET(z, n, text) w_get(n,TLat),
constexpr ftype _w[QN] = {
  BOOST_PP_REPEAT( QN, TO_SEQ_ELEM_WGET, w_get )
};

