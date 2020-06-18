#define DIM 3
#define QN 125
#define QN_IN_DIM 1000   // QN*2^DIM

// 1 and 3 is ZOT lattice
const int ec1=1;
const int ec2=3;
const int ec1d=ec1*ec1;
const int ec2d=ec2*ec2;

//#include "REFERENCE_TEMP_AND_WEIGHTS-SET-5.h"

//---- Reference Temperature
constexpr const ftype _cs2 = ( -sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) + ec1d + ec2d )/10.;
constexpr const ftype cs2 = _cs2;
const ftype TLat=cs2;

constexpr ftype W0get(const ftype Tv) { return (3*Tv*Tv-(ec1d+ec2d)*Tv+ec1d*ec2d)/(ec1d*ec2d); }
constexpr ftype W1get(const ftype Tv) { return Tv*(ec2d-3*Tv)/(2*ec1d*(ec2d-ec1d)); }
constexpr ftype W2get(const ftype Tv) { return Tv*(ec1d-3*Tv)/(2*ec2d*(ec1d-ec2d)); }

const ftype Tmin = 1./3.;
const ftype Tmax = 3.0;

#define NON_ISOTHERMAL_RELAXATION

constexpr const int3 _e[QN] = {
  {0   , 0   , 0 }, {0   , 0   , +ec1 }, {0   , 0   , -ec1}, {0   , 0   , +ec2}, {0   , 0   , -ec2}, // 0--4
  {0   , +ec1, 0 }, {0   , +ec1, +ec1 }, {0   , +ec1, -ec1}, {0   , +ec1, +ec2}, {0   , +ec1, -ec2}, // 5--9
  {0   , -ec1, 0 }, {0   , -ec1, +ec1 }, {0   , -ec1, -ec1}, {0   , -ec1, +ec2}, {0   , -ec1, -ec2}, // 10--14
  {0   , +ec2, 0 }, {0   , +ec2, +ec1 }, {0   , +ec2, -ec1}, {0   , +ec2, +ec2}, {0   , +ec2, -ec2}, // 15--19
  {0   , -ec2, 0 }, {0   , -ec2, +ec1 }, {0   , -ec2, -ec1}, {0   , -ec2, +ec2}, {0   , -ec2, -ec2}, // 20--24
  {+ec1, 0   , 0 }, {+ec1, 0   , +ec1 }, {+ec1, 0   , -ec1}, {+ec1, 0   , +ec2}, {+ec1, 0   , -ec2}, // 25--29
  {+ec1, +ec1, 0 }, {+ec1, +ec1, +ec1 }, {+ec1, +ec1, -ec1}, {+ec1, +ec1, +ec2}, {+ec1, +ec1, -ec2}, // 30--34
  {+ec1, -ec1, 0 }, {+ec1, -ec1, +ec1 }, {+ec1, -ec1, -ec1}, {+ec1, -ec1, +ec2}, {+ec1, -ec1, -ec2}, // 35--39
  {+ec1, +ec2, 0 }, {+ec1, +ec2, +ec1 }, {+ec1, +ec2, -ec1}, {+ec1, +ec2, +ec2}, {+ec1, +ec2, -ec2}, // 40--44
  {+ec1, -ec2, 0 }, {+ec1, -ec2, +ec1 }, {+ec1, -ec2, -ec1}, {+ec1, -ec2, +ec2}, {+ec1, -ec2, -ec2}, // 45--49
  {-ec1, 0   , 0 }, {-ec1, 0   , +ec1 }, {-ec1, 0   , -ec1}, {-ec1, 0   , +ec2}, {-ec1, 0   , -ec2}, // 50--54
  {-ec1, +ec1, 0 }, {-ec1, +ec1, +ec1 }, {-ec1, +ec1, -ec1}, {-ec1, +ec1, +ec2}, {-ec1, +ec1, -ec2}, // 55--59
  {-ec1, -ec1, 0 }, {-ec1, -ec1, +ec1 }, {-ec1, -ec1, -ec1}, {-ec1, -ec1, +ec2}, {-ec1, -ec1, -ec2}, // 60--64
  {-ec1, +ec2, 0 }, {-ec1, +ec2, +ec1 }, {-ec1, +ec2, -ec1}, {-ec1, +ec2, +ec2}, {-ec1, +ec2, -ec2}, // 65--69
  {-ec1, -ec2, 0 }, {-ec1, -ec2, +ec1 }, {-ec1, -ec2, -ec1}, {-ec1, -ec2, +ec2}, {-ec1, -ec2, -ec2}, // 70--74
  {+ec2, 0   , 0 }, {+ec2, 0   , +ec1 }, {+ec2, 0   , -ec1}, {+ec2, 0   , +ec2}, {+ec2, 0   , -ec2}, // 75--79
  {+ec2, +ec1, 0 }, {+ec2, +ec1, +ec1 }, {+ec2, +ec1, -ec1}, {+ec2, +ec1, +ec2}, {+ec2, +ec1, -ec2}, // 80--84
  {+ec2, -ec1, 0 }, {+ec2, -ec1, +ec1 }, {+ec2, -ec1, -ec1}, {+ec2, -ec1, +ec2}, {+ec2, -ec1, -ec2}, // 85--89
  {+ec2, +ec2, 0 }, {+ec2, +ec2, +ec1 }, {+ec2, +ec2, -ec1}, {+ec2, +ec2, +ec2}, {+ec2, +ec2, -ec2}, // 90--94
  {+ec2, -ec2, 0 }, {+ec2, -ec2, +ec1 }, {+ec2, -ec2, -ec1}, {+ec2, -ec2, +ec2}, {+ec2, -ec2, -ec2}, // 95--99
  {-ec2, 0   , 0 }, {-ec2, 0   , +ec1 }, {-ec2, 0   , -ec1}, {-ec2, 0   , +ec2}, {-ec2, 0   , -ec2}, // 100--104
  {-ec2, +ec1, 0 }, {-ec2, +ec1, +ec1 }, {-ec2, +ec1, -ec1}, {-ec2, +ec1, +ec2}, {-ec2, +ec1, -ec2}, // 105--109
  {-ec2, -ec1, 0 }, {-ec2, -ec1, +ec1 }, {-ec2, -ec1, -ec1}, {-ec2, -ec1, +ec2}, {-ec2, -ec1, -ec2}, // 110--114
  {-ec2, +ec2, 0 }, {-ec2, +ec2, +ec1 }, {-ec2, +ec2, -ec1}, {-ec2, +ec2, +ec2}, {-ec2, +ec2, -ec2}, // 115--119
  {-ec2, -ec2, 0 }, {-ec2, -ec2, +ec1 }, {-ec2, -ec2, -ec1}, {-ec2, -ec2, +ec2}, {-ec2, -ec2, -ec2}, // 120--124
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
  W0*W0*W0, W0*W0*W1, W0*W0*W1, W0*W0*W2, W0*W0*W2,
  W0*W1*W0, W0*W1*W1, W0*W1*W1, W0*W1*W2, W0*W1*W2,
  W0*W1*W0, W0*W1*W1, W0*W1*W1, W0*W1*W2, W0*W1*W2,
  W0*W2*W0, W0*W2*W1, W0*W2*W1, W0*W2*W2, W0*W2*W2,
  W0*W2*W0, W0*W2*W1, W0*W2*W1, W0*W2*W2, W0*W2*W2,
  W1*W0*W0, W1*W0*W1, W1*W0*W1, W1*W0*W2, W1*W0*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W0*W0, W1*W0*W1, W1*W0*W1, W1*W0*W2, W1*W0*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W2*W0*W0, W2*W0*W1, W2*W0*W1, W2*W0*W2, W2*W0*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W0*W0, W2*W0*W1, W2*W0*W1, W2*W0*W2, W2*W0*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  };
  return arr[i];
}
#define TO_SEQ_ELEM_WGET(z, n, text) w_get(n,TLat),
constexpr ftype _w[QN] = {
  BOOST_PP_REPEAT( QN, TO_SEQ_ELEM_WGET, w_get )
};

/*constexpr ftype _w[QN] = {
  W0*W0*W0, W0*W0*W1, W0*W0*W1, W0*W0*W2, W0*W0*W2,
  W0*W1*W0, W0*W1*W1, W0*W1*W1, W0*W1*W2, W0*W1*W2,
  W0*W1*W0, W0*W1*W1, W0*W1*W1, W0*W1*W2, W0*W1*W2,
  W0*W2*W0, W0*W2*W1, W0*W2*W1, W0*W2*W2, W0*W2*W2,
  W0*W2*W0, W0*W2*W1, W0*W2*W1, W0*W2*W2, W0*W2*W2,
  W1*W0*W0, W1*W0*W1, W1*W0*W1, W1*W0*W2, W1*W0*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W0*W0, W1*W0*W1, W1*W0*W1, W1*W0*W2, W1*W0*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W2*W0*W0, W2*W0*W1, W2*W0*W1, W2*W0*W2, W2*W0*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W0*W0, W2*W0*W1, W2*W0*W1, W2*W0*W2, W2*W0*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
};*/


