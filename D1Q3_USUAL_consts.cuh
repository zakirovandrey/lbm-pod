#define DIM 1
#define QN 3
#define QN_IN_DIM 6   // QN*2^DIM

const int ec1=1;
const int ec1d=ec1*ec1;

constexpr const ftype _cs2 = 1./3.;
constexpr const ftype cs2 = _cs2;
const ftype TLat=cs2;

constexpr ftype W0get(const ftype Tv) { return 2./3.; }
constexpr ftype W1get(const ftype Tv) { return 1./6.; }
constexpr ftype W2get(const ftype Tv) { return 0; }

const ftype Tmin = 0;
const ftype Tmax = 100.0;


//#define NON_ISOTHERMAL_RELAXATION

constexpr const int3 _e[QN] = {
  {0   , 0   , 0 }, {+ec1, 0   , 0 }, {-ec1, 0   , 0 }
};

#define WEIGHTS_MANUALLY

#define W0 W0get(Tc)
#define W1 W1get(Tc)
#define W2 W2get(Tc)
#define W3 W3get(Tc)
#define W4 W4get(Tc)
#define W5 W5get(Tc)

constexpr ftype w_get(const int i, const ftype Tc) {
  const ftype arr[QN] = {W0, W1, W2};
  return arr[i];
}
#define TO_SEQ_ELEM_WGET(z, n, text) w_get(n,TLat),
constexpr ftype _w[QN] = {
  BOOST_PP_REPEAT( QN, TO_SEQ_ELEM_WGET, w_get )
};

