#define DIM 2
#define QN 5
#define QN_IN_DIM 20   // QN*2^DIM
constexpr const int3 _e[QN] = {
  /* 0 */ { 0, 0, 0}, 
  /* 1 */ { 1, 0, 0}, 
  /* 2 */ { 0, 1, 0}, 
  /* 3 */ {-1, 0, 0}, 
  /* 4 */ { 0,-1, 0}, 
};

constexpr const ftype cs2 = 1./3.;
const ftype TLat=cs2;
constexpr ftype W0get(const ftype Tv) { return 1./2.; }
constexpr ftype W1get(const ftype Tv) { return 1./8.; }
constexpr ftype W2get(const ftype Tv) { return 0; }

const ftype Tmin = 0;
const ftype Tmax = 100;


//#define NON_ISOTHERMAL_RELAXATION

#define WEIGHTS_MANUALLY

#define W0 W0get(Tc)
#define W1 W1get(Tc)
#define W2 W2get(Tc)
#define W3 W3get(Tc)
#define W4 W4get(Tc)
#define W5 W5get(Tc)

constexpr ftype w_get(const int i, const ftype Tc) {
  const ftype arr[QN] = {
  W0, W1,W1,W1,W1
  };
  return arr[i];
}
#define TO_SEQ_ELEM_WGET(z, n, text) w_get(n,TLat),
constexpr ftype _w[QN] = {
  BOOST_PP_REPEAT( QN, TO_SEQ_ELEM_WGET, w_get )
};

