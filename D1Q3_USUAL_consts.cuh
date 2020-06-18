#define DIM 1
#define QN 3
#define QN_IN_DIM 6   // QN*2^DIM
constexpr const int3 _e[QN] = { { 0, 0, 0}, { 1, 0, 0}, { -1, 0, 0} };
constexpr const ftype _cs2 = dx*dx/3.;
const ftype cs2 = _cs2;
const ftype TLat=cs2;

const ftype W0=2./3;
const ftype W1=1./6;
const ftype W2=0;
const ftype W3=0;


//#define NON_ISOTHERMAL_RELAXATION

#define WEIGHTS_MANUALLY
constexpr ftype _w[QN] = {
  W0, W1, W1
};
