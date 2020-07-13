#define DIM 3
#define QN 19
#define QN_IN_DIM 152   // QN*2^DIM
constexpr const int3 _e[QN] = {
 { 0 , 0, 0},
 { 1 , 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, {0, 0, 1}, {0, 0, -1},
 { -1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 1, 0},
 { -1, 0,-1}, { 1, 0,-1}, {-1, 0, 1}, { 1, 0, 1}, 
 {  0,-1,-1}, { 0, 1,-1}, { 0,-1, 1}, { 0, 1, 1}
};
constexpr const ftype cs2 = 1./3.;
const ftype TLat=cs2;
constexpr ftype W0get(const ftype Tv) { return 1./3.; }
constexpr ftype W1get(const ftype Tv) { return 1./18.; }
constexpr ftype W2get(const ftype Tv) { return 1./36.; }

const ftype Tmin = 0;
const ftype Tmax = 1;

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
  W0,
  W1,W1,W1,W1,W1,W1,
  W2,W2,W2,W2,W2,W2, W2,W2,W2,W2,W2,W2
  };
  return arr[i];
}
#define TO_SEQ_ELEM_WGET(z, n, text) w_get(n,TLat),
constexpr ftype _w[QN] = {
  BOOST_PP_REPEAT( QN, TO_SEQ_ELEM_WGET, w_get )
};
#define E_MATRIX_X(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f1 +f2-f3 +f5-f6   +f8-f9+f11-f12+f14
#define E_MATRIX_Y(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f3 -f4-f5+f6+f7+f8 -f15+f16-f17+f18
#define E_MATRIX_Z(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f9-f10-f11+f12+f13+f14-f15-f16+f17+f18


