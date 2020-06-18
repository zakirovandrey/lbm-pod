#define DIM 3
#define QN 13
#define QN_IN_DIM 104   // QN*2^DIM
constexpr const int3 _e[QN] = {
 { 0, 0, 0},
 { -1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 1, 0},
 { -1, 0,-1}, { 1, 0,-1}, {-1, 0, 1}, { 1, 0, 1}, 
 {  0,-1,-1}, { 0, 1,-1}, { 0,-1, 1}, { 0, 1, 1}
};
const ftype W0=1./3;
const ftype W1=1./18;
const ftype W2=1./36;
const ftype W3=0;
const ftype cs2 = dx*dx/3.;

//#define NON_ISOTHERMAL_RELAXATION

#define E_MATRIX_X(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f1 +f2-f3 +f5-f6   +f8-f9+f11-f12+f14
#define E_MATRIX_Y(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f3 -f4-f5+f6+f7+f8 -f15+f16-f17+f18
#define E_MATRIX_Z(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f9-f10-f11+f12+f13+f14-f15-f16+f17+f18


