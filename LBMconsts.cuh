#pragma once
#include "params.h"
#include <boost/preprocessor.hpp>
namespace LBMconsts {

enum {Void=1,Fluid=2,Iface=3,Solid=4,Box=5, I_fluid=6, I_void=7, notIface=8};
//enum {Void=1,Fluid=2,Iface=3,I_fluid=4,I_void=5,Solid=6,Box=7};
//static const unsigned IfaceAll = (1<<Iface)|(1<<I_fluid)|(1<<I_void); 
namespace Change{
    enum {I2F=1,I2V};
};

constexpr const ftype dx=1.0;
constexpr const ftype dy=1.0;
constexpr const ftype dz=1.0;
constexpr const ftype dt=1.0;

//#define D3Q125_POD
#define D3Q39_POD
//#define D3Q27_EXP_POD
//#define D3Q27_POD
//#define D3Q19_POD
//#define D3Q15_POD
//#define D3Q13_POD
//#define D2Q36_POD
//#define D2Q25_POD
//#define D2Q17_POD
//#define D2Q16_POD
//#define D2Q12_POD
//#define D2Q9_POD
//#define D2Q7_POD
//#define D2Q6_POD
//#define D1Q6_POD
//#define D1Q5_POD
//#define D1Q3_POD
//#define D1Q1_POD

//#define D3Q125
//#define D3Q64
//#define D3Q27

//#define D2Q81
//#define D2Q49
//#define D2Q25
//#define D2Q9

//#define D3Q19
//#define D3Q13
// #define D1Q9
// #define D1Q7
// #define D1Q5
// #define D1Q3
// #define D2Q9
// #define D2Q5

template<int B, int ...Btail> constexpr void debug_consts_assert() { static_assert(B, ""); }

#ifdef __CUDA_ARCH__
#define HOST_DEV_CHOOSE __device__
#else
#define HOST_DEV_CHOOSE constexpr
#endif

#define TO_SEQ_ELEM(z, n, text) text[n],

#define TO_SEQ_ELEM_CHECK(z, n, text) int(text.checkval[n]),
#define SEQ_LIST_COLLECTION_CHECK(arg, n) BOOST_PP_REPEAT(n, TO_SEQ_ELEM_CHECK, arg)
   #define SEQ_LIST_COLLECTION(arg, n) BOOST_PP_REPEAT(n, TO_SEQ_ELEM, arg)

#ifdef D1Q3
  #include "D1Q3_USUAL_consts.cuh"
#elif defined D1Q5
  #include "D1Q5_USUAL_consts.cuh"
#elif defined D1Q7
  #include "D1Q7_USUAL_consts.cuh"
#elif defined D1Q9
  #include "D1Q9_USUAL_consts.cuh"

#elif defined D3Q39_POD
  #include "D3Q39_POD_consts.cuh"
#elif defined D3Q27_EXP_POD
  #include "D3Q27_EXP_POD_consts.cuh"
#elif defined D3Q13_POD
  #include "D3Q13_POD_consts.cuh"
#elif defined D3Q125_POD
  #include "D3Q125_POD_consts.cuh"
#elif defined D3Q27_POD
  #include "D3Q27_POD_consts.cuh"
#elif defined D3Q19_POD
  #include "D3Q19_POD_consts.cuh"
#elif defined D3Q15_POD
  #include "D3Q15_POD_consts.cuh"
#elif defined D2Q36_POD
  #include "D2Q36_POD_consts.cuh"
#elif defined D2Q25_POD
  #include "D2Q25_POD_consts.cuh"
#elif defined D2Q17_POD
  #include "D2Q17_POD_consts.cuh"
#elif defined D2Q16_POD
  #include "D2Q16_POD_consts.cuh"
#elif defined D2Q12_POD
  #include "D2Q12_POD_consts.cuh"
#elif defined D2Q9_POD
  #include "D2Q9_POD_consts.cuh"
#elif defined D2Q7_POD
  #include "D2Q7_POD_consts.cuh"
#elif defined D2Q6_POD
  #include "D2Q6_POD_consts.cuh"
#elif defined D1Q5_POD
  #include "D1Q5_POD_consts.cuh"
#elif defined D1Q6_POD
  #include "D1Q6_POD_consts.cuh"
#elif defined D1Q3_POD
  #include "D1Q3_POD_consts.cuh"
#elif defined D1Q1_POD
  #include "D1Q1_POD_consts.cuh"

#elif defined D3Q27
  #include "D3Q27_USUAL_consts.cuh"
#elif defined D3Q64
  #include "D3Q64_USUAL_consts.cuh"
#elif defined D3Q125
  #include "D3Q125_USUAL_consts.cuh"

#elif defined D2Q25
  #include "D2Q25_USUAL_consts.cuh"
#elif defined D2Q49
  #include "D2Q49_USUAL_consts.cuh"
#elif defined D2Q81
  #include "D2Q81_USUAL_consts.cuh"

#elif defined D3Q19
  #include "D3Q19_USUAL_consts.cuh"
#elif defined D3Q15
  #include "D3Q15_USUAL_consts.cuh"
#elif defined D2Q5
  #include "D2Q5_USUAL_consts.cuh"
#elif defined D2Q9
  #include "D2Q9_USUAL_consts.cuh"
#endif

const int Qn=QN;

template<class T, const int Narr=Qn> struct Consts_Wrap {
  constexpr Consts_Wrap():arr() {}
  constexpr T operator[](int iq) const { return arr[iq]; }
  T arr[Narr];
};

#ifndef SEMI_LAGRANGIAN 
HOST_DEV_CHOOSE const int3 e[Qn] = { SEQ_LIST_COLLECTION(_e,QN) };

#ifndef WEIGHTS_MANUALLY
constexpr const struct _W: public Consts_Wrap<ftype> {
  constexpr _W():Consts_Wrap<ftype>() {
    for(int i=0; i<Qn; i++) {
      const auto sum = abs(_e[i].x)+abs(_e[i].y)+abs(_e[i].z);
      const ftype wchoices[] = {W0, W1, W2, W3};
      arr[i] = wchoices[sum];
    }
  }
} _w;
#endif
template<int DIRX,int DIRY,int DIRZ> struct _Reverse: public Consts_Wrap<char> {
  constexpr _Reverse():Consts_Wrap<char>() {
    for(int i=0; i<Qn; i++) for(int j=0; j<Qn; j++) {
      const bool matchX = _e[i].x==_e[j].x && DIRX==0 || _e[i].x==-_e[j].x && DIRX;
      const bool matchY = _e[i].y==_e[j].y && DIRY==0 || _e[i].y==-_e[j].y && DIRY;
      const bool matchZ = _e[i].z==_e[j].z && DIRZ==0 || _e[i].z==-_e[j].z && DIRZ;
      if(matchX && matchY && matchZ) { arr[i]=j; break; }
    }
  }
};
#else
template<int DIRX,int DIRY,int DIRZ> struct _Reverse: public Consts_Wrap<char> {
  constexpr _Reverse():Consts_Wrap<char>() {
    for(int i=0; i<Qn; i++) for(int j=0; j<Qn; j++) {
      const bool matchX = _ef[i].x==_ef[j].x && DIRX==0 || _ef[i].x==-_ef[j].x && DIRX;
      const bool matchY = _ef[i].y==_ef[j].y && DIRY==0 || _ef[i].y==-_ef[j].y && DIRY;
      const bool matchZ = _ef[i].z==_ef[j].z && DIRZ==0 || _ef[i].z==-_ef[j].z && DIRZ;
      if(matchX && matchY && matchZ) { arr[i]=j; break; }
    }
  }
};
#endif

constexpr const _Reverse<1,0,0> _reverse100;
constexpr const _Reverse<0,1,0> _reverse010;
constexpr const _Reverse<0,0,1> _reverse001;
constexpr const _Reverse<1,1,1> _reverse111;

//TODO remake everything to std::array
//#include "utility"
//auto ints = std::make_index_sequence<int, Qn>{};

HOST_DEV_CHOOSE const ftype w[Qn] = { SEQ_LIST_COLLECTION(_w, QN) } ;
HOST_DEV_CHOOSE const int reverseX[Qn]   = { SEQ_LIST_COLLECTION(_reverse100, QN) } ;
HOST_DEV_CHOOSE const int reverseY[Qn]   = { SEQ_LIST_COLLECTION(_reverse010, QN) } ;
HOST_DEV_CHOOSE const int reverseZ[Qn]   = { SEQ_LIST_COLLECTION(_reverse001, QN) } ;
HOST_DEV_CHOOSE const int reverseXYZ[Qn] = { SEQ_LIST_COLLECTION(_reverse111, QN) } ;

//const long int No=(Qn+1)/2;
const long int No=1;

constexpr const ftype dcs2 = 1./cs2;
const ftype dcs4 = 1./(cs2*cs2);
const ftype dcs6 = 1./(cs2*cs2*cs2);

static_assert(!isnan(dcs2));
static_assert(cs2!=0);

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

// #ifdef __CUDA_ARCH__
// #define e_c e_const
// #define w_c w_const
// #else
// #define e_c e_host
// #define w_c w_host
// #endif
};
