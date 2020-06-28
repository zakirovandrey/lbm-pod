//---- See the publication------// 
//     Chikatamarla, S. & Karlin, Iliya. (2009).
//     Lattices for the lattice Boltzmann method. 
//     Physical review. E, Statistical, nonlinear, and soft matter physics. 79. 046701.
//       10.1103/PhysRevE.79.046701.
//-----------------------------//

static_assert(ec1d==1);
static_assert(ec2d==2*2);
static_assert(ec3d==3*3);

namespace CK{
  constexpr const double p=ec4;
  static_assert(p>4);
  constexpr const double p2=ec4d;
  constexpr const double p4=p2*p2; constexpr const double p6=p2*p2*p2; constexpr const double p8=p4*p4;
  constexpr const double p10=p2*p2*p2*p2*p2; constexpr const double p12=p6*p6;
  constexpr const double D2 = 245*p12 - 10780*p10 + 137004*p8 - 164696*p6 - 3442481*p4 - 14677656*p2 + 3523824;
  constexpr const double D1 = sprout::math::cbrt(35*p6-231*p4-7203*p2-3*sprout::math::sqrt(6.0)*sprout::math::sqrt(D2)+21979);
  constexpr double get_Tpredefined() {
    if(p==5) return 0.75608085259426858; else
    if(p==6) return 0.73202304233490170; else
    if(p==7) return 0.72116103412433775; else
    if(p==8) return 0.71498624076925002;
    else return 0;
  }
  constexpr const ftype _cs2 = get_Tpredefined();
    /*
    p2/36+7./8.
    -1./(36*sprout::math::sqrt(35.0))*sprout::math::sqrt(
        35*(p2+14)*(p2+14)-840*(2*p2+7)+ 12*sprout::math::cbrt(35.0)*D1+ (sprout::math::pow(1235,2./3.)*(-7*p4+110*p2+203)) *D1
    )
    +1./(18*sprout::math::sqrt(70.))* sprout::math::sqrt( 35*(p2+14)*(p2+14)-840*(2*p2+7)-6*sprout::math::cbrt(35.0)*D1+ sprout::math::pow(635,2./3.)*(7*p4-110*p2-203)/D1
    +sprout::math::sqrt(35.0)*(-35*p6+1050*p4-8232*p2+4112)/sprout::math::sqrt(35*(p2+14)*(p2+14)-840*(2*p2+7)+12*sprout::math::cbrt(35.0)*D1+1./D1*(sprout::math::pow(1235,2./3.)*(-7*p4+110*p2+203))) )
  ;*/
}
constexpr const ftype cs2 = CK::_cs2;
const ftype TLat=cs2;

constexpr ftype W0get(const ftype Tv) {
  using namespace CK;
  return ( (Tv*(3*(14-5*Tv)*Tv-49)+36)*p2+3*Tv*(7*Tv*(5*(Tv-2)*Tv+7)-12) ) / (36*p2);
}
constexpr ftype W1get(const ftype Tv) { 
  using namespace CK;
  const ftype Tv2=Tv*Tv;
  const ftype Tv3=Tv*Tv*Tv;
  return Tv*((Tv*(5*Tv-13)+12)*p2+Tv*(5*(13-7*Tv)*Tv-36)) / ( 16*(p2-1) );
}
constexpr ftype W2get(const ftype Tv) {
  using namespace CK;
  const ftype Tv2=Tv*Tv;
  const ftype Tv3=Tv*Tv*Tv;
  return Tv*((-5*(Tv-2)*Tv-3)*p2+Tv*(5*Tv*(7*Tv-10)+9)) / ( 40*(p2-4) );
}
constexpr ftype W3get(const ftype Tv) {
  using namespace CK;
  const ftype Tv2=Tv*Tv;
  const ftype Tv3=Tv*Tv*Tv;
  return Tv*(p2*(15*(Tv-1)*Tv+4)-3*Tv*(5*Tv*(7*Tv-5)+4)) / ( 720*(p2-9) );
}
constexpr ftype W4get(const ftype Tv) {
  using namespace CK;
  const ftype Tv2=Tv*Tv;
  const ftype Tv3=Tv*Tv*Tv;
  return 3*Tv*(7*Tv*(5*(Tv-2)*Tv+7)-12) / ( 2*p2*(p2*(p2-7)*(p2-7)-36)  );
}
