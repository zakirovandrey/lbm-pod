//---- See the publication------// 
//     Chikatamarla, S. & Karlin, Iliya. (2009).
//     Lattices for the lattice Boltzmann method. 
//     Physical review. E, Statistical, nonlinear, and soft matter physics. 79. 046701.
//       10.1103/PhysRevE.79.046701.
//-----------------------------//

static_assert(ec1d==1);

namespace CK{
  constexpr const double m=ec2;  constexpr const double n=ec3;
  constexpr const double m2=ec2d;  constexpr const double n2=ec3d;
  constexpr const double m4=m2*m2;    constexpr const double n4=n2*n2;
  constexpr const double m6=m2*m2*m2; constexpr const double n6=n2*n2*n2;
  constexpr const ftype _cs2 = 
    1./210.*(10*(m2+n2+1) - cbrtCR(4.0)*sprout::math::cbrt(-250*m6+825*(n2+1)*m4+75*(11*n4-104*n2+11)*m2-25*(n2+1)*(10*n4-43*n2+10)
      +1./27.*sprout::math::sqrt(4*sprout::math::pow(945*((n2+1)*m2+n2)-225*(m2+n2+1)*(m2+n2+1),3)+455625*sprout::math::pow(10*m6-33*(n2+1)*m4+(-33*n4+312*n2-33)*m2
      +(n2+1)*(10*n4-43*n2+10),2))) - 2*sprout::math::cbrt(50.0)*(5*m4-11*(n2+1)*m2+5*n4-11*n2+5)*sprout::math::pow(-50*m6+165*(n2+1)*m4
      +15*(11*n4-104*n2+11)*m2-5*(n2+1)*(10*n4-43*n2+10) + 1./135.*sprout::math::sqrt(4*sprout::math::pow(945*((n2+1)*m2+n2)-225*(m2+n2+1)*(m2+n2+1),3)
      +455625*sprout::math::pow(10*m6-33*(n2+1)*m4+(-33*n4+312*n2-33)*m2+(n2+1)*(10*n4-43*n2+10),2)),(-1./3.)))
    ;
}
constexpr const ftype cs2 = CK::_cs2;
const ftype TLat=cs2;

constexpr ftype W0get(const ftype Tv) {
  using namespace CK;
  return ( -15*Tv*Tv*Tv + 3*(m2+n2+1)*Tv*Tv - ((n2+1)*m2+n2)*Tv + m2*n2 )/(m2*n2);
}
constexpr ftype W1get(const ftype Tv) { 
  using namespace CK;
  const ftype Tv2=Tv*Tv;
  const ftype Tv3=Tv*Tv*Tv;
  return
    1./(2*(m2-1)*((Tv-1)*m2-3*Tv2+Tv)) * ((Tv*(m2*n2+3*Tv2-(m2+n2)*Tv)*(m2+15*Tv2-3*(m2+1)*Tv)*(-15*Tv3+3*(m2+n2
    +1)*Tv2-((n2+1)*m2+n2)*Tv+m2*n2))*1.0/(n2*(n2-1)*(((Tv-1)*n2-3*Tv2+Tv)*m2+Tv*(n2+15*Tv2-3*(n2+1)*Tv)))
    -(m2-3*Tv)*Tv*(-15*Tv3+3*(m2+n2+1)*Tv2-((n2+1)*m2+n2)*Tv+m2*n2)/n2)
  ;
}
constexpr ftype W2get(const ftype Tv) {
  using namespace CK;
  const ftype Tv2=Tv*Tv;
  const ftype Tv3=Tv*Tv*Tv;
  return
    (Tv*(n2+15*Tv2-3*(n2+1)*Tv)*(-15*Tv3+3*(m2+n2+1)*Tv2-((m2+1)*n2+m2)*Tv+m2*n2))
    *1.0/(2*(n-m)*m2*(m+n)*(m2-1)*(((Tv-1)*m2-3*Tv2+Tv)*n2+Tv*(m2+15*Tv2-3*(m2+1)*Tv))) 
  /*return
    -( Tv*(3*Tv-1)*(-15*Tv3+3*(m2+n2+1)*Tv2-((n2+1)*m2+n2)*Tv+m2*n2) / (m2*n2)
    + (Tv*((Tv-1)*n2-3*Tv2+Tv)*(m2+15*Tv2-3*(m2+1)*Tv)*(-15*Tv3+3*(m2+n2+1)*Tv2-((n2+1)*m2+n2)*Tv+m2*n2))
    * ((m-n)*n2*(m+n)*(((Tv-1)*n2-3*Tv2+Tv)*m2+Tv*(n2+15*Tv2-3*(n2+1)*Tv))))
    * 1.0/(2*(m2-1)*((Tv-1)*m2-3*Tv2+Tv))*/
    ;
}
constexpr ftype W3get(const ftype Tv) {
  using namespace CK;
  const ftype Tv2=Tv*Tv;
  const ftype Tv3=Tv*Tv*Tv;
  return
    (Tv*(m2+15*Tv2-3*(m2+1)*Tv)*(-15*Tv3+3*(m2+n2+1)*Tv2-((n2+1)*m2+n2)*Tv+m2*n2))
    *1.0/(2*(m-n)*n2*(m+n)*(n2-1)*(((Tv-1)*n2-3*Tv2+Tv)*m2+Tv*(n2+15*Tv2-3*(n2+1)*Tv)))
  ;
}

