//---- See the publication------// 
//     Chikatamarla, S. & Karlin, Iliya. (2009).
//     Lattices for the lattice Boltzmann method. 
//     Physical review. E, Statistical, nonlinear, and soft matter physics. 79. 046701.
//       10.1103/PhysRevE.79.046701.
//-----------------------------//

constexpr const ftype _cs2 = ( -sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) + ec1d + ec2d )/10.;
const ftype cs2 = _cs2;

const ftype TLat=cs2;

const ftype W0 = double( -ec1d*ec1d-ec2d*ec2d + 18*ec1d*ec2d + (ec1d  +  ec2d)*sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) )/(25*ec1d*ec2d);
const ftype W1 = double( 3*ec1d*ec1d-2*ec2d*ec2d-9*ec1d*ec2d - (3*ec1d-2*ec2d)*sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) )/(100*ec1d*(ec1d-ec2d));
const ftype W2 = double( 3*ec2d*ec2d-2*ec1d*ec1d-9*ec1d*ec2d - (3*ec2d-2*ec1d)*sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) )/(100*ec2d*(ec2d-ec1d));

#ifdef NON_GAUSSIAN_EQUILIBRIA

constexpr const ftype _cs2 = ( +sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) + ec1d + ec2d )/10.;
const ftype cs2 = _cs2;

const ftype TLat=cs2;

const ftype W0 = double( -ec1d*ec1d-ec2d*ec2d + 18*ec1d*ec2d - (ec1d  +  ec2d)*sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) )/(25*ec1d*ec2d);
const ftype W1 = double( 3*ec1d*ec1d-2*ec2d*ec2d-9*ec1d*ec2d + (3*ec1d-2*ec2d)*sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) )/(100*ec1d*(ec1d-ec2d));
const ftype W2 = double( 3*ec2d*ec2d-2*ec1d*ec1d-9*ec1d*ec2d + (3*ec2d-2*ec1d)*sqrtNR(double(ec1d*ec1d+ec2d*ec2d-14/3.*ec1d*ec2d)) )/(100*ec2d*(ec2d-ec1d));

#endif

