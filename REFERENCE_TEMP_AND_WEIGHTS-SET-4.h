//---- See the publication------// 
//     Chikatamarla, S. & Karlin, Iliya. (2009).
//     Lattices for the lattice Boltzmann method. 
//     Physical review. E, Statistical, nonlinear, and soft matter physics. 79. 046701.
//       10.1103/PhysRevE.79.046701.
//-----------------------------//

constexpr const ftype _cs2 = ( sqrtNR(double(ec1d*ec1d+ec2d*ec2d-10*ec1d*ec2d)) + ec1d + ec2d )/6.;
const ftype cs2 = _cs2;

const ftype TLat=cs2;

const ftype W1 = double( ec1d-5*ec2d + sqrtNR(double(ec1d*ec1d+ec2d*ec2d-10*ec1d*ec2d)) )/(12*(ec1d-ec2d));
const ftype W2 = double( 5*ec1d-ec2d - sqrtNR(double(ec1d*ec1d+ec2d*ec2d-10*ec1d*ec2d)) )/(12*(ec1d-ec2d));

