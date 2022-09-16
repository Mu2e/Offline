#include "Offline/Mu2eKinKal/inc/KKConstantBField.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  using Grad = ROOT::Math::SMatrix<double,3>;
  using SVEC3 = KinKal::SVEC3;

  void KKConstantBField::print(std::ostream& os) const {
    os << "KKConstantBField , value " << bfield_ << std::endl;
  }

  Grad KKConstantBField::fieldGrad(VEC3 const& position) const { return Grad(); }
}
