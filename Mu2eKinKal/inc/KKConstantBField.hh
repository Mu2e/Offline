#ifndef Mu2eKinKal_KKConstantBField_hh
#define Mu2eKinKal_KKConstantBField_hh
//
//  Wrapper to Mu2e BField map for KinKal
//
// Mu2e includes
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// KinKal includes
#include "KinKal/General/BFieldMap.hh"

namespace mu2e
{
  using VEC3 = KinKal::VEC3;
  class KKConstantBField : public KinKal::BFieldMap {
    public:
      using Grad = ROOT::Math::SMatrix<double,3>; // field gradient: ie dBi/d(x,y,z)
      KKConstantBField(VEC3 const& bf) : bfield_(bf) {}
      virtual ~KKConstantBField() {}
      VEC3 fieldVect(VEC3 const& position) const override { return bfield_; }
      // return BFieldMap gradient = dB_i/dx_j, at a given point
      Grad fieldGrad(VEC3 const& position) const override;
      // return the BFieldMap derivative at a given point along a given velocity, WRT time
      VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override { return VEC3(); }
      // is the point inside the region defined by this map?
      bool inRange(VEC3 const& position) const override { return true; }
      void print(std::ostream& os ) const override;
    private:
      VEC3 bfield_;
  };
}
#endif
