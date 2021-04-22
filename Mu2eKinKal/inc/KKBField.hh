#ifndef Mu2eKinKal_KKBField_hh
#define Mu2eKinKal_KKBField_hh
//
//  Wrapper to Mu2e BField map for KinKal
//
// Mu2e includes 
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// KinKal includes
#include "KinKal/General/BFieldMap.hh"

namespace mu2e 
{
  using VEC3 = KinKal::VEC3;
  class KKBField : public KinKal::BFieldMap {
    public:
      using Grad = ROOT::Math::SMatrix<double,3>; // field gradient: ie dBi/d(x,y,z)
    // construct from BField object and system translator.  
    // This should be a single BField map valid in the detector system, to avoid making continuous translations FIXME!
      KKBField(BFieldManager const& bfmgr, DetectorSystem const& det) : bfmgr_(bfmgr), det_(det) {}
      virtual ~KKBField() {}
      // KinKal BField interface
      // return value of the field at a poin
      virtual VEC3 fieldVect(VEC3 const& position) const override;
      // return BFieldMap gradient = dB_i/dx_j, at a given point
      virtual Grad fieldGrad(VEC3 const& position) const override;
      // return the BFieldMap derivative at a given point along a given velocity, WRT time
      virtual VEC3 fieldDeriv(VEC3 const& position, VEC3 const& velocity) const override;
    private:
      BFieldManager const& bfmgr_;
      DetectorSystem const& det_;
  };
}
#endif





