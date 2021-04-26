#include "Mu2eKinKal/inc/KKBField.hh"
namespace mu2e {
  using Grad = ROOT::Math::SMatrix<double,3>;
  using SVEC3 = KinKal::SVEC3;

  VEC3 KKBField::fieldVect(VEC3 const& position) const {
    // change coordinates to mu2e; the map should be native in detector coordinates FIXME!
    CLHEP::Hep3Vector vpoint(position.x(),position.y(),position.z());
    CLHEP::Hep3Vector vpoint_mu2e = det_.toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr_.getBField(vpoint_mu2e);
    return VEC3(field.x(),field.y(),field.z());
  }
      
  Grad KKBField::fieldGrad(VEC3 const& position) const {
    Grad retval;
    auto dBdx = fieldDeriv(position,VEC3(1.0,0.0,0.0));
    auto dBdy = fieldDeriv(position,VEC3(0.0,1.0,0.0));
    auto dBdz = fieldDeriv(position,VEC3(0.0,0.0,1.0));
    SVEC3 dBdxv(dBdx.X(),dBdx.Y(), dBdx.Z());
    SVEC3 dBdyv(dBdy.X(),dBdy.Y(), dBdy.Z());
    SVEC3 dBdzv(dBdz.X(),dBdz.Y(), dBdz.Z());
    retval.Place_in_row(dBdxv,0,0);
    retval.Place_in_row(dBdyv,1,0);
    retval.Place_in_row(dBdzv,2,0);
    return retval;
  }
// numerical derivatives for now: TODO!
  VEC3 KKBField::fieldDeriv(VEC3 const& position, VEC3 const& velocity) const {
    static double dt(0.01); // 10 psec
    VEC3 start = fieldVect(position);
    VEC3 end = fieldVect(position + velocity*dt);
    return (end-start)/dt;
  }

}
