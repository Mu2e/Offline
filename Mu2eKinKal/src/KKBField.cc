#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  using Grad = ROOT::Math::SMatrix<double,3>;
  using SVEC3 = KinKal::SVEC3;

  VEC3 KKBField::fieldVect(VEC3 const& position) const {
    // change coordinates to mu2e; the map should be native in detector coordinates FIXME!
    CLHEP::Hep3Vector vpoint(position.x(),position.y(),position.z());
    CLHEP::Hep3Vector vpoint_mu2e = det_.toMu2e(vpoint);
    CLHEP::Hep3Vector field;
    //    = bfmgr_.getBField(vpoint_mu2e);
    if(bfmgr_.getBFieldWithStatus(vpoint_mu2e,field))
      return VEC3(field);
    // Outside all field maps the real field is ~0 (detector hall, far from the solenoids), but a
    // CentralHelix evaluated at IDENTICALLY zero field is singular (omega->0 -> 1/omega -> NaN). A
    // bfcorr=true fit's domain walk can transiently sample just outside the maps, so return a tiny
    // NON-zero axial placeholder instead: physically negligible, but it keeps omega finite so the
    // fit can't NaN.
    static const VEC3 nullfield(0.0,0.0,1.0e-3);
    return nullfield;
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
    static double dt(0.1); // 100 psec, ~3cm.  this is arbitrary and should be set according to the field sampling TODO
    // fieldVect is defined everywhere (tiny non-zero placeholder outside the maps), so the derivative
    // can always be evaluated; out of the maps it is ~0. No out-of-range throw: a bfcorr=true fit whose
    // domain walk steps just outside the maps must not be aborted.
    VEC3 start = fieldVect(position);
    VEC3 end = fieldVect(position + velocity*dt);
    return (end-start)/dt;
  }
  bool KKBField::inRange(VEC3 const& position) const {
    // clumsy conversion to CLHEP
    CLHEP::Hep3Vector vpoint(position.x(),position.y(),position.z());
    CLHEP::Hep3Vector vpoint_mu2e = det_.toMu2e(vpoint);
    bool retval (false);
    for(auto const& bfmap : bfmgr_.getInnerMaps()){
      retval |=
        vpoint_mu2e.x() > bfmap->xmin() && vpoint_mu2e.x() < bfmap->xmax() &&
        vpoint_mu2e.y() > bfmap->ymin() && vpoint_mu2e.y() < bfmap->ymax() &&
        vpoint_mu2e.z() > bfmap->zmin() && vpoint_mu2e.z() < bfmap->zmax();
    }
    return retval;
  }

  void KKBField::print(std::ostream& os) const {
    os << "KKBField based on ";
    bfmgr_.print(os);
  }

}
