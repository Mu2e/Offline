#include "Offline/Mu2eKinKal/inc/KKBField.hh"
#include "Offline/BFieldGeom/inc/BFieldConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "cetlib_except/exception.h"
namespace mu2e {
  using Grad = ROOT::Math::SMatrix<double,3>;

  KKBField::KKBField(BFieldManager const& bfmgr, DetectorSystem const& det) :
    KinKal::BFieldMap(GeomHandle<BFieldConfig>()->minUsableField()), bfmgr_(bfmgr), det_(det) {}

  VEC3 KKBField::fieldVect(VEC3 const& position) const {
    // change coordinates to mu2e; the map should be native in detector coordinates FIXME!
    CLHEP::Hep3Vector vpoint(position.x(),position.y(),position.z());
    CLHEP::Hep3Vector vpoint_mu2e = det_.toMu2e(vpoint);
    CLHEP::Hep3Vector field;
    //    = bfmgr_.getBField(vpoint_mu2e);
    if(bfmgr_.getBFieldWithStatus(vpoint_mu2e,field))
      return VEC3(field);
// see if there's no maps; that says this is the no-field case
// FIXME need to deal with case when outside all maps
    static const VEC3 nullfield(0.0,0.0,0.0);
    return nullfield;
  }

  Grad KKBField::fieldGrad(VEC3 const& position) const {
    // central differences of the field map: Grad(i,j) = dB_i/dx_j.  The step is small compared to the map grid spacing (25 mm)
    static const double step(5.0); // mm
    Grad retval;
    CLHEP::Hep3Vector vpoint_mu2e = det_.toMu2e(CLHEP::Hep3Vector(position.x(),position.y(),position.z()));
    for(int jdim=0;jdim<3;jdim++){
      CLHEP::Hep3Vector dx; dx[jdim] = step;
      CLHEP::Hep3Vector bplus, bminus;
      // outside the map, or too close to its edge, return a null gradient
      if(!(bfmgr_.getBFieldWithStatus(vpoint_mu2e+dx,bplus) && bfmgr_.getBFieldWithStatus(vpoint_mu2e-dx,bminus))) return Grad();
      auto dbdx = (bplus-bminus)/(2.0*step);
      for(int idim=0;idim<3;idim++) retval(idim,jdim) = dbdx[idim];
    }
    return retval;
  }
  // numerical derivatives for now: TODO!
  VEC3 KKBField::fieldDeriv(VEC3 const& position, VEC3 const& velocity) const {
    static double dt(0.1); // 100 psec, ~3cm.  this is arbitrary and should be set according to the field sampling TODO
    if(!inRange(position))throw std::runtime_error("Position out of range");
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
