#ifndef GeometryService_KinKalGeomMaker_hh
#define GeometryService_KinKalGeomMaker_hh
//
// Create KinKalGeom objects. These depend on other objects served by GeometryService so this must be external to KinKalGeom itself
// Original author: Dave Brown (LBNL) 4/2026
//
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <string>
namespace mu2e {
  class DetectorSystem;
  class KinKalGeomMaker {
    public:
      explicit KinKalGeomMaker(int debug) : debug_(debug) {}
      std::unique_ptr<KinKalGeom>& makeKKG();
    private:
      void makeTracker();
      void makeDS();
      void makeTarget();
      void makeCRV();
      void makePassiveMaterials();
      // append one averaged rectangular concrete passive-material plane (see .cc)
      void addConcretePlane(DetectorSystem const& det, int normalAxis,
          CLHEP::Hep3Vector const& centerMu2e, double hw1, double hw2,
          double halfThickness, std::string const& material);
      // run-1 building hatch concrete when ExtShieldDownstream concrete is zeroed
      void addBuildingHatchConcrete(DetectorSystem const& det);
      std::unique_ptr<KinKalGeom> kkg_;
      int debug_ = 0;
  };
}
#endif
