#ifndef GeometryService_KinKalGeomMaker_hh
#define GeometryService_KinKalGeomMaker_hh
//
// Create KinKalGeom objects. These depend on other objects served by GeometryService so this must be external to KinKalGeom itself
// Original author: Dave Brown (LBNL) 4/2026
//
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
namespace mu2e {
  class KinKalGeomMaker {
    public:
      std::unique_ptr<KinKalGeom>& makeKKG();
    private:
      void makeTracker();
      void makeDS();
      void makeTarget();
      void makeCalo();
      void makeTCRV();
      std::unique_ptr<KinKalGeom> kkg_;
  };
}
#endif
