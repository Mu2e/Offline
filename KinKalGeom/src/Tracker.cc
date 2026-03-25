#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
namespace mu2e {
  namespace KinKalGeom {
    using KinKal::VEC3;
    using KinKal::Cylinder;
    using KinKal::Disk;
    // currently use hard-coded geometry.  Note: these are only comparable to the MC virtual detector positions
    // at the ~10 um level, as the G4 virtual detectors are that thick, and the steps are semi-random (step tolerance)
    Tracker::Tracker() {}

    void Tracker::check_init() const {
      if(!initialized_){
        // cast off const
        (const_cast<Tracker*>(this))->initialize();
        initialized_ = true;
      }
    }
    void Tracker::initialize() {
      // surfaces need to match with virtual detectors. The following is extracted from VirtualDetectorMaker and needs to be updated if that changes.
      double vdHL(0.01); // from geom.txt
      auto tracker = *(GeomHandle<mu2e::Tracker>());
      GeomHandle<Beamline> bg;
      GeomHandle<DetectorSystem> det;
      double solenoidOffset = bg->solenoidOffset();
      CLHEP::Hep3Vector ttOffset(-solenoidOffset,0.,tracker.g4Tracker()->z0());
      double zFrontGlobal = tracker.g4Tracker()->mother().position().z()-tracker.g4Tracker()->mother().tubsParams().zHalfLength()-vdHL;
      double zBackGlobal  = tracker.g4Tracker()->mother().position().z()+tracker.g4Tracker()->mother().tubsParams().zHalfLength()+vdHL;
      double zFrontLocal  = zFrontGlobal - tracker.g4Tracker()->z0();
      double zBackLocal   = zBackGlobal  - tracker.g4Tracker()->z0();
      double zMidLocal = 0.5*(zFrontLocal + zBackLocal);
      double halfLen = 0.5*(zBackLocal-zFrontLocal);
      CLHEP::Hep3Vector ttOffset_det = det->toDetector(ttOffset);

//      auto ttorigin = det->getOrigin();
      std::cout << "solff " << solenoidOffset <<" zfronloc " << zFrontLocal << " zbackloc " << zBackLocal << " ttoff " << ttOffset_det << std::endl;

      // cylinders are defined by TT_outer (_inner) virtual detectors
      // Disks are defined to match TT_front (mid, back) virtual detectors
      outer_  = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),850.11,halfLen);
      inner_ = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),376.9,halfLen);
      // expand the disk radii to meet the DS
      front_ = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zFrontLocal),950.);
      mid_ = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,10.09),950.);
      back_ = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zBackLocal),950.);
    }
  }
}
