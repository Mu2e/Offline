//
// Create KinKalGeom objects. These depend on other objects served by GeometryService so this must be external to KinKalGeom itself
// Original author: Dave Brown (LBNL) 4/2026
//
#include "Offline/GeometryService/inc/KinKalGeomMaker.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/KinKalGeom/inc/StoppingTarget.hh"
#include "Offline/KinKalGeom/inc/TestCRV.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

namespace mu2e {
  using KinKal::VEC3;
  using KinKal::Cylinder;
  using KinKal::Disk;
  using KinKal::Surface;
  using KinKal::Frustrum;
  using KinKal::Annulus;
  using KinKal::Rectangle;
  using CylPtr = std::shared_ptr<KinKal::Cylinder>;
  using DiskPtr = std::shared_ptr<KinKal::Disk>;
  using RecPtr = std::shared_ptr<KinKal::Rectangle>;
  using AnnPtr = std::shared_ptr<KinKal::Annulus>;
  using FruPtr = std::shared_ptr<KinKal::Frustrum>;
  using SurfacePtr = std::shared_ptr<KinKal::Surface>;
  using KKGMap = std::multimap<SurfaceId,SurfacePtr>;

  std::unique_ptr<KinKalGeom>& KinKalGeomMaker::makeKKG() {
    kkg_ = std::make_unique<KinKalGeom>();
    makeTracker();
    makeDS();
    makeTarget();
    makeTCRV();
    return kkg_;
  }

  void KinKalGeomMaker::makeTracker() {
      // surfaces need to match with virtual detectors. The following is extracted from VirtualDetectorMaker and needs to be updated if that changes.
    // Note that these are placed at the center of the VDs, which have half-thickness of 0.01mm. Since the VD hits are recorded where the SimParticle
    // enters the volume, the reco track will be sampled at a different position depending on the track direction by that amount. This is a fundamental
    // discrepancy between reco and sim data
    auto const& tracker = *(GeomHandle<mu2e::Tracker>());
    GeomHandle<Beamline> bg;
    GeomHandle<DetectorSystem> det;
    double vdHL(0.01); // hardcoded in VirtualDetectorMaker line 56
    // below are from VirtualDetectorMaker lnes 241-244
    double zFrontGlobal = tracker.g4Tracker()->mother().position().z()-tracker.g4Tracker()->mother().tubsParams().zHalfLength()-vdHL;
    double zBackGlobal  = tracker.g4Tracker()->mother().position().z()+tracker.g4Tracker()->mother().tubsParams().zHalfLength()+vdHL;
    // the 0.4 below comes from offsets in the mother volume nesting.
    double zFrontLocal  = zFrontGlobal - tracker.g4Tracker()->z0() + 0.4;
    double zBackLocal   = zBackGlobal  - tracker.g4Tracker()->z0() - 0.4;
    double zMidLocal = 10.1; // 10.1 is hard-coded in VirtualDetectorMaker line 224
    double halfLen = 0.5*(zBackLocal-zFrontLocal);
    double orvd = tracker.g4Tracker()->mother().tubsParams().outerRadius();
    double irvd = tracker.g4Tracker()->getInnerTrackerEnvelopeParams().innerRadius();

    // std::cout << " zfronloc " << zFrontLocal << " zmidloc " << zMidLocal << " zbackloc " << zBackLocal << " orvd " << orvd << " irvd " << irvd << std::endl;

    // cylinders are defined by TT_outer (_inner) virtual detectors
    // Disks are defined to match TT_front (mid, back) virtual detectors
    auto outer = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),orvd,halfLen);
    auto inner = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),irvd,halfLen);
    // expand the disk radii to meet the DS
    auto front = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zFrontLocal),950.);
    auto mid = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zMidLocal),950.);
    auto back = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zBackLocal),950.);
    // add all these to the map
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Front),std::static_pointer_cast<Surface>(front)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Mid),std::static_pointer_cast<Surface>(mid)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Back),std::static_pointer_cast<Surface>(back)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Inner),std::static_pointer_cast<Surface>(inner)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TT_Outer),std::static_pointer_cast<Surface>(outer)));
    // construct the tracker object
    kkg_->tracker_ = std::make_unique<KKGeom::Tracker>(outer,inner,front,mid,back);
  }

  void KinKalGeomMaker::makeDS() {
    auto inner= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),950,5450);
    auto outer= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),1328,5450); // bounding surfaces
    auto front= std::make_shared<Disk>(outer->frontDisk());
    auto back= std::make_shared<Disk>(outer->backDisk());
    auto ipa= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-2770),300.0,500.0);
    auto ipafront= std::make_shared<Disk>(ipa->frontDisk());
    auto ipaback= std::make_shared<Disk>(ipa->backDisk());
    auto opa= std::make_shared<Frustrum>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-3766),454.0,728.4,2125.0); // inner surface
    auto tsda= std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-5967),235.0,525.0); // back surface

    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Front),std::static_pointer_cast<Surface>(front)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Back),std::static_pointer_cast<Surface>(back)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Inner),std::static_pointer_cast<Surface>(inner)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Outer),std::static_pointer_cast<Surface>(outer)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA),std::static_pointer_cast<Surface>(ipa)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA_Front),std::static_pointer_cast<Surface>(ipafront)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA_Back),std::static_pointer_cast<Surface>(ipaback)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::OPA),std::static_pointer_cast<Surface>(opa)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TSDA),std::static_pointer_cast<Surface>(tsda)));

    kkg_->ds_ = std::make_unique<KKGeom::DetectorSolenoid>(inner, outer, front, back, ipa, ipafront, ipaback, opa, tsda);
  }

  void KinKalGeomMaker::makeTarget() {
    // currently use hard-coded geometry
    auto outer = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-4300),75,400.0);
    auto inner= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-4300),21.5,400.0);
    auto front= std::make_shared<Disk>(outer->frontDisk());
    auto back= std::make_shared<Disk>(outer->backDisk());
    double startz = -4700;
    double endz = -3900;
    double dz = (endz-startz)/36.0;

    std::vector<AnnPtr> foils;
    for(int ifoil=0;ifoil < 37; ++ifoil){
      double zpos = startz + ifoil*dz;
      foils.push_back(std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zpos),21.5,75));
    }

    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Front),std::static_pointer_cast<Surface>(front)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Back),std::static_pointer_cast<Surface>(back)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Inner),std::static_pointer_cast<Surface>(inner)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Outer),std::static_pointer_cast<Surface>(outer)));
    for(size_t ifoil=0;ifoil < foils.size();++ifoil){
      kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Foils,ifoil),std::static_pointer_cast<Surface>(foils[ifoil])));
    }

    kkg_->st_ = std::make_unique<KKGeom::StoppingTarget>(outer,inner,front,back,foils);
  }

  void KinKalGeomMaker::makeTCRV() {
    // currently use hard-coded geometry
    auto ex1= std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(0.0,4775,-438),3000,1675); // layer widths are approximate FIXME
    auto t1= std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,1.0), VEC3(0.0,4625,-438),1185,850);
    auto t2= std::make_shared<Rectangle>(VEC3(0.0,1.0,0.0),VEC3(1.0,0.0,0.0), VEC3(0.0,4925,-438),1600,850);

    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TCRV,0),std::static_pointer_cast<Surface>(t1)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TCRV,1),std::static_pointer_cast<Surface>(ex1)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TCRV,2),std::static_pointer_cast<Surface>(t2)));

    kkg_->tcrv_ = std::make_unique<KKGeom::TestCRV>(ex1,t1,t2);
  }
}
