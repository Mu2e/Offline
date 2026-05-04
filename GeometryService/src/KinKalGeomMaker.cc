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
#include "Offline/KinKalGeom/inc/Calo.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"

namespace mu2e {
  using KinKal::VEC3;
  using KinKal::Cylinder;
  // NOTE: Do NOT use `using KinKal::Disk` because CalorimeterGeom also defines mu2e::Disk
  // Use explicit KinKal::Disk instead
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
    makeCalo();
    return kkg_;
  }

  void KinKalGeomMaker::makeCalo() {
    // construct calorimeter geometry from GeometryService
    auto const& calo_det = *(GeomHandle<mu2e::Calorimeter>());
    auto const& tracker = *(GeomHandle<mu2e::Tracker>());

    // Extract geometry from the first two disks (D0 and D1)
    auto const& disk0 = calo_det.disk(0);
    auto const& disk1 = calo_det.disk(1);
    auto const& geom0 = disk0.geomInfo();
    auto const& geom1 = disk1.geomInfo();

    // Extract Z positions and radii from geometry (in global coordinates)
    double z0_global = geom0.origin().z();
    double z1_global = geom1.origin().z();
    double r0_inner = geom0.innerEnvelopeR();
    double r0_outer = geom0.outerEnvelopeR();
    double r1_inner = geom1.innerEnvelopeR();
    double r1_outer = geom1.outerEnvelopeR();

    // Use true disk geometry: compute symmetric front/back from disk center
    double diskHalfLen_0 = 0.5 * geom0.size().z();
    double diskHalfLen_1 = 0.5 * geom1.size().z();

    double z0_front_global = z0_global - diskHalfLen_0;
    double z0_back_global  = z0_global + diskHalfLen_0;
    double z1_front_global = z1_global - diskHalfLen_1;
    double z1_back_global  = z1_global + diskHalfLen_1;

    // Convert from global to local coordinates (relative to tracker center)
    double tracker_z0 = tracker.g4Tracker()->z0();
    double z0 = z0_global - tracker_z0;
    double z1 = z1_global - tracker_z0;
    double z0_front = z0_front_global - tracker_z0;
    double z0_back  = z0_back_global - tracker_z0;
    double z1_front = z1_front_global - tracker_z0;
    double z1_back  = z1_back_global - tracker_z0;

    // Construct Calo with geometry in local coordinates
    kkg_->calo_ = std::make_unique<KKGeom::Calo>(z0, z1, r0_inner, r0_outer, r1_inner, r1_outer, z0_front, z0_back, z1_front, z1_back);
    auto const& calo = *kkg_->calo_;

    // Add all calo surfaces to the map
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_Outer),std::static_pointer_cast<Surface>(calo.EMC_Disk_0_OuterPtr())));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_Inner),std::static_pointer_cast<Surface>(calo.EMC_Disk_0_InnerPtr())));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_Inner),std::static_pointer_cast<Surface>(calo.EMC_Disk_1_InnerPtr())));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_Outer),std::static_pointer_cast<Surface>(calo.EMC_Disk_1_OuterPtr())));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_Front),std::static_pointer_cast<Surface>(calo.EMC_Disk_0_FrontPtr())));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_Front),std::static_pointer_cast<Surface>(calo.EMC_Disk_1_FrontPtr())));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_0_Back),std::static_pointer_cast<Surface>(calo.EMC_Disk_0_BackPtr())));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::EMC_Disk_1_Back),std::static_pointer_cast<Surface>(calo.EMC_Disk_1_BackPtr())));
  }

  void KinKalGeomMaker::makeTracker() {
      // surfaces need to match with virtual detectors. The following is extracted from VirtualDetectorMaker and needs to be updated if that changes.
    // Note that these are placed at the center of the VDs, which have half-thickness of 0.01mm. Since the VD hits are recorded where the SimParticle
    // enters the volume, the reco track will be sampled at a different position depending on the track direction by that amount. This is a fundamental
    // discrepancy between reco and sim data
    auto const& tracker = *(GeomHandle<mu2e::Tracker>());
    GeomHandle<Beamline> bg;
    GeomHandle<DetectorSystem> det;
    GeomHandle<DetectorSolenoid> ds;
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
    double irds = ds->rIn1();
    // cylinders are defined by TT_outer (_inner) virtual detectors
    // Disks are defined to match TT_front (mid, back) virtual detectors
    auto outer = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),orvd,halfLen);
    auto inner = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),irvd,halfLen);
    // expand the disk radii to the DS
    auto front = std::make_shared<KinKal::Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zFrontLocal),irds);
    auto mid = std::make_shared<KinKal::Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zMidLocal),irds);
    auto back = std::make_shared<KinKal::Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zBackLocal),irds);
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
    // currently use hard-coded geometry
    auto inner= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),950,5450);
    auto outer= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),1328,5450); // bounding surfaces
    auto front= std::make_shared<KinKal::Disk>(outer->frontDisk());
    auto back= std::make_shared<KinKal::Disk>(outer->backDisk());
    auto ipa= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-2770),300.0,500.0);
    auto ipafront= std::make_shared<KinKal::Disk>(ipa->frontDisk());
    auto ipaback= std::make_shared<KinKal::Disk>(ipa->backDisk());
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
    auto front= std::make_shared<KinKal::Disk>(outer->frontDisk());
    auto back= std::make_shared<KinKal::Disk>(outer->backDisk());
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
