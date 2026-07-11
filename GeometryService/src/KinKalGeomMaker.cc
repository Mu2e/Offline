//
// Create KinKalGeom objects. These depend on other objects served by GeometryService so this must be external to KinKalGeom itself
// Original author: Dave Brown (LBNL) 4/2026
//
#include "Offline/GeometryService/inc/KinKalGeomMaker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "Offline/ExternalShieldingGeom/inc/ExtShieldDownstream.hh"
#include "Offline/Mu2eHallGeom/inc/Mu2eHall.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/KinKalGeom/inc/StoppingTarget.hh"
#include "Offline/KinKalGeom/inc/CRV.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <map>

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
  using mu2e::KKGeom::KKCRVSector;


  std::unique_ptr<KinKalGeom>& KinKalGeomMaker::makeKKG() {
    kkg_ = std::make_unique<KinKalGeom>();
    makeTracker();
    makeDS();
    makeTarget();
    makeCRV();
    makePassiveMaterials();
    return kkg_;
  }

  // sort by transverse distance
  struct sortCRVSectors {
    bool operator () (KKCRVSector const& sect1, KKCRVSector const& sect2) {
      return sect1.sector_->center().Rho() > sect2.sector_->center().Rho(); // put largest distance first as cosmic rays (generally) go outside-in (downwards)
    }
  }crvsectorsort;

  void KinKalGeomMaker::makeTracker() {
      // surfaces need to match with virtual detectors. The following is extracted from VirtualDetectorMaker and needs to be updated if that changes.
    // Note that these are placed at the center of the VDs, which have half-thickness of 0.01mm. Since the VD hits are recorded where the SimParticle
    // enters the volume, the reco track will be sampled at a different position depending on the track direction by that amount. This is a fundamental
    // discrepancy between reco and sim data
    auto const& tracker = *(GeomHandle<mu2e::Tracker>());
    auto const& g4tmom = tracker.g4Tracker()->mother();
    auto const& ds = *(GeomHandle<DetectorSolenoid>());
    double vdHL(0.01); // hardcoded in VirtualDetectorMaker line 56
    // below are from VirtualDetectorMaker lnes 241-244
    double zFrontGlobal = g4tmom.position().z()-g4tmom.tubsParams().zHalfLength()-vdHL;
    double zBackGlobal  = g4tmom.position().z()+g4tmom.tubsParams().zHalfLength()+vdHL;
    // the 0.4 below comes from offsets in the mother volume nesting.
    double zFrontLocal  = zFrontGlobal - tracker.g4Tracker()->z0() + 0.4;
    double zBackLocal   = zBackGlobal  - tracker.g4Tracker()->z0() - 0.4;
    double zMidLocal = 10.1; // 10.1 is hard-coded in VirtualDetectorMaker line 224
    double halfLen = 0.5*(zBackLocal-zFrontLocal);
    double orvd = g4tmom.tubsParams().outerRadius();
    double irvd = tracker.g4Tracker()->getInnerTrackerEnvelopeParams().innerRadius();
    double irds = ds.rIn1();
    // cylinders are defined by TT_outer (_inner) virtual detectors
    // Disks are defined to match TT_front (mid, back) virtual detectors
    auto outer = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),orvd,halfLen);
    auto inner = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,zMidLocal),irvd,halfLen);
    // expand the disk radii to the DS
    auto front = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zFrontLocal),irds);
    auto mid = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zMidLocal),irds);
    auto back = std::make_shared<Disk>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,zBackLocal),irds);
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
    GeomHandle<DetectorSystem> det;
    GeomHandle<DetectorSolenoid> ds;
//    std::cout << "DS Cryo or " << ds->rOut1() << "," << ds->rOut2() << " ir " << ds->rIn1()<<","<< ds->rIn2() << " halfl " << ds->halfLength()
//      << " zpos " << ds->position().z() << " material " << ds->material() << std::endl;
//    std::cout << "DS shield or " << ds->shield_rOut1() << "," << ds->shield_rOut2() << " ir " << ds->shield_rIn1()<<","<< ds->shield_rIn2() << " halfl " << ds->shield_halfLength() << " material " << ds->shield_material() << std::endl;
//    std::cout << "DS ncoils " << ds->nCoils() << std::endl;
//    for(size_t icoil = 0; icoil < static_cast<size_t>(ds->nCoils()); icoil++){
//      std::cout << "DS coil ir " << ds->coil_rIn() << " or " << ds->coil_rOut()[icoil] << " length " << ds->coil_zLength()[icoil] << " zpos " << ds->coil_zPosition()[icoil]
//        << " material " << ds->coil_materials()[icoil] << std::endl;
//    }
    //DS Cryo or 1303,1328 ir 950,969.05 halfl 5450 zpos 8689 material StainlessSteel
    //DS shield or 1237.3,1250 ir 1010,1022.7 halfl 5287.7 material G4_Al
    //DS ncoils 11
    //DS coil ir 1050 or 1091 length 419.75 zpos 3539 material DS1CoilMix
    //DS coil ir 1050 or 1091 length 419.75 zpos 3964 material DS1CoilMix
    //DS coil ir 1050 or 1091 length 419.75 zpos 4389 material DS1CoilMix
    //DS coil ir 1050 or 1091 length 419.75 zpos 5042 material DS1CoilMix
    //DS coil ir 1050 or 1091 length 362.25 zpos 5699 material DS1CoilMix
    //DS coil ir 1050 or 1091 length 362.25 zpos 6369 material DS1CoilMix
    //DS coil ir 1050 or 1091 length 362.25 zpos 7176 material DS1CoilMix
    //DS coil ir 1050 or 1070.5 length 1777.5 zpos 7949 material DS2CoilMix
    //DS coil ir 1050 or 1070.5 length 1777.5 zpos 9761 material DS2CoilMix
    //DS coil ir 1050 or 1070.5 length 1777.5 zpos 11544 material DS2CoilMix
    //DS coil ir 1050 or 1091 length 362.25 zpos 13326 material DS1CoilMix
    //
    //
    //
    //
    auto inner= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),ds->rIn1(),ds->halfLength());
    auto outer= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-1482),ds->rOut2(),ds->halfLength()); // bounding surfaces
    auto front= std::make_shared<Disk>(outer->frontDisk());
    auto back= std::make_shared<Disk>(outer->backDisk());
    KKGeom::DetectorSolenoid::MaterialCylinderCollection materialCylinders;
    auto toDetectorZ = [&det](double zmu2e) {
      return VEC3(det->toDetector(CLHEP::Hep3Vector(0.0,0.0,zmu2e))).Z();
    };
    auto addMaterialCylinder = [this,&materialCylinders,&toDetectorZ](SurfaceId const& sid,
        double rin, double rout, double zcenter, double halfLength, std::string const& material) {
      auto cylinder = std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),
          VEC3(0.0,0.0,toDetectorZ(zcenter)),0.5*(rin+rout),halfLength);
      materialCylinders.emplace_back(sid,cylinder,material,rout-rin);
      kkg_->map_.emplace(std::make_pair(sid,std::static_pointer_cast<Surface>(cylinder)));
    };
    // KinKal DetMaterial name for the DS thermal shield (aluminium). The G4 geometry calls this
    // material "G4_Al"; here we use a DS-specific MatEnv material so the passive plane self-documents.
    std::string const dsShieldMaterial = "DSShieldAl";
    addMaterialCylinder(SurfaceId(SurfaceIdEnum::DS_CryoInner),ds->rIn1(),ds->rIn2(),ds->position().z(),ds->halfLength(),ds->material());
    addMaterialCylinder(SurfaceId(SurfaceIdEnum::DS_CryoOuter),ds->rOut1(),ds->rOut2(),ds->position().z(),ds->halfLength(),ds->material());
    addMaterialCylinder(SurfaceId(SurfaceIdEnum::DS_ShieldInner),ds->shield_rIn1(),ds->shield_rIn2(),ds->position().z(),ds->shield_halfLength(),dsShieldMaterial);
    addMaterialCylinder(SurfaceId(SurfaceIdEnum::DS_ShieldOuter),ds->shield_rOut1(),ds->shield_rOut2(),ds->position().z(),ds->shield_halfLength(),dsShieldMaterial);
    // Coils: collapse the 11 short coil cylinders into ONE long cylinder spanning the coil z-range.
    // NB: the DSCoilMix density (TrackerConditions/data/MaterialsList.data) assumes the effective thickness below;
    // both follow the run-2 coil geometry.
    {
      size_t const ncoils = static_cast<size_t>(ds->nCoils());
      double zmin = 1e9, zmax = -1e9, tsum = 0.0, rsum = 0.0, zlentot = 0.0;
      for(size_t icoil = 0; icoil < ncoils; icoil++){
        double const zlen = ds->coil_zLength().at(icoil);
        double const z0   = ds->coil_zPosition().at(icoil);
        double const rout = ds->coil_rOut().at(icoil);
        zmin = std::min(zmin, z0); zmax = std::max(zmax, z0 + zlen);
        tsum += zlen * (rout - ds->coil_rIn());          // z-weighted radial thickness
        rsum += zlen * 0.5*(ds->coil_rIn() + rout);       // z-weighted mid radius
        zlentot += zlen;
      }
      double const teff = tsum / zlentot;                 // effective radial thickness (~28.08 mm)
      double const reff = rsum / zlentot;                 // effective mid radius
      addMaterialCylinder(SurfaceId(SurfaceIdEnum::DS_Coil),
          reff - 0.5*teff, reff + 0.5*teff, 0.5*(zmin + zmax), 0.5*(zmax - zmin), "DSCoilMix");
    }

    // DS cryostat + thermal-shield END WALLS (z-normal annuli). Tracks bound for the upstream/
    // downstream CRV sectors (CRV-U / CRV-D) leave the solenoid through its z-faces and cross the
    // cryostat and shield end walls; the radial shells above only catch tracks that exit radially.
    // Geometry from DetectorSolenoid: (shield_)endWallHalfLength().
    {
      double const dsz = ds->position().z();
      auto addEndWall = [this,&toDetectorZ](SurfaceId const& sid, double zmu2e, double rin,
          double rout, double halfThick, std::string const& material) {
        auto annulus = std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),
            VEC3(0.0,0.0,toDetectorZ(zmu2e)),rin,rout);
        kkg_->passiveMaterialPlanes_.emplace_back(sid,annulus,material,halfThick);
        kkg_->map_.emplace(std::make_pair(sid,std::static_pointer_cast<Surface>(annulus)));
      };
      double const cOff = ds->halfLength()        - ds->endWallHalfLength();        // cryo wall center inset
      double const sOff = ds->shield_halfLength() - ds->shield_endWallHalfLength(); // shield wall center inset
      // cryostat end walls (StainlessSteel)
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Front,1), dsz-cOff, ds->rIn1(), ds->rOut2(),
          ds->endWallHalfLength(), ds->material());
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Back,1),  dsz+cOff, ds->rIn1(), ds->rOut2(),
          ds->endWallHalfLength(), ds->material());
      // thermal-shield end walls (aluminium; see dsShieldMaterial above)
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Front,2), dsz-sOff, ds->shield_rIn1(), ds->shield_rOut2(),
          ds->shield_endWallHalfLength(), dsShieldMaterial);
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Back,2),  dsz+sOff, ds->shield_rIn1(), ds->shield_rOut2(),
          ds->shield_endWallHalfLength(), dsShieldMaterial);
    }

    // hard-coded for now
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

    kkg_->ds_ = std::make_unique<KKGeom::DetectorSolenoid>(inner, outer, front, back, ipa, ipafront, ipaback, opa, tsda, materialCylinders);
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

  void KinKalGeomMaker::makeCRV() {
    GeomHandle<CosmicRayShield> CRS;
    GeomHandle<DetectorSystem> det;
    auto const& shields = CRS->getCRSScintillatorShields();
    std::vector<KKCRVSector> sectors;
    // Strongback (aluminium support plate) thickness from the CRV geometry
    double strongBackThickness = 0.0;
    for(auto const& shield : shields){
      int const td = shield.getFirstBar().getBarDetail().getThicknessDirection();
      for(int im = 0; im < shield.nModules(); ++im){
        auto const& mod = shield.getModule(im);
        if(mod.nAluminumSheets() > 0){ strongBackThickness = 2.0*mod.getAluminumSheet(0).getHalfLengths().at(td); break; }
      }
      if(strongBackThickness != 0.0) break;
    }
    // loop over the shields (= sectors)
    for (auto const& shield : shields) {
      //
      // First find this shield's orientation; the first bar is enough for that
      //
      auto const& firstbar = shield.getFirstBar();
      auto fbarpos = VEC3(det->toDetector(firstbar.getPosition())); // convert to detector (tracker) coordinates and root vectors
      auto bardet = firstbar.getBarDetail();
      // normal (w) direction is the thickness direction. Make sure it points away from the tracker
      VEC3 wdir;
      switch(bardet.getThicknessDirection()) {
        case 0:
          wdir = VEC3(copysign(1.0,fbarpos.X()),0.0,0.0);
          break;
        case 1:
          wdir = VEC3(0.0,copysign(1.0,fbarpos.Y()),0.0);
          break;
        case 2:
          wdir = VEC3(0.0,0.0,copysign(1.0,fbarpos.Z()));
          break;
        default:
          throw cet::exception("Service")<<"invalid direction "<< bardet.getThicknessDirection() << std::endl;
          break;
      }
      // u direction points along the bars (length direction). Sign is unimportant.
      VEC3 udir;
      switch(bardet.getLengthDirection()) {
        case 0:
          udir = VEC3(1.0,0.0,0.0);
          break;
        case 1:
          udir = VEC3(0.0,1.0,0.0);
          break;
        case 2:
          udir = VEC3(0.0,0.0,1.0);
          break;
        default:
          throw cet::exception("Service")<<"invalid direction "<< bardet.getLengthDirection() << std::endl;
          break;
      }
      // v points along bar width
      VEC3 vdir;
      switch(bardet.getWidthDirection()) {
        case 0:
          vdir = VEC3(1.0,0.0,0.0);
          break;
        case 1:
          vdir = VEC3(0.0,1.0,0.0);
          break;
        case 2:
          vdir = VEC3(0.0,0.0,1.0);
          break;
        default:
          throw cet::exception("Service")<<"invalid direction "<< firstbar.getBarDetail().getWidthDirection() << std::endl;
          break;
      }
      // next compute the average position. All the bars have the same position along their length
      double upos = fbarpos.Dot(udir);
      double uhw = firstbar.getHalfLength();
      // average first and last layers to get the w position and half-width
      auto const& firstmod = shield.getModule(0);
      auto flaypos = VEC3(det->toDetector(firstmod.getLayer(0).getPosition()));
      auto llaypos = VEC3(det->toDetector(firstmod.getLayer(firstmod.nLayers()-1).getPosition()));
      double wpos = 0.5*(flaypos+llaypos).Dot(wdir);
      // wdir can point opposite the layer-stacking order, depending on the sector.
      double whw = 0.5*std::abs((llaypos-flaypos).Dot(wdir)) + firstbar.getHalfThickness();
      // include the layer stagger when computing the position and width perp to the bars
      auto nlay = firstmod.nLayers();
      auto nbar = firstmod.getLayer(0).nBars();
      auto const& lastmod = shield.getModule(shield.nModules()-1);
      auto vf0 = VEC3(det->toDetector(firstmod.getLayer(0).getBar(0).getPosition())).Dot(vdir);
      auto vf3 = VEC3(det->toDetector(firstmod.getLayer(nlay-1).getBar(0).getPosition())).Dot(vdir);
      auto vl0 = VEC3(det->toDetector(lastmod.getLayer(0).getBar(nbar-1).getPosition())).Dot(vdir);
      auto vl3 = VEC3(det->toDetector(lastmod.getLayer(nlay-1).getBar(nbar-1).getPosition())).Dot(vdir);
      double vpos = 0.25*(vf0+vf3+vl0+vl3);
      double vmin = std::min({vf0,vf3,vl0,vl3});
      double vmax = std::max({vf0,vf3,vl0,vl3});
      double vhw = 0.5*(vmax-vmin)+ firstbar.getHalfWidth();
      VEC3 midpoint = upos*udir + vpos*vdir + wpos*wdir;
      // create the rectangle
      KKCRVSector sector;
      sector.sname_ = shield.getName();
      sector.sector_ = std::make_shared<KinKal::Rectangle>(wdir,udir,midpoint,uhw,vhw);
      sector.whw_ = whw;
      // Each CRV module carries a ~12.7 mm aluminium strongback (support plate) on the tracker-facing side
      // of the scintillator stack. Thickness (strongBackThickness, above) and material now come from the
      // CRV geometry rather than crs.strongBackThickness / crs.aluminumSheetMaterialName.
      // KinKal DetMaterial name for the strongback (aluminium). The CRV geometry's aluminium-sheet
      // material is the G4 name "G4_Al"; we use a CRV-specific MatEnv material so the plane self-documents.
      std::string const strongBackMaterial = "CRVStrongbackAl";
      auto const sbcenter = midpoint - (whw + 0.5*strongBackThickness)*wdir;
      auto const sbplane  = std::make_shared<KinKal::Rectangle>(wdir,udir,sbcenter,uhw,vhw);
      SurfaceId sbsid(SurfaceIdEnum::CRV_StrongBack, static_cast<int>(kkg_->passiveMaterialPlanes_.size()));
      kkg_->passiveMaterialPlanes_.emplace_back(sbsid, sbplane, strongBackMaterial, 0.5*strongBackThickness);
      kkg_->map_.emplace(std::make_pair(sbsid, std::static_pointer_cast<Surface>(sbplane)));
      sectors.push_back(sector);
    }
    // sort the sectors according to their transverse distance (largest first), to optimize searching for downward going tracks.
    std::sort(sectors.begin(),sectors.end(),crvsectorsort);
    if(debug_ > 0){
      for(auto const& sector : sectors){
        std::cout << "CRV sector " <<  sector.sname_;
        auto const& sectptr = sector.sector_;
        std::cout << " midpoint " << sectptr->center() << " wdir " << sectptr->normal() << " udir " << sectptr->uDirection() << " vdir " << sectptr->vDirection()
          << " uhw " << sectptr->uHalfLength() << " vhw " << sectptr->vHalfLength() <<  " whw " << sector.whw_ << std::endl;
      }
    }

    kkg_->crv_ = std::make_unique<KKGeom::CRV>(sectors);
    // fill map
    unsigned isect(0);
    for(auto const& sector : kkg_->crv_->sectors()){
      // Tolerate CRV sector names with no SurfaceId enum entry: historical
      // geometries (e.g. geom_common_MDC2020 -> crv_counters_v09, sector
      // CRV_C3) predate the current enum. Such sectors remain part of the
      // CRV geometry but are not addressable as KinKal extrapolation
      // surfaces.
      SurfaceIdEnum sid(sector.sname_, false, false);
      if(sid.id() == SurfaceIdDetail::unknown){
        if(debug_ > 0) std::cout << "KinKalGeomMaker: CRV sector " << sector.sname_
          << " has no SurfaceId enum entry; skipping KinKal surface" << std::endl;
        continue;
      }
      kkg_->map_.emplace(std::make_pair(SurfaceId(sid.id()),std::static_pointer_cast<Surface>(sector.sector_)));
      isect++;
    }
  }

  // Append one averaged rectangular concrete passive-material plane for an ExtShield region.
  // 'halfThickness' is the half material depth a track crosses along the normal
  // (used for the energy-loss/scattering path length). All planes use the broad
  // DS_HatchConcrete SurfaceId; a running index distinguishes the individual planes.
  void KinKalGeomMaker::addConcretePlane(DetectorSystem const& det,
      int normalAxis, CLHEP::Hep3Vector const& centerMu2e,
      double hw1, double hw2, double halfThickness, std::string const& material) {
    // Build the plane's local frame. The Rectangle ctor takes (normal,udir,center,
    // uHalfLen,vHalfLen) with vdir = normal x udir (right-handed). We pick udir/vdir
    // as the two global axes orthogonal to the normal so the half-extents map cleanly.
    VEC3 norm, udir;
    double uhw, vhw;
    switch(normalAxis) {
      case 0: // x-normal (side walls): u=y, v=z
        norm = VEC3(1.0,0.0,0.0); udir = VEC3(0.0,1.0,0.0); uhw = hw1; vhw = hw2; break;
      case 1: // y-normal (roof): u=x, v=z (matches the Type-2 roof convention below)
        norm = VEC3(0.0,1.0,0.0); udir = VEC3(1.0,0.0,0.0); uhw = hw1; vhw = hw2; break;
      case 2: // z-normal (downstream endcap): u=x, v=y
        norm = VEC3(0.0,0.0,1.0); udir = VEC3(1.0,0.0,0.0); uhw = hw1; vhw = hw2; break;
      default:
        throw cet::exception("Service") << "invalid concrete-plane normal axis " << normalAxis << std::endl;
    }
    // Split the crossing into two Gauss-point planes along the normal, at center +/- halfThickness/sqrt3,
    // each carrying material half-depth 'halfThickness' (so the two together = the full slab depth
    // 2*halfThickness). This integrates the bulk multiple scattering with a 2-point Gauss rule (exact
    // for a uniform slab: reproduces the continuous angle-position covariance) while preserving the total
    // crossed material, so energy loss is correct.
    double const invSqrt3 = 1.0/std::sqrt(3.0);
    auto const center0 = VEC3(det.toDetector(centerMu2e));
    for(int iplane = 0; iplane < 2; ++iplane) {
      double const sign = iplane == 0 ? -1.0 : 1.0;
      auto const center = center0 + (sign*halfThickness*invSqrt3)*norm;
      auto const plane = std::make_shared<Rectangle>(norm, udir, center, uhw, vhw);
      SurfaceId sid(SurfaceIdEnum::DS_HatchConcrete,
                    static_cast<int>(kkg_->passiveMaterialPlanes_.size()));
      kkg_->passiveMaterialPlanes_.emplace_back(sid, plane, material, halfThickness);
      kkg_->map_.emplace(std::make_pair(sid, std::static_pointer_cast<Surface>(plane)));
    }
  }

  void KinKalGeomMaker::makePassiveMaterials() {
    GeomHandle<DetectorSystem> det;
    // External shielding (concrete) geometry now comes from the ExtShieldDownstream GeometryService object
    art::ServiceHandle<GeometryService> geomService;
    if(!geomService->hasElement<ExtShieldDownstream>()) return;
    GeomHandle<ExtShieldDownstream> esd;
    if(esd->getNumberOfBoxTypes() <= 0) return;
    // Group block indices by their 1-based type number. The type identity is what SimpleConfig
    // supplied and what outline/material/orientation alone cannot recover (several types share
    // identical geometry); getBoxTypes() restores it. Vectors below are parallel per block.
    auto const& esdTypes    = esd->getBoxTypes();
    auto const& esdOutlines = esd->getOutlines();     // per block: list of {u,v} vertices (mm)
    auto const& esdLengths  = esd->getLengths();       // per block: HALF length along W (mm)
    auto const& esdMats     = esd->getMaterialNames();
    auto const& esdCenters  = esd->getCentersOfBoxes();
    std::map<int,std::vector<std::size_t>> boxesByType;
    for(std::size_t i = 0; i < esdTypes.size(); ++i) boxesByType[esdTypes[i]].push_back(i);
    auto firstBox = [&](int t)->long {
      auto it = boxesByType.find(t);
      return (it == boxesByType.end() || it->second.empty()) ? -1L : static_cast<long>(it->second.front());
    };
    auto nBoxOf   = [&](int t){ auto it = boxesByType.find(t); return it == boxesByType.end() ? 0 : static_cast<int>(it->second.size()); };
    auto materialOf = [&](int t){ long const b = firstBox(t); return b >= 0 ? esdMats[b] : std::string(); };
    // Mu2e-frame center {x,y,z} of the ibox-th (1-based) block of a type (mirrors centerType{t}Box{ibox}).
    auto centerOf = [&](int t, int ibox)->std::vector<double> {
      auto const& c = esdCenters[boxesByType[t][ibox-1]];
      return { c.x(), c.y(), c.z() };
    };
    // --- ExtShieldDownstream box dimensions. A box "type" is an extruded U-V outline of half-length
    // esdLengths; the helpers return its half-extents. The {U,V,W=length} -> {x,y,z} orientation per
    // region is not carried here, so it stays explicit at the call sites (as before).
    // Helpers return 0 for absent types so an optional region is simply skipped (its nBox guard fires).
    auto vertsOf = [&](int t, bool uaxis){
      std::vector<double> v; long const b = firstBox(t);
      if(b >= 0) for(auto const& uv : esdOutlines[b]) v.push_back(uaxis ? uv[0] : uv[1]);
      return v;
    };
    auto maxAbs  = [](std::vector<double> const& v){ double m=0.0; for(double x : v) m=std::max(m,std::fabs(x)); return m; };
    auto spanHalf= [](std::vector<double> const& v){ if(v.empty()) return 0.0; auto e=std::minmax_element(v.begin(),v.end()); return 0.5*(*e.second - *e.first); };
    auto uHalf   = [&](int t){ return maxAbs(vertsOf(t,true)); };    // outline U half-width (max |U|)
    auto vHalf   = [&](int t){ return maxAbs(vertsOf(t,false)); };   // outline V half-height (max |V|)
    auto vSpanH  = [&](int t){ return spanHalf(vertsOf(t,false)); }; // V half-extent (handles 0-based outlines)
    auto lenHalf = [&](int t){ long const b = firstBox(t); return b >= 0 ? esdLengths[b] : 0.0; };  // already half
    // V-weighted mean half-thickness of a T cross-section along its U (thin) axis: used where a
    // T-block is collapsed to a single slab whose normal is U (the side walls), giving the
    // energy-loss-correct mean depth without per-component Gauss planes. Reproduces the run-2 449.7.
    auto tMeanUHalf = [&](int t){
      auto uV=vertsOf(t,true), vV=vertsOf(t,false);
      if(uV.empty() || vV.empty()) return 0.0;
      double const uo=maxAbs(uV); double ui=uo;
      for(double u : uV) if(u > 0.0 && u < uo-1.0) ui=std::min(ui,u);   // inner (stem) U half-width
      auto const ve=std::minmax_element(vV.begin(),vV.end());
      double const vmin=*ve.first, vmax=*ve.second; double vsh=0.5*(vmin+vmax);
      for(double v : vV) if(v > vmin+1.0 && v < vmax-1.0){ vsh=v; break; } // shoulder V
      double const cb=vsh-vmin, st=vmax-vsh, span=vmax-vmin;             // crossbar / stem V extents
      return span>0.0 ? (uo*cb + ui*st)/span : uo;
    };

    // Model the ExtShieldDownstream Type 2 regular concrete roof T-blocks.
    // These blocks sit between the DS outer cryostat (y~1328 mm) and the CRV top sectors T1/T2 (y~2653 mm).
    // With orientations "010"/"012" (both map U->z, V->y, W->x):
    //   U (±uhw_outer mm) -> ±z footprint of the wide crossbar per block
    //   V range [vmin, vmax]  -> y; vshoulder separates crossbar from stem
    //   W (length/2 = xhw mm) -> x (transverse, same for both components)
    //
    // The T cross-section is decomposed into two non-overlapping rectangles:
    //   Crossbar: V in [vmin, vshoulder], U in [-uhw_outer, +uhw_outer]
    //   Stem    : V in [vshoulder, vmax], U in [-uhw_inner, +uhw_inner]
    // Each is modelled with two Gauss-point planes so tracks through the crossbar-
    // only region (the ~2/3 of z-width outside the stem) correctly see half the
    // concrete compared to tracks through the full T height.
    // Wrapped in an IIFE so a degenerate Type-2 (no boxes / empty outline) skips only this block
    // via `return`, without abandoning the side walls / roof / endcap / hatch handled below.
    [&]{
      int const nType2 = nBoxOf(2);
      if(nType2 <= 0) return;
      double const invSqrt3 = 1.0/std::sqrt(3.0);
      std::vector<double> uVerts = vertsOf(2,true), vVerts = vertsOf(2,false);
      if(uVerts.empty() || vVerts.empty()) return;

      auto const uExt = std::minmax_element(uVerts.begin(), uVerts.end());
      auto const vExt = std::minmax_element(vVerts.begin(), vVerts.end());
      double const xhw       = lenHalf(2);
      auto const material    = materialOf(2);

      double const vmin_t    = *vExt.first;   // = -452.2 mm (bottom of crossbar)
      double const vmax_t    = *vExt.second;  // = +452.2 mm (top of stem)
      double const uhw_outer = *uExt.second;  // = 680.8 mm  (crossbar z half-width per block)

      // Inner U half-width (stem): smallest positive U vertex that differs from uhw_outer
      double uhw_inner = uhw_outer;
      for(double u : uVerts)
        if(u > 0 && u < uhw_outer - 1.0) uhw_inner = std::min(uhw_inner, u);

      // Shoulder V (where T width steps inward): V value strictly between vmin and vmax
      double vshoulder = 0.0;
      for(double v : vVerts)
        if(v > vmin_t + 1.0 && v < vmax_t - 1.0) { vshoulder = v; break; }

      // Crossbar component: V in [vmin_t, vshoulder]
      double const yhalf_cb     = 0.5*(vshoulder - vmin_t);   // = 223.6 mm
      double const vcenter_cb   = 0.5*(vmin_t + vshoulder);   // = -228.6 mm (local V offset)

      // Stem component: V in [vshoulder, vmax_t]
      double const yhalf_stem   = 0.5*(vmax_t - vshoulder);   // = 228.6 mm
      double const vcenter_stem = 0.5*(vshoulder + vmax_t);   // = +223.6 mm (local V offset)

      // Accumulate x/y center (common to both components) and z bounding boxes
      // (different per component because stem is narrower in z than crossbar).
      double xcenter = 0.0, ycenter = 0.0;
      double zmin_cb = 1.0e9, zmax_cb = -1.0e9;
      double zmin_st = 1.0e9, zmax_st = -1.0e9;
      int nfound = 0;
      for(int ibox = 1; ibox <= nType2; ++ibox) {
        std::vector<double> ctr = centerOf(2, ibox);
        xcenter += ctr[0];
        ycenter += ctr[1];
        zmin_cb = std::min(zmin_cb, ctr[2] - uhw_outer);
        zmax_cb = std::max(zmax_cb, ctr[2] + uhw_outer);
        zmin_st = std::min(zmin_st, ctr[2] - uhw_inner);
        zmax_st = std::max(zmax_st, ctr[2] + uhw_inner);
        ++nfound;
      }
      if(nfound == 0) return;
      xcenter /= nfound;
      ycenter /= nfound;

      // Global Mu2e Y centers for each T component (ycenter = mean box center = 2006.1 mm)
      double const ycenter_cb   = ycenter + vcenter_cb;    // ≈ 1777.5 mm
      double const ycenter_stem = ycenter + vcenter_stem;  // ≈ 2229.7 mm

      double const zcenters[2] = { 0.5*(zmin_cb + zmax_cb), 0.5*(zmin_st + zmax_st) };
      double const zhws[2]     = { 0.5*(zmax_cb - zmin_cb), 0.5*(zmax_st - zmin_st) };
      double const yctrs[2]    = { ycenter_cb,   ycenter_stem  };
      double const yhalfs[2]   = { yhalf_cb,     yhalf_stem    };

      // Four planes total: two Gauss-point planes for each T component.
      for(int icomp = 0; icomp < 2; ++icomp) {
        for(int iplane = 0; iplane < 2; ++iplane) {
          double const ysign = iplane == 0 ? -1.0 : 1.0;
          auto const center = VEC3(det->toDetector(
              CLHEP::Hep3Vector(xcenter, yctrs[icomp] + ysign*yhalfs[icomp]*invSqrt3, zcenters[icomp])));
          auto const plane = std::make_shared<Rectangle>(
              VEC3(0.0,1.0,0.0), VEC3(1.0,0.0,0.0), center, xhw, zhws[icomp]);
          SurfaceId sid(SurfaceIdEnum::DS_HatchConcrete,
                        static_cast<int>(kkg_->passiveMaterialPlanes_.size()));
          kkg_->passiveMaterialPlanes_.emplace_back(sid, plane, material, yhalfs[icomp]);
          kkg_->map_.emplace(std::make_pair(sid, std::static_pointer_cast<Surface>(plane)));
        }
      }
    }();

    // Collapse all type-<boxType> ExtShieldDownstream concrete boxes of a region into one
    // averaged plane. 'normalAxis' is the slab normal (0=x,1=y,2=z); the in-plane
    // half-extents are (center span)/2 + the block half-size along each in-plane
    // axis; 'halfThickness' is supplied by the caller (block depth along normal).
    // Boxes whose center lies far from the region cluster (different sub-region,
    // e.g. the type-1 north vs south wall) are separated by an explicit center cut.
    auto buildAveragedRegion = [&](int boxType, int normalAxis,
        double halfU_block, double halfV_block, double halfThickness,
        std::function<bool(std::vector<double> const&)> select) {
      int const nbox = nBoxOf(boxType);
      if(nbox <= 0) return;
      auto const material = materialOf(boxType);
      // axis indices orthogonal to the normal, in increasing order
      int a1 = (normalAxis == 0) ? 1 : 0;
      int a2 = (normalAxis == 2) ? 1 : 2;
      double nsum = 0.0;              // accumulate the normal-axis center (slab position)
      double min1 = 1e9, max1 = -1e9; // in-plane axis 1 center range
      double min2 = 1e9, max2 = -1e9; // in-plane axis 2 center range
      int nfound = 0;
      for(int ibox = 1; ibox <= nbox; ++ibox) {
        std::vector<double> ctr = centerOf(boxType, ibox);
        if(select && !select(ctr)) continue;
        nsum += ctr[normalAxis];
        min1 = std::min(min1, ctr[a1]); max1 = std::max(max1, ctr[a1]);
        min2 = std::min(min2, ctr[a2]); max2 = std::max(max2, ctr[a2]);
        ++nfound;
      }
      if(nfound == 0) return;
      // assemble the center from the per-axis values (normal axis = mean slab
      // position; in-plane axes = midpoint of the center bounding box)
      double cc[3] = {0.0,0.0,0.0};
      cc[normalAxis] = nsum / nfound;
      cc[a1] = 0.5*(min1 + max1);
      cc[a2] = 0.5*(min2 + max2);
      CLHEP::Hep3Vector center(cc[0],cc[1],cc[2]);
      double const hw1 = 0.5*(max1 - min1) + halfU_block;
      double const hw2 = 0.5*(max2 - min2) + halfV_block;
      this->addConcretePlane(*det, normalAxis, center, hw1, hw2, halfThickness, material);
    };

    // --- Side walls. Type 1 regular-concrete wall T-blocks. Geometry from
    // ExtShieldDownstream (outlineType1 U/V verts + lengthType1): the T cross-section
    // lies in the U-V plane, extruded over lengthType1; for the walls the orientation
    // maps U->x (wall thickness), V->z (footprint), length->y (vertical extent).
    //   U verts +/-680.8 / +/-223.6  -> crossbar 1361.6 mm thick, stem 447.2 mm thick
    //   V verts +/-452.2 (span 904.4); crossbar occupies V in [-452.2,-5] = 447.2 mm,
    //                                  stem     occupies V in [-5, 452.2] = 457.2 mm
    // The wall is one x-normal slab per side whose x half-thickness is the z-weighted
    // mean of the T thickness, so the mean crossed concrete (energy loss) is exact
    // without four Gauss planes per wall:
    //   2*xHalf = (1361.6*447.2 + 447.2*457.2)/904.4 = 899.3 mm  ->  xHalfThick = 449.7
    // Box centers split into the north wall (x ~ -1897) and south wall (x ~ -5911)
    // about the DS axis (x = -3904); each is emitted as a single slab.
    {
      double const yHalfBlock = lenHalf(1);     // length (W) -> y vertical half-height
      double const zHalfBlock = vSpanH(1);      // V outline -> z footprint (452.2 in run-2)
      double const xHalfThick = tMeanUHalf(1);  // U-outline T-mean -> x half-thickness (449.7 in run-2)
      // Separate the two symmetric side walls at the DS axis = the mean x of the type-1 box
      // centers (-3904 in run-2), derived from config so it tracks the geometry instead of a literal.
      double splitX = 0.0; int nWall = 0;
      int const nBox1 = nBoxOf(1);
      for(int ibox = 1; ibox <= nBox1; ++ibox) {
        std::vector<double> c = centerOf(1, ibox);
        splitX += c[0]; ++nWall;
      }
      if(nWall > 0) {
        splitX /= nWall;
        buildAveragedRegion(1, 0, yHalfBlock, zHalfBlock, xHalfThick,
            [splitX](std::vector<double> const& c){ return c[0] > splitX; }); // north
        buildAveragedRegion(1, 0, yHalfBlock, zHalfBlock, xHalfThick,
            [splitX](std::vector<double> const& c){ return c[0] < splitX; }); // south
      }
    }

    // --- Other roof concrete types (same y-normal plane family as Type-2). These
    // sit at y ~ 2006 between the DS and the top CRV. We add
    // them as averaged y-normal slabs using the roof T half-height (the V outline
    // ~452.2 mm maps to the vertical thickness for the roof). Types: 4 (barite roof
    // T, 6 boxes), 9 (barite roof L, 2), 10 (concrete roof filler, 1), 23 (barite
    // roof T around ST, 2), 29 (endcap top-back, 1), 30 (endcap top-front, 1).
    // For these the wide U crossbar (680.8) and the W length (4918 for the long
    // roof types) span x and z; we use the conservative full-width block half-sizes
    // so the averaged footprint covers the tiled blocks. y half-thickness = 452.2.
    {
      double const roofYhalfThick = vHalf(2);   // V outline half -> roof vertical depth (452.2)
      double const roofXhalfBlock = lenHalf(2); // length/2 spans x (~2459)
      double const roofZhalfBlock = uHalf(2);   // U crossbar half spans z (680.8)
      // These slabs borrow Type-2's dimensions; if Type-2 is zeroed but these types aren't, skip
      // rather than emit zero-thickness planes (a partial config the run-1 nBoxOf(2)<=0 path misses).
      if(roofYhalfThick > 0.0) {
      // type 4: barite roof T blocks (upstream)
      buildAveragedRegion(4, 1, roofXhalfBlock, roofZhalfBlock, roofYhalfThick, nullptr);
      // type 23: barite roof T blocks around the stopping target
      buildAveragedRegion(23, 1, roofXhalfBlock, roofZhalfBlock, roofYhalfThick, nullptr);
      // type 9: barite roof L blocks on the ends of the upstream roof
      buildAveragedRegion(9, 1, roofXhalfBlock, roofZhalfBlock, roofYhalfThick, nullptr);
      // type 10: concrete roof filler block
      buildAveragedRegion(10, 1, roofXhalfBlock, roofZhalfBlock, roofYhalfThick, nullptr);
      // types 29/30: new-endcap top back / front roof pieces. y-normal slabs; their own wide U
      // outline spans x and their short length spans z, modelled at the roof slab depth.
      buildAveragedRegion(29, 1, uHalf(29), lenHalf(29), roofYhalfThick, nullptr);
      buildAveragedRegion(30, 1, uHalf(30), lenHalf(30), roofYhalfThick, nullptr);
      }
    }

    // --- Downstream endcap concrete (z-normal). Types 27 (3 boxes) and 28 (1 box)
    // form the wall closing the DS<->CRV gap at z ~ 18542, beyond the DS, at the
    // CRV-D sectors. They have large x/y outlines (U up to 2463.8, V up to 3246.2,
    // the V outline runs 0..3246.2 so the per-block V half-extent is ~1623 mm) and
    // a short z length (~1117 / ~907 mm), i.e. a slab normal to z. We emit one
    // averaged z-normal plane per type, with x half from the U outline, y half from
    // the centers' bounding box plus the block V half, and z half-thickness from
    // length/2. NOTE: the asymmetric (0-based) V outline means the precise y
    // placement of the face would need the per-box "300" rotation to be decoded;
    // here we cover the whole stacked x-y face conservatively, which is adequate for
    // the broad energy-loss correction this region needs for CRV-D-bound tracks.
    {
      // Endcap face (z-normal): the type-27 U outline -> x, its V outline span -> y; the slab
      // z half-thickness is each type's own length/2. (Type 28 shares the type-27 x-y face, as before.)
      double const ecXhalf = uHalf(27);   // U outline -> x half (2463.8)
      double const ecYhalf = vSpanH(27);  // V outline span -> y half (~1623)
      // Both endcap types use Type-27's x-y face; skip if Type-27 is zeroed so Type-28 doesn't
      // emit a zero-area plane.
      if(ecXhalf > 0.0 && ecYhalf > 0.0) {
      // type 27: 3 stacked endcap blocks (cover most of the x-y face)
      buildAveragedRegion(27, 2, ecXhalf, ecYhalf, lenHalf(27), nullptr);
      // type 28: endcap bottom block
      buildAveragedRegion(28, 2, ecXhalf, ecYhalf, lenHalf(28), nullptr);
      }
    }

    // Run-1 geometries zero the ExtShieldDownstream concrete T-blocks but retain the building
    // dsAreaHatchBlock concrete (CONCRETE_MARS) between tracker and CRV-T. Extracted geometries
    // have no ExtShieldDownstream object and returned above. Reaching here means the shield is
    // declared (numberOfBoxTypes>0); the run-1 fallback fires when its roof T-blocks are zeroed.
    if(nBoxOf(2) <= 0) {
      addBuildingHatchConcrete(*det);
    }
  }

  void KinKalGeomMaker::addBuildingHatchConcrete(DetectorSystem const& det) {
    // The dsAreaHatchBlock concrete is a building ExtrudedSolid served by Mu2eHall, replacing the
    // former building.dsArea.hatchblock.* / yOfFloorSurface SimpleConfig reads. The solid's offset
    // already folds the floor-surface Y into its Mu2e-origin offset. The 2D outline lies in the
    // (x, z) plane: vertex x -> Mu2e x, vertex y -> Mu2e z; the slab normal is y (a roof plane).
    GeomHandle<Mu2eHall> hall;
    auto const& solids = hall->getBldgSolids();
    auto spanHalf = [](std::vector<CLHEP::Hep2Vector> const& v, bool xaxis) {
      double mn = 1e30, mx = -1e30;
      for(auto const& p : v){ double const c = xaxis ? p.x() : p.y(); mn = std::min(mn,c); mx = std::max(mx,c); }
      return 0.5*(mx - mn);
    };
    auto mid = [](std::vector<CLHEP::Hep2Vector> const& v, bool xaxis) {
      double mn = 1e30, mx = -1e30;
      for(auto const& p : v){ double const c = xaxis ? p.x() : p.y(); mn = std::min(mn,c); mx = std::max(mx,c); }
      return 0.5*(mn + mx);
    };
    auto addSlab = [&](std::string const& name) {
      auto const it = solids.find(name);
      if(it == solids.end()) return;
      auto const& solid = it->second;
      auto const& verts = solid.getVertices();
      if(verts.empty()) return;
      auto const& off = solid.getOffsetFromMu2eOrigin();
      CLHEP::Hep3Vector center(off.x() + mid(verts,true), off.y(), off.z() + mid(verts,false));
      addConcretePlane(det, 1, center, spanHalf(verts,true), spanHalf(verts,false),
                       solid.getYhalfThickness(), solid.getMaterial());
    };
    addSlab("dsAreaHatchBlock");
    addSlab("dsAreaHatchBlockE");
  }
}
