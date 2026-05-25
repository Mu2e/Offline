//
// Create KinKalGeom objects. These depend on other objects served by GeometryService so this must be external to KinKalGeom itself
// Original author: Dave Brown (LBNL) 4/2026
//
#include "Offline/GeometryService/inc/KinKalGeomMaker.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/Tracker.hh"
#include "Offline/KinKalGeom/inc/DetectorSolenoid.hh"
#include "Offline/KinKalGeom/inc/StoppingTarget.hh"
#include "Offline/KinKalGeom/inc/CRV.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "cetlib_except/exception.h"
#include <cmath>
#include <algorithm>

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
    return kkg_;
  }

  // sort by transverse distance
  struct sortCRVSectors {
    bool operator () (KKCRVSector const& sect1, KKCRVSector const& sect2) {
      return sect1.sector_->center().Rho() < sect2.sector_->center().Rho();
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

  void KinKalGeomMaker::makeCRV() {
    GeomHandle<CosmicRayShield> CRS;
    GeomHandle<DetectorSystem> det;
    auto const& shields = CRS->getCRSScintillatorShields();
    std::vector<KKCRVSector> sectors;
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
      double whw = 0.5*(llaypos-flaypos).Dot(wdir) + firstbar.getHalfThickness();
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
      sectors.push_back(sector);
    }
    // sort the sectors according to their transverse distance
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
      kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::CRV,isect),std::static_pointer_cast<Surface>(sector.sector_)));
      isect++;
    }
  }
}
