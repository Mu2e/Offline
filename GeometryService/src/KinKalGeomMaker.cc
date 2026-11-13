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
#include "Offline/KinKalGeom/inc/Calo.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "Offline/GeneralUtilities/inc/OrientationResolver.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <functional>

namespace mu2e {
  using KinKal::VEC3;
  using KinKal::Cylinder;
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
    makeCalo();
    return kkg_;
  }

  void KinKalGeomMaker::makeCalo() {
    // construct calorimeter geometry from GeometryService
    art::ServiceHandle<GeometryService> geomService;
    if(!geomService->hasElement<Calorimeter>()) return;
    auto const& calo_det = *(GeomHandle<mu2e::Calorimeter>());
    GeomHandle<DetectorSystem> det;
    auto toDetectorZ = [&det](double zmu2e) {
      return VEC3(det->toDetector(CLHEP::Hep3Vector(0.0,0.0,zmu2e))).Z();
    };

    // Extract geometry from the first two disks (D0 and D1)
    auto const& geom0 = calo_det.disk(0).diskInfo();
    auto const& geom1 = calo_det.disk(1).diskInfo();
    double r0_inner = geom0.innerEnvelopeR();
    double r0_outer = geom0.outerEnvelopeR();
    double r1_inner = geom1.innerEnvelopeR();
    double r1_outer = geom1.outerEnvelopeR();

    // disk front/back faces, converted from Mu2e to detector coordinates
    double z0_front = toDetectorZ(geom0.origin().z() - 0.5*geom0.size().z());
    double z0_back  = toDetectorZ(geom0.origin().z() + 0.5*geom0.size().z());
    double z1_front = toDetectorZ(geom1.origin().z() - 0.5*geom1.size().z());
    double z1_back  = toDetectorZ(geom1.origin().z() + 0.5*geom1.size().z());

    // Construct Calo with geometry in detector coordinates
    kkg_->calo_ = std::make_unique<KKGeom::Calo>(r0_inner, r0_outer, r1_inner, r1_outer, z0_front, z0_back, z1_front, z1_back);
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

    // Crystal case inner support (see caloBuildCase in Mu2eG4/src/constructDiskCalorimeter.cc): the case runs from the
    // crystal caps to the back of the crystals, and holds solid aluminium rings at its front, middle and back, between
    // the case inner radius and the carbon fiber ring that separates them from the crystals
    auto const& g4info = calo_det.G4Info();
    double const capthick = g4info.get<double>("crystalCapZLength");
    double const caselen = g4info.get<double>("diskCaseZLength");
    double const alrin = g4info.get<double>("diskInAlRingRIn");
    double const alzlen = g4info.get<double>("diskInAlRingZLength");
    double const cfrin = g4info.get<double>("diskInCFRingRIn");
    double const cfrout = g4info.get<double>("diskInCFRingROut");
    // KinKal material for the aluminium (G4_Al) rings, named per subsystem like the front panel pipes
    std::array<std::string,KKGeom::CaloInnerSupport::nmaterials> const ringmats = {"CaloRingAl",g4info.get<std::string>("innerCFRingMaterial")};
    static const SurfaceIdEnum::enum_type ringsid[2] = {SurfaceIdEnum::EMC_Disk_0_InnerRing, SurfaceIdEnum::EMC_Disk_1_InnerRing};
    for(int idisk=0; idisk < 2; ++idisk){
      using Support = KKGeom::CaloInnerSupport;
      double const zfront = toDetectorZ(calo_det.disk(idisk).diskInfo().frontFaceCenter().z() - capthick);
      double const zback = zfront + caselen;
      double const zmid = 0.5*(zfront+zback);
      auto const center = det->toDetector(calo_det.disk(idisk).diskInfo().origin());
      Support support(ringsid[idisk],VEC3(center.x(),center.y(),zfront),alrin,cfrin,cfrout,zfront,zback,ringmats);
      support.addRing(Support::Ring{alrin,cfrin,zfront,zfront+alzlen,Support::aluminium});
      support.addRing(Support::Ring{alrin,cfrin,zmid-0.5*alzlen,zmid+0.5*alzlen,Support::aluminium});
      support.addRing(Support::Ring{alrin,cfrin,zback-alzlen,zback,Support::aluminium});
      support.addRing(Support::Ring{cfrin,cfrout,zfront,zback,Support::carbon});
      kkg_->map_.emplace(std::make_pair(support.surfaceId(Support::alFront),std::static_pointer_cast<Surface>(support.front())));
      kkg_->map_.emplace(std::make_pair(support.surfaceId(Support::alHole),std::static_pointer_cast<Surface>(support.inner())));
      kkg_->map_.emplace(std::make_pair(support.surfaceId(Support::cfRing),std::static_pointer_cast<Surface>(support.cfInner())));
      kkg_->map_.emplace(std::make_pair(support.surfaceId(Support::cfFront),std::static_pointer_cast<Surface>(support.front())));
      kkg_->map_.emplace(std::make_pair(support.surfaceId(Support::alRing),std::static_pointer_cast<Surface>(support.cfInner())));
      kkg_->caloInnerSupports_.push_back(std::move(support));
    }

    // Front panel passive material: a carbon fiber / foam / carbon fiber sandwich directly upstream of the
    // crystal caps (see caloBuildFrontPlate in Mu2eG4/src/constructDiskCalorimeter.cc). Each layer is an
    // annulus carrying its full thickness; index 0 is the upstream carbon skin, 1 the foam, 2 the downstream skin.
    // The cooling pipes embedded in the foam (index 3) are described by their centerlines, see makeCaloPipes
    if(!g4info.get<bool>("hasFrontPanel")) return;
    double const fprin = g4info.get<double>("FPInnerRadius");
    double const fprout = g4info.get<double>("FPOuterRadius");
    double const carbonthick = g4info.get<double>("FPCarbonZLength");
    double const foamthick = g4info.get<double>("FPFoamZLength");
    std::string const carbonmat = g4info.get<std::string>("FPCarbonMaterial");
    std::string const foammat = g4info.get<std::string>("FPFoamMaterial");
    static const SurfaceIdEnum::enum_type fpsid[2] = {SurfaceIdEnum::EMC_Disk_0_FrontPanel, SurfaceIdEnum::EMC_Disk_1_FrontPanel};
    for(int idisk=0; idisk < 2; ++idisk){
      // the front panel ends where the crystal caps begin
      double const fpback = calo_det.disk(idisk).diskInfo().frontFaceCenter().z() - capthick;
      double const zlayer[3] = {fpback - 1.5*carbonthick - foamthick, fpback - carbonthick - 0.5*foamthick, fpback - 0.5*carbonthick};
      double const thick[3] = {carbonthick, foamthick, carbonthick};
      std::string const* mat[3] = {&carbonmat, &foammat, &carbonmat};
      for(int ilayer=0; ilayer < 3; ++ilayer){
        SurfaceId sid(fpsid[idisk],ilayer);
        auto annulus = std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,toDetectorZ(zlayer[ilayer])),fprin,fprout);
        kkg_->caloMaterialPlanes_.emplace_back(sid,annulus,*mat[ilayer],thick[ilayer]);
        kkg_->map_.emplace(std::make_pair(sid,std::static_pointer_cast<Surface>(annulus)));
      }
      // the pipe centerlines lie one pipe radius downstream of the upstream face of the foam
      double const zpipe = fpback - carbonthick - foamthick + g4info.get<double>("pipeRadius");
      auto const center = det->toDetector(calo_det.disk(idisk).diskInfo().origin());
      makeCaloPipes(g4info,SurfaceId(fpsid[idisk],3),VEC3(center.x(),center.y(),toDetectorZ(zpipe)));
    }
  }

  // Cooling pipes of a calorimeter front panel, reproducing caloBuildFrontPlate in Mu2eG4/src/constructDiskCalorimeter.cc:
  // for each pipe a large torus arc (placed twice, rotated by pi), a small bending torus and a straight section running out
  // to the manifold (each placed 4 times, mirrored in x and y), plus the outer manifold torus. G4 rotations by pi about
  // z, y and x map (x,y) to (-x,-y), (-x,y) and (x,-y) respectively. Pipe interiors are foam, already in the foam layer.
  void KinKalGeomMaker::makeCaloPipes(CaloG4Info const& g4info, SurfaceId const& sid, VEC3 const& center) {
    using Pipes = KKGeom::CaloFrontPanelPipes;
    double const deg = M_PI/180.0;
    int const npipes = g4info.get<int>("nPipes");
    double const prout = g4info.get<double>("pipeRadius");
    double const prin = prout - g4info.get<double>("pipeThickness");
    double const radSmTor = g4info.get<double>("radSmTor");
    double const xsmall = g4info.get<double>("xsmall");
    double const xdistance = g4info.get<double>("xdistance");
    auto const torRadius = g4info.get<std::vector<double>>("pipeTorRadius");
    auto const largeTorPhi = g4info.get<std::vector<double>>("largeTorPhi");
    auto const smallTorPhi = g4info.get<std::vector<double>>("smallTorPhi");
    auto const yposition = g4info.get<std::vector<double>>("yposition");
    auto const straightEndPhi = g4info.get<std::vector<double>>("straightEndPhi");
    double const coolTorRadius = g4info.get<double>("FPCoolPipeTorRadius");
    double const coolrout = g4info.get<double>("FPCoolPipeRadius");
    double const coolrin = coolrout - g4info.get<double>("FPCoolPipeThickness");
    double const fprin = g4info.get<double>("FPInnerRadius");
    double const fprout = g4info.get<double>("FPOuterRadius");
    double const cx = center.X(), cy = center.Y();

    auto plane = std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),center,fprin,coolTorRadius+coolrout);
    // KinKal material for the aluminium (G4_Al) pipe walls, named per subsystem like the CRV strongback
    Pipes pipes(sid,plane,"CaloPipeAl");
    for(int ipipe=0; ipipe < npipes; ++ipipe){
      double const angle = 0.5*largeTorPhi.at(ipipe)*deg;
      double const sAngle = smallTorPhi.at(ipipe)*deg;
      double const sx = xsmall + xdistance*ipipe;
      double const sy = yposition.at(ipipe);
      double const sephi = straightEndPhi.at(ipipe)*deg;
      // large torus
      for(double sign : {1.0,-1.0})
        pipes.addArc(Pipes::Arc{cx,cy,sign,sign,torRadius.at(ipipe),-angle*0.9975,2*angle*0.9975,prout,prin});
      // small bending torus
      for(auto const& [xs,ys] : {std::pair(1.0,1.0),std::pair(-1.0,-1.0),std::pair(-1.0,1.0),std::pair(1.0,-1.0)})
        pipes.addArc(Pipes::Arc{cx+xs*sx,cy+ys*sy,xs,ys,radSmTor,angle-sAngle+M_PI,sAngle,prout,prin});
      // straight section, shortened by 2 mm at each end as in G4
      double const xend = (coolTorRadius - coolrout)*std::sin(sephi);
      double const yend = (coolTorRadius - coolrout)*std::cos(sephi);
      double const xstart = sx - radSmTor*std::cos(sephi);
      double const ystart = sy + radSmTor*std::sin(sephi);
      double const slen = std::sqrt((xend-xstart)*(xend-xstart) + (yend-ystart)*(yend-ystart));
      double const ux = (xend-xstart)/slen, uy = (yend-ystart)/slen;
      for(auto const& [xs,ys] : {std::pair(1.0,1.0),std::pair(-1.0,1.0),std::pair(-1.0,-1.0),std::pair(1.0,-1.0)})
        pipes.addLine(Pipes::Line{cx+xs*(xstart+2.0*ux),cy+ys*(ystart+2.0*uy),cx+xs*(xend-2.0*ux),cy+ys*(yend-2.0*uy),prout,prin});
    }
    // manifold torus around the panel, placed rotated by pi about y
    double const pipeTotalSep = torRadius.back() - torRadius.front();
    double const angMax = 0.5*M_PI - std::asin((g4info.get<double>("pipeInitSeparation")+pipeTotalSep)/fprout) - 0.1;
    pipes.addArc(Pipes::Arc{cx,cy,-1.0,1.0,coolTorRadius,angMax,2*M_PI-2*angMax,coolrout,coolrin});
    kkg_->map_.emplace(std::make_pair(sid,std::static_pointer_cast<Surface>(plane)));
    kkg_->caloFrontPanelPipes_.push_back(std::move(pipes));
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
    auto front= std::make_shared<KinKal::Disk>(outer->frontDisk());
    auto back= std::make_shared<KinKal::Disk>(outer->backDisk());
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
    // Each annulus carries the full wall thickness and spans the gap between the inner and outer shells (constructDS).
    {
      double const dsz = ds->position().z();
      auto addEndWall = [this,&toDetectorZ](SurfaceId const& sid, double zmu2e, double rin,
          double rout, double thick, std::string const& material) {
        auto annulus = std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),
            VEC3(0.0,0.0,toDetectorZ(zmu2e)),rin,rout);
        kkg_->passiveMaterialPlanes_.emplace_back(sid,annulus,material,thick);
        kkg_->map_.emplace(std::make_pair(sid,std::static_pointer_cast<Surface>(annulus)));
      };
      double const cOff = ds->halfLength()        - ds->endWallHalfLength();        // cryo wall center inset
      double const sOff = ds->shield_halfLength() - ds->shield_endWallHalfLength(); // shield wall center inset
      // cryostat end walls (StainlessSteel)
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Front,1), dsz-cOff, ds->rIn2(), ds->rOut1(),
          2.0*ds->endWallHalfLength(), ds->material());
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Back,1),  dsz+cOff, ds->rIn2(), ds->rOut1(),
          2.0*ds->endWallHalfLength(), ds->material());
      // thermal-shield end walls (aluminium; see dsShieldMaterial above)
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Front,2), dsz-sOff, ds->shield_rIn2(), ds->shield_rOut1(),
          2.0*ds->shield_endWallHalfLength(), dsShieldMaterial);
      addEndWall(SurfaceId(SurfaceIdEnum::DS_Back,2),  dsz+sOff, ds->shield_rIn2(), ds->shield_rOut1(),
          2.0*ds->shield_endWallHalfLength(), dsShieldMaterial);
    }

    // hard-coded for now
    auto ipa= std::make_shared<Cylinder>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-2770),300.0,500.0);
    auto ipafront= std::make_shared<KinKal::Disk>(ipa->frontDisk());
    auto ipaback= std::make_shared<KinKal::Disk>(ipa->backDisk());
    auto opa= std::make_shared<Frustrum>(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,-3766),454.0,728.4,2125.0); // inner surface
    auto tsda= std::make_shared<Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-5967),235.0,525.0); // back surface
    // TS5 face: a disk filling the TSDA hole, offset just downstream of the TSDA back surface so an upstream-going track
    // crosses it before the TSDA. A track crossing it leaves the DS into the TS, which is rare for anything but cosmics
    static constexpr double ts5FaceOffset = 1.0; // mm
    auto ts5face= std::make_shared<KinKal::Disk>(tsda->normal(),tsda->uDirection(),tsda->center()+VEC3(0.0,0.0,ts5FaceOffset),tsda->innerRadius());

    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Front),std::static_pointer_cast<Surface>(front)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Back),std::static_pointer_cast<Surface>(back)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Inner),std::static_pointer_cast<Surface>(inner)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::DS_Outer),std::static_pointer_cast<Surface>(outer)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA),std::static_pointer_cast<Surface>(ipa)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA_Front),std::static_pointer_cast<Surface>(ipafront)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::IPA_Back),std::static_pointer_cast<Surface>(ipaback)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::OPA),std::static_pointer_cast<Surface>(opa)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TSDA),std::static_pointer_cast<Surface>(tsda)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::TS5_Face),std::static_pointer_cast<Surface>(ts5face)));

    kkg_->ds_ = std::make_unique<KKGeom::DetectorSolenoid>(inner, outer, front, back, ipa, ipafront, ipaback, opa, tsda, ts5face, materialCylinders);
  }

  void KinKalGeomMaker::makeTarget() {
    // build the foils from the stopping target served by GeometryService, so that any target geometry
    // (foil count, positions, tilts, radii, holes and thicknesses) is used in the fit and extrapolation
    GeomHandle<StoppingTarget> target;
    if(target->nFoils() <= 0) throw cet::exception("GEOM") << "KinKalGeomMaker: the stopping target has no foils" << std::endl;
    std::vector<AnnPtr> foils;
    std::vector<double> thicknesses;
    foils.reserve(target->nFoils());
    thicknesses.reserve(target->nFoils());
    // the bounding surfaces enclose all the foils: outer radius, inner (hole) radius and z extent
    double rout(0.0), rin(std::numeric_limits<double>::max());
    double zmin(std::numeric_limits<double>::max()), zmax(std::numeric_limits<double>::lowest());
    VEC3 xdir(1.0,0.0,0.0);
    for(int ifoil=0; ifoil < target->nFoils(); ++ifoil){
      auto const& tfoil = target->foil(ifoil);
      VEC3 center(tfoil.centerInDetectorSystem());
      VEC3 norm = VEC3(tfoil.normal()).Unit();
      if(std::fabs(norm.X()) > 0.99) throw cet::exception("GEOM") << "KinKalGeomMaker: stopping target foil " << ifoil
        << " normal " << norm << " lies along x" << std::endl;
      // U is the x axis projected onto the foil plane; for an untilted foil it is the x axis
      VEC3 udir = (xdir - xdir.Dot(norm)*norm).Unit();
      foils.push_back(std::make_shared<Annulus>(norm,udir,center,tfoil.rIn(),tfoil.rOut()));
      thicknesses.push_back(2.0*tfoil.halfThickness());
      // a foil off the axis reaches further out, and its hole edge comes closer to the axis
      double roff = center.Rho();
      rout = std::max(rout, tfoil.rOut() + roff);
      rin = std::min(rin, tfoil.rIn() - roff);
      // a tilted foil extends in z beyond its center by its outer radius times the sine of the tilt
      double dz = tfoil.halfThickness() + tfoil.rOut()*std::sqrt(std::max(0.0, 1.0 - norm.Z()*norm.Z()));
      zmin = std::min(zmin, center.Z() - dz);
      zmax = std::max(zmax, center.Z() + dz);
    }
    VEC3 zdir(0.0,0.0,1.0);
    VEC3 zmid(0.0,0.0,0.5*(zmin+zmax));
    double halflen = 0.5*(zmax-zmin);
    auto outer = std::make_shared<Cylinder>(zdir,zmid,rout,halflen);
    // no inner boundary when the hole of some foil reaches the axis (a foil without a hole)
    std::shared_ptr<Cylinder> inner;
    if(rin > 0.0) inner = std::make_shared<Cylinder>(zdir,zmid,rin,halflen);
    auto front= std::make_shared<KinKal::Disk>(outer->frontDisk());
    auto back= std::make_shared<KinKal::Disk>(outer->backDisk());

    if(debug_ > 0){
      std::cout << "KinKalGeomMaker: stopping target with " << foils.size() << " foils, z " << zmin << " to " << zmax
        << " outer radius " << rout << " inner radius ";
      if(inner) std::cout << rin; else std::cout << "none";
      std::cout << std::endl;
      for(size_t ifoil=0; ifoil < foils.size(); ++ifoil){
        auto const& foil = *foils[ifoil];
        std::cout << "  foil " << ifoil << " center " << foil.center() << " normal " << foil.normal()
          << " radii " << foil.innerRadius() << " " << foil.outerRadius() << " thickness " << thicknesses[ifoil] << std::endl;
      }
    }

    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Front),std::static_pointer_cast<Surface>(front)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Back),std::static_pointer_cast<Surface>(back)));
    if(inner) kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Inner),std::static_pointer_cast<Surface>(inner)));
    kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Outer),std::static_pointer_cast<Surface>(outer)));
    for(size_t ifoil=0;ifoil < foils.size();++ifoil){
      kkg_->map_.emplace(std::make_pair(SurfaceId(SurfaceIdEnum::ST_Foils,ifoil),std::static_pointer_cast<Surface>(foils[ifoil])));
    }

    kkg_->st_ = std::make_unique<KKGeom::StoppingTarget>(outer,inner,front,back,foils,thicknesses);
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
      // each layer's bar spans across the bars, one per module, so the material of a crossing can follow the layer stagger
      auto const svdir = sector.sector_->vDirection(); // the axis layerFraction uses
      for(int il = 0; il < nlay; ++il){
        mu2e::KKGeom::KKCRVLayer layer;
        layer.w_ = (VEC3(det->toDetector(firstmod.getLayer(il).getPosition())) - midpoint).Dot(sector.sector_->normal());
        for(int im = 0; im < shield.nModules(); ++im){
          auto const& lay = shield.getModule(im).getLayer(il);
          double const va = (VEC3(det->toDetector(lay.getBar(0).getPosition())) - midpoint).Dot(svdir);
          double const vb = (VEC3(det->toDetector(lay.getBar(lay.nBars()-1).getPosition())) - midpoint).Dot(svdir);
          layer.v_.emplace_back(std::min(va,vb) - firstbar.getHalfWidth(), std::max(va,vb) + firstbar.getHalfWidth());
        }
        sector.layers_.push_back(layer);
      }
      // Strongback: one plane over aluminium sheet 0 (layer 0) of all the sector's modules, with the full sheet thickness
      if(firstmod.nAluminumSheets() > 0){
        int const td = bardet.getThicknessDirection();
        int const ld = bardet.getLengthDirection();
        int const wd = bardet.getWidthDirection();
        CLHEP::Hep3Vector lo(1.0e30,1.0e30,1.0e30), hi(-1.0e30,-1.0e30,-1.0e30);
        for(int im = 0; im < shield.nModules(); ++im){
          auto const& sheet = shield.getModule(im).getAluminumSheet(0);
          for(int k = 0; k < 3; ++k){
            lo[k] = std::min(lo[k], sheet.getPosition()[k] - sheet.getHalfLengths().at(k));
            hi[k] = std::max(hi[k], sheet.getPosition()[k] + sheet.getHalfLengths().at(k));
          }
        }
        std::string const strongBackMaterial = "CRVStrongbackAl";
        auto const sbcenter = VEC3(det->toDetector(0.5*(lo+hi)));
        auto const sbplane = std::make_shared<KinKal::Rectangle>(wdir,udir,sbcenter,0.5*(hi[ld]-lo[ld]),0.5*(hi[wd]-lo[wd]));
        SurfaceId sbsid(SurfaceIdEnum::CRV_StrongBack, static_cast<int>(kkg_->passiveMaterialPlanes_.size()));
        kkg_->passiveMaterialPlanes_.emplace_back(sbsid, sbplane, strongBackMaterial, hi[td]-lo[td]);
        kkg_->map_.emplace(std::make_pair(sbsid, std::static_pointer_cast<Surface>(sbplane)));
      }
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
      kkg_->map_.emplace(std::make_pair(SurfaceId(sector.sname_),std::static_pointer_cast<Surface>(sector.sector_)));
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

  namespace {
    // axis-aligned box in Mu2e coordinates
    struct MuBox {
      CLHEP::Hep3Vector lo{1.0e30,1.0e30,1.0e30}, hi{-1.0e30,-1.0e30,-1.0e30};
      bool empty() const { return lo.x() > hi.x(); }
      void extend(MuBox const& b) { for(int k = 0; k < 3; ++k){ lo[k] = std::min(lo[k],b.lo[k]); hi[k] = std::max(hi[k],b.hi[k]); } }
      CLHEP::Hep3Vector center() const { return 0.5*(lo+hi); }
      double half(int k) const { return 0.5*(hi[k]-lo[k]); }
    };
    // Mu2e box of a block-local box, placed as in constructExternalShielding (local point p goes to center + R^-1 p)
    MuBox placeBox(double u0, double u1, double v0, double v1, double w0, double w1,
        CLHEP::Hep3Vector const& center, CLHEP::HepRotation const& rinv) {
      MuBox box;
      for(int icorner = 0; icorner < 8; ++icorner){
        CLHEP::Hep3Vector const global = center + rinv*CLHEP::Hep3Vector((icorner&1) ? u1 : u0, (icorner&2) ? v1 : v0, (icorner&4) ? w1 : w0);
        for(int k = 0; k < 3; ++k){ box.lo[k] = std::min(box.lo[k],global[k]); box.hi[k] = std::max(box.hi[k],global[k]); }
      }
      return box;
    }
    CLHEP::HepRotation inverseRotation(std::string const& orientation) {
      CLHEP::HepRotation rot(CLHEP::HepRotation::IDENTITY);
      OrientationResolver().getRotationFromOrientation(rot,orientation);
      return rot.inverse();
    }
    // Mu2e cuboids of a block: its U-V outline cut into U strips, extruded along W (no notches, holes or tolerances)
    std::vector<MuBox> blockCuboids(std::vector<std::vector<double>> const& outline, double halfLength,
        CLHEP::Hep3Vector const& center, std::string const& orientation) {
      std::vector<MuBox> boxes;
      size_t const nv = outline.size();
      if(nv < 3) return boxes;
      auto const rinv = inverseRotation(orientation);
      std::vector<double> ucuts;
      for(auto const& uv : outline) ucuts.push_back(uv[0]);
      std::sort(ucuts.begin(),ucuts.end());
      ucuts.erase(std::unique(ucuts.begin(),ucuts.end(),[](double a, double b){ return std::fabs(a-b) < 1.0e-6; }),ucuts.end());
      for(size_t is = 0; is+1 < ucuts.size(); ++is){
        double const um = 0.5*(ucuts[is]+ucuts[is+1]);
        std::vector<double> vcross; // outline crossings of the line U = um
        for(size_t iv = 0; iv < nv; ++iv){
          auto const& a = outline[iv];
          auto const& b = outline[(iv+1)%nv];
          if((a[0]-um)*(b[0]-um) < 0.0) vcross.push_back(a[1] + (b[1]-a[1])*(um-a[0])/(b[0]-a[0]));
        }
        std::sort(vcross.begin(),vcross.end());
        for(size_t ic = 0; ic+1 < vcross.size(); ic += 2)
          boxes.push_back(placeBox(ucuts[is],ucuts[is+1],vcross[ic],vcross[ic+1],-halfLength,halfLength,center,rinv));
      }
      return boxes;
    }
    // Mu2e envelope of a block, grown by its G4 tolerances (du, dv, dw)
    MuBox blockEnvelope(std::vector<std::vector<double>> const& outline, double halfLength, std::vector<double> const& tol,
        CLHEP::Hep3Vector const& center, std::string const& orientation) {
      if(outline.empty()) return MuBox();
      double u0 = 1.0e30, u1 = -1.0e30, v0 = 1.0e30, v1 = -1.0e30;
      for(auto const& uv : outline){ u0 = std::min(u0,uv[0]); u1 = std::max(u1,uv[0]); v0 = std::min(v0,uv[1]); v1 = std::max(v1,uv[1]); }
      double const du = tol.size() > 0 ? tol[0] : 0.0, dv = tol.size() > 1 ? tol[1] : 0.0, dw = tol.size() > 2 ? tol[2] : 0.0;
      return placeBox(u0-du,u1+du,v0-dv,v1+dv,-halfLength-dw,halfLength+dw,center,inverseRotation(orientation));
    }
  }

  void KinKalGeomMaker::makePassiveMaterials() {
    GeomHandle<DetectorSystem> det;
    // External shielding (concrete) geometry comes from the ExtShieldDownstream GeometryService object
    art::ServiceHandle<GeometryService> geomService;
    if(!geomService->hasElement<ExtShieldDownstream>()) return;
    GeomHandle<ExtShieldDownstream> esd;
    if(esd->getNumberOfBoxTypes() <= 0) return;
    // Decode every block (vectors parallel per block; getBoxTypes() gives the 1-based type, lengths are halves)
    auto const& esdTypes    = esd->getBoxTypes();
    auto const& esdOutlines = esd->getOutlines();
    auto const& esdLengths  = esd->getLengths();
    auto const& esdMats     = esd->getMaterialNames();
    auto const& esdCenters  = esd->getCentersOfBoxes();
    auto const& esdOrients  = esd->getOrientations();
    auto const& esdTols     = esd->getTolerances();
    struct Block { int type; std::string const* material; std::vector<MuBox> cuboids; MuBox envelope; };
    std::vector<Block> blocks;
    for(size_t ib = 0; ib < esdTypes.size(); ++ib){
      Block block{esdTypes[ib],&esdMats[ib],blockCuboids(esdOutlines[ib],esdLengths[ib],esdCenters[ib],esdOrients[ib]),
        blockEnvelope(esdOutlines[ib],esdLengths[ib],esdTols[ib],esdCenters[ib],esdOrients[ib])};
      if(!block.envelope.empty()) blocks.push_back(std::move(block));
    }
    auto hasType = [](std::vector<int> const& types, int t){ return std::find(types.begin(),types.end(),t) != types.end(); };
    // Envelope of the blocks of the given types (and material) passing the selection
    auto envelope = [&](std::vector<int> const& types, std::string const& material, std::function<bool(Block const&)> const& select){
      MuBox env;
      for(auto const& block : blocks)
        if(hasType(types,block.type) && *block.material == material && (!select || select(block))) env.extend(block.envelope);
      return env;
    };
    auto materials = [&](std::vector<int> const& types){
      std::vector<std::string> mats;
      for(auto const& block : blocks)
        if(hasType(types,block.type) && std::find(mats.begin(),mats.end(),*block.material) == mats.end()) mats.push_back(*block.material);
      return mats;
    };
    // One slab (a Gauss-plane pair) filling a box, normal along 'normalAxis'; in-plane half widths in addConcretePlane order
    auto addSlab = [&](int normalAxis, MuBox const& box, std::string const& material){
      if(box.empty()) return;
      int const a1 = (normalAxis == 0) ? 1 : 0;
      int const a2 = (normalAxis == 2) ? 1 : 2;
      addConcretePlane(*det, normalAxis, box.center(), box.half(a1), box.half(a2), box.half(normalAxis), material);
    };
    // Each row of interlocking blocks is one slab filling the envelope of its blocks.
    // Roof, y-normal: main roof (types 2, 10, 23 over the DS, 29, 30 over the endcap) and upstream row (types 4, 9).
    std::vector<int> const mainRoof = {2,10,23,29,30};
    std::vector<int> const upstreamRoof = {4,9};
    double roofBottom = 1.0e30;
    for(auto const& mat : materials(mainRoof)){
      auto const env = envelope(mainRoof,mat,nullptr);
      addSlab(1,env,mat);
      if(!env.empty()) roofBottom = std::min(roofBottom,env.lo.y());
    }
    for(auto const& mat : materials(upstreamRoof)) addSlab(1,envelope(upstreamRoof,mat,nullptr),mat);
    // Side walls (type 1), x-normal, one per side of the DS axis (the mean x of the wall blocks)
    std::vector<int> const walls = {1};
    std::vector<MuBox> wallEnvelopes;
    double xsum = 0.0; int nwall = 0;
    for(auto const& block : blocks) if(block.type == 1){ xsum += block.envelope.center().x(); ++nwall; }
    if(nwall > 0){
      double const xsplit = xsum/nwall;
      for(auto const& mat : materials(walls)){
        wallEnvelopes.push_back(envelope(walls,mat,[xsplit](Block const& b){ return b.envelope.center().x() > xsplit; }));
        wallEnvelopes.push_back(envelope(walls,mat,[xsplit](Block const& b){ return b.envelope.center().x() < xsplit; }));
        for(auto it = wallEnvelopes.end()-2; it != wallEnvelopes.end(); ++it) addSlab(0,*it,mat);
      }
    }
    // Downstream endcap (types 27, 28 stacked in y), z-normal: one slab per outline strip, under the main roof and
    // beyond the side walls (whose last blocks are half thick where the side columns start)
    std::vector<int> const endcap = {27,28};
    for(auto const& mat : materials(endcap)){
      std::vector<MuBox> strips;
      for(auto const& block : blocks){
        if(!hasType(endcap,block.type) || *block.material != mat) continue;
        for(auto const& cub : block.cuboids){
          auto same = [&cub](MuBox const& s){ for(int k : {0,2}) if(std::fabs(s.lo[k]-cub.lo[k]) > 1.0 || std::fabs(s.hi[k]-cub.hi[k]) > 1.0) return false; return true; };
          auto it = std::find_if(strips.begin(),strips.end(),same);
          if(it == strips.end()) strips.push_back(cub); else it->extend(cub);
        }
      }
      for(auto& strip : strips){
        strip.hi[1] = std::min(strip.hi[1],roofBottom);
        for(auto const& wall : wallEnvelopes)
          if(!wall.empty() && strip.lo.x() < wall.hi.x() && wall.lo.x() < strip.hi.x() && strip.lo.y() < wall.hi.y() && wall.lo.y() < strip.hi.y())
            strip.lo[2] = std::max(strip.lo[2],wall.hi.z());
        if(strip.hi[1] > strip.lo[1] && strip.hi[2] > strip.lo[2]) addSlab(2,strip,mat);
      }
    }

    // Run-1 geometries zero the ExtShieldDownstream concrete T-blocks but retain the building
    // dsAreaHatchBlock concrete (CONCRETE_MARS) between tracker and CRV-T. Extracted geometries
    // have no ExtShieldDownstream object and returned above. Reaching here means the shield is
    // declared (numberOfBoxTypes>0); the run-1 fallback fires when its roof T-blocks are absent.
    if(std::none_of(blocks.begin(),blocks.end(),[](Block const& b){ return b.type == 2; })) {
      addBuildingHatchConcrete(*det);
    }
  }

  void KinKalGeomMaker::addBuildingHatchConcrete(DetectorSystem const& det) {
    // dsAreaHatchBlock(E): Mu2eHall building solids; the hall rotation (constructHall) maps vertex x to Mu2e z, y to x
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
      CLHEP::Hep3Vector center(off.x() + mid(verts,false), off.y(), off.z() + mid(verts,true));
      addConcretePlane(det, 1, center, spanHalf(verts,false), spanHalf(verts,true),
                       solid.getYhalfThickness(), solid.getMaterial());
    };
    addSlab("dsAreaHatchBlock");
    addSlab("dsAreaHatchBlockE");
  }
}
