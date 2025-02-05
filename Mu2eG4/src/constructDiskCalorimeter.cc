//
// Free function to create the Disk calorimeter
//
// Original author Bertrand Echenard
//
// Note: Virtual detectors are added inside the disk unit
//

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/Mu2eG4/inc/constructDiskCalorimeter.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"
#include "Offline/Mu2eG4/inc/constructDS.hh"

#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Torus.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4UnionSolid.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4VisAttributes.hh"
#include "Geant4/G4Colour.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4UnitsTable.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4String.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/TwoVector.h"
#include "boost/regex.hpp"



namespace mu2e {



  VolumeInfo constructDiskCalorimeter(const VolumeInfo& mother, const SimpleConfig& config)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry( config, "calorimeterEnvelope", "calorimeter.envelope");
    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper            = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry& reg           = (art::ServiceHandle<Mu2eG4Helper>())->antiLeakRegistry();
    const DiskCalorimeter& cal      = *(GeomHandle<DiskCalorimeter>());
    const DetectorSolenoid& ds      = *(GeomHandle<DetectorSolenoid>());

    G4Material* vacuumMaterial      = materialFinder.get("calorimeter.vacuumMaterial");

    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("calorimeterEnvelope");
    const int  verbosity            = config.getInt("calorimeter.verbosityLevel",1);

    const unsigned nDisks           = cal.nDisks();
    const double caloDiskRadiusIn   = cal.caloInfo().getDouble("caloDiskRadiusIn");
    const double caloDiskRadiusOut  = cal.caloInfo().getDouble("caloDiskRadiusOut");
    const double caloFEBRadiusOut   = cal.caloInfo().getDouble("caloFEBRadiusOut");
    const double mother_z0          = cal.caloInfo().getDouble("caloMotherZ0");
    const double mother_z1          = cal.caloInfo().getDouble("caloMotherZ1");
    const double mother_zlength     = mother_z1-mother_z0;
    const double mother_zCenter     = (mother_z1+mother_z0)/2.0;
    const double FEBOffset          = cal.caloInfo().getDouble("FEBToDiskZOffset");
    const auto   FEBPhiMinMax       = calcFEBPhiRange(cal);
    const bool   hasCrates          = cal.caloInfo().getBool("hasCrates");
    const bool   hasCable           = ds.hasCableRunCal() && hasCrates;



    //--------------------------------------
    // Construct calorimeter disks mother volume
    //
    caloVolInfG4.clear();

    const auto& posDS3InWorld        = mother.centerInWorld;
    const auto& posDS3               = mother.centerInMu2e();
    G4ThreeVector posDiskMother      = G4ThreeVector(posDS3.x(), 0, mother_zCenter);
    G4ThreeVector posDiskMotherInDS  = posDiskMother - posDS3;


    // Disk and FEB mother volumes
    VolumeInfo caloMotherInfo("CalorimeterMother");
    VolumeInfo caloDiskInfo("CaloDiskMother");
    VolumeInfo caloFEBInfo("CaloFEBMother");

    caloDiskInfo.solid           = new G4Tubs(caloDiskInfo.name, caloDiskRadiusIn,  caloDiskRadiusOut, mother_zlength/2.0, 0.0, CLHEP::twopi);
    caloFEBInfo.solid            = new G4Tubs(caloFEBInfo.name,  caloDiskRadiusOut, caloFEBRadiusOut,  mother_zlength/2.0, FEBPhiMinMax[0], FEBPhiMinMax[1]-FEBPhiMinMax[0]);
    caloMotherInfo.solid         = new G4UnionSolid(caloMotherInfo.name, caloDiskInfo.solid, caloFEBInfo.solid);
    caloMotherInfo.logical       = caloLogical(caloMotherInfo, vacuumMaterial, 0, G4Color::Black(), 0, 0);
    caloMotherInfo.physical      = caloPlacement(caloMotherInfo, mother, 0, posDiskMotherInDS, false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(caloFEBInfo);
    helper.addVolInfo(caloDiskInfo);
    helper.addVolInfo(caloMotherInfo);

    if ( verbosity > 0) {
      double zhl = dynamic_cast<G4Tubs*>(caloMotherInfo.solid)->GetZHalfLength();
      double calorimeterOffsetInMu2eZ = caloMotherInfo.centerInMu2e()[CLHEP::Hep3Vector::Z];
      G4cout << __func__ << " Calorimeter mother center in Mu2e   : " << caloMotherInfo.centerInMu2e() << G4endl;
      G4cout << __func__ << " Calorimeter mother Z extent in Mu2e : " << calorimeterOffsetInMu2eZ - zhl << ", " << calorimeterOffsetInMu2eZ + zhl << G4endl;
    }


    // Disk and FEB volumes
    std::vector<VolumeInfo> calorimeterDisk(nDisks), calorimeterFEB(nDisks);
    for (size_t idisk=0;idisk<nDisks;++idisk) {
      G4ThreeVector posDisk           = cal.disk(idisk).geomInfo().origin() - posDiskMother;
      calorimeterDisk[idisk]          = caloBuildDisk(config,idisk);

      G4RotationMatrix* rot           = reg.add(new G4RotationMatrix(cal.disk(idisk).geomInfo().rotation()));
      //auto rot                        = const_cast<G4RotationMatrix*>(&cal.disk(idisk).geomInfo().rotation());
      calorimeterDisk[idisk].physical = caloPlacement(calorimeterDisk[idisk], caloMotherInfo, rot, posDisk, false, 0, config, doSurfaceCheck, verbosity);
      helper.addVolInfo(calorimeterDisk[idisk]);

      if (hasCrates) {
        calorimeterFEB[idisk]          = caloBuildFEB(config,idisk);
        G4Tubs*  FEB                   = dynamic_cast<G4Tubs*>(calorimeterFEB[idisk].solid);
        G4Tubs*  disk                  = dynamic_cast<G4Tubs*>(calorimeterDisk[idisk].solid);
        G4ThreeVector posFEB           = posDisk + G4ThreeVector(0,0,FEB->GetZHalfLength() - disk->GetZHalfLength() - FEBOffset);
        calorimeterFEB[idisk].physical = caloPlacement(calorimeterFEB[idisk], caloMotherInfo, rot, posFEB, false, 0, config, doSurfaceCheck, verbosity);
        helper.addVolInfo(calorimeterFEB[idisk]);

        if (hasCable) {
          VolumeInfo CableRun    = caloBuildCable(config, idisk, calorimeterFEB[idisk]);
          G4Tubs*    Cable       = dynamic_cast<G4Tubs*>(CableRun.solid);
          G4ThreeVector posCable = posFEB + G4ThreeVector(0,0, FEB->GetZHalfLength() + Cable->GetZHalfLength());
          CableRun.physical      = caloPlacement(CableRun, caloMotherInfo, rot, posCable, false, 0, config, doSurfaceCheck, verbosity);
          helper.addVolInfo(CableRun);
        }
      }

      if (verbosity > 0) {
        G4Tubs* disk      = dynamic_cast<G4Tubs*>(calorimeterDisk[idisk].solid);
        G4double zHalftot = disk->GetZHalfLength();
        G4cout << __func__ << " CalorimeterDisk center in Mu2e    : " << calorimeterDisk[idisk].centerInMu2e() << G4endl;
        G4cout << __func__ << " CalorimeterDisk Z extent in Mu2e  : " << calorimeterDisk[idisk].centerInMu2e()[CLHEP::Hep3Vector::Z] - zHalftot
                           << ", " << caloDiskInfo.centerInMu2e()[CLHEP::Hep3Vector::Z] + zHalftot << G4endl;
      }
    }


    //Fix the "world posiiton" of all wolumes
    for (auto& kv : caloVolInfG4) {
      CLHEP::Hep3Vector posWorld(0,0,0);
      std::string mother(kv.first);
      while (caloVolInfG4.find(mother) != caloVolInfG4.end()) {
        posWorld += helper.locateVolInfo(mother).centerInParent;
        mother = caloVolInfG4.find(mother)->second;
      }
      helper.locateVolInfo(kv.first).centerInWorld = posWorld + posDS3InWorld;
    }

    if (verbosity >1) {
      if (verbosity==99) {
        std::cout<<"List of all volumes in the calorimeter and FEB"<<std::endl;
        browseCaloSolids(caloDiskInfo.logical);
        browseCaloSolids(caloFEBInfo.logical);
      }

      for (auto& kv : caloVolInfG4){
         const auto& vol = helper.locateVolInfo(kv.first);
         std::cout<<vol.name<<"   posInWorld="<<vol.centerInWorld<<"    posInMu2e="<<vol.centerInMu2e()<<std::endl;
      }
    }

    if (verbosity > 0) G4cout << __func__ << " Calorimeter constructed "<< G4endl;

    return caloDiskInfo;
  }



  //--------------------------------------------------------------------------------------------------------------------------------
  // Full disk with front plate, crystal casing and backplate
  //
  VolumeInfo caloBuildDisk(const SimpleConfig& config, unsigned idisk)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry(config, "calorimeterEnvelope", "calorimeter.envelope");
    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper       = *(art::ServiceHandle<Mu2eG4Helper>());
    const DiskCalorimeter& cal = *(GeomHandle<DiskCalorimeter>());

    const bool doSurfaceCheck  = geomOptions->doSurfaceCheck("calorimeterEnvelope");
    const int  verbosity       = config.getInt("calorimeter.verbosityLevel",1);

    const bool hasFrontPanel   = cal.caloInfo().getBool("hasFrontPanel");
    const bool hasBackPanel    = cal.caloInfo().getBool("hasBackPanel");
    const double vdThickness   = cal.caloInfo().getDouble("vdThickness");
    const double R0disk        = cal.caloInfo().getDouble("caloDiskRadiusIn");
    const double R1disk        = cal.caloInfo().getDouble("caloDiskRadiusOut");

    G4Material* vacuumMaterial = materialFinder.get("calorimeter.vacuumMaterial");

    VolumeInfo crystalCase     = caloBuildCase(config,idisk);
    VolumeInfo frontPlate      = caloBuildFrontPlate(config,idisk);
    VolumeInfo backPlate       = caloBuildBackPlate(config,idisk);

    G4double zHalfDisk         = dynamic_cast<G4Tubs*>(crystalCase.solid)->GetZHalfLength();
    G4double zHalfFP           = (hasFrontPanel) ? dynamic_cast<G4Tubs*>(frontPlate.solid)->GetZHalfLength() : 0;
    G4double zHalfBP           = (hasBackPanel)  ? dynamic_cast<G4Tubs*>(backPlate.solid)->GetZHalfLength()  : 0;
    G4double zHalftot          = zHalfFP + zHalfDisk + zHalfBP + vdThickness;


    VolumeInfo fullDisk("CaloDisk_"+std::to_string(idisk));
    fullDisk.solid   = new G4Tubs(fullDisk.name,R0disk,R1disk,zHalftot,0,CLHEP::twopi);
    fullDisk.logical = caloLogical(fullDisk, vacuumMaterial, 0, G4Color::Black(), 0, 0);

    auto crystalCasePos  = G4ThreeVector(0,0,-zHalftot+vdThickness+2.0*zHalfFP+zHalfDisk);
    crystalCase.physical = caloPlacement(crystalCase, fullDisk, 0, crystalCasePos, false, 0, config, doSurfaceCheck, verbosity);
    helper.addVolInfo(crystalCase);

    if (hasFrontPanel){
      auto frontPlatePos  = G4ThreeVector(0,0,-zHalftot+vdThickness+zHalfFP);
      frontPlate.physical = caloPlacement(frontPlate, fullDisk, 0, frontPlatePos, false, 0, config, doSurfaceCheck, verbosity);
      helper.addVolInfo(frontPlate);
    }

    if (hasBackPanel) {
      auto backPlatePos  = G4ThreeVector(0,0,+zHalftot-vdThickness-zHalfBP);
      backPlate.physical = caloPlacement(backPlate, fullDisk, 0, backPlatePos, false, 0, config, doSurfaceCheck, verbosity);
      helper.addVolInfo(backPlate);
    }


    //------------------------------------------------------------
    // Perform a few cross-checks to make sure the geometry description matches with DiskCAlorimeterMaker
    //
    std::vector<const G4LogicalVolume*> volumeNodes;
    auto volumePtr = findCaloSolid(fullDisk.logical,"CaloCrystalCsI_"+std::to_string(idisk), volumeNodes);
    volumeNodes.push_back(fullDisk.logical);

    //navigate the nodes to calculate the total translation, and add half crystal size to match the originToCrystalOrigin value in geomInfo
    double  zTranslation(0);
    for (size_t i=1;i<volumeNodes.size();++i) {
      for (size_t j=0;j<volumeNodes[i]->GetNoDaughters();++j) {
        if (volumeNodes[i]->GetDaughter(j)->GetLogicalVolume() == volumeNodes[i-1]) {
          zTranslation += volumeNodes[i]->GetDaughter(j)->GetTranslation().z();
          break;
        }
      }
    }
    if (volumePtr) zTranslation -= dynamic_cast<G4Box*>(volumePtr->GetSolid())->GetZHalfLength();

    double delta1 = 2*zHalftot-cal.disk(idisk).geomInfo().size().z();
    double delta2 = zTranslation-cal.disk(idisk).geomInfo().originToCrystalOrigin().z();

    if (std::abs(delta1) > 1e-3)  G4cout << __func__ <<"PANIC..... geometry description in Geant4 and DiskMaker do NOT match  - disk size: "
                                         << 2*zHalftot<<" vs "<<cal.disk(idisk).geomInfo().size().z()<<G4endl;

    if (std::abs(delta2) > 1e-3)  G4cout << __func__ <<"PANIC..... geometry description in Geant4 and DiskMaker do NOT match  - originToCrystalOrigin: "
                                         << zTranslation<<" vs "<<cal.disk(idisk).geomInfo().originToCrystalOrigin().z()<<G4endl;

    if (verbosity){
      G4cout << __func__ <<" Compare disk size             Geant4 / CaloInfo "<<2*zHalftot<<" / "<<cal.disk(idisk).geomInfo().size().z()<<G4endl;
      G4cout << __func__ <<" Compare originToCrystalOrigin Geant4 / CaloInfo "<<zTranslation<<" / "<<cal.disk(idisk).geomInfo().originToCrystalOrigin().z()<<G4endl;
    }

    return fullDisk;
  }


  //--------------------------------------------------------------------------------------------------------------------------------
  // Front plate: a sandwich of carbon fiberplate - foam plane - carbon fiber. The pipes are inside the foam layer, closest to the beam.
  // the big pipe is aligned with the small pipes (centered on the same z). Big cooling pipe needs to be larger than
  // carbon thickness+small pipe radius and smaller than Carbon thick + foam thick - small pipe radius or the model needs to be updated.
  //
  VolumeInfo caloBuildFrontPlate(const SimpleConfig& config, unsigned idisk)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry(config, "calorimeterEnvelope", "calorimeter.envelope");

    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper                     = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry& reg                    = helper.antiLeakRegistry();
    const DiskCalorimeter& cal               = *(GeomHandle<DiskCalorimeter>());

    const bool isPipeVisible                 = geomOptions->isVisible          ("calorimeterPipe");
    const bool isPipeSolid                   = geomOptions->isSolid            ("calorimeterPipe");
    const bool forceEdge                     = geomOptions->forceAuxEdgeVisible("calorimeterPipe");
    const bool doSurfaceCheck                = geomOptions->doSurfaceCheck     ("calorimeterPipe");
    const int  verbosity                     = config.getInt                   ("calorimeter.verbosityLevel",1);

    const double FPInnerRadius               = cal.caloInfo().getDouble("FPInnerRadius");
    const double FPOuterRadius               = cal.caloInfo().getDouble("FPOuterRadius");
    const double FPCarbonDZ                  = cal.caloInfo().getDouble("FPCarbonZLength")/2.0;
    const double FPFoamDZ                    = cal.caloInfo().getDouble("FPFoamZLength")/2.0;
    const double FPCoolPipeTorRadius         = cal.caloInfo().getDouble("FPCoolPipeTorRadius");
    const double FPCoolPipeRadius            = cal.caloInfo().getDouble("FPCoolPipeRadius");
    const double FPCoolPipeThickness         = cal.caloInfo().getDouble("FPCoolPipeThickness");
    const double FPCoolPipeRadiusIn          = FPCoolPipeRadius-FPCoolPipeThickness;

    const int nPipes                         = cal.caloInfo().getInt    ("nPipes");
    const double pipeRadius                  = cal.caloInfo().getDouble ("pipeRadius");
    const double pipeThickness               = cal.caloInfo().getDouble ("pipeThickness");
    const double pipeInitSeparation          = cal.caloInfo().getDouble ("pipeInitSeparation");
    const std::vector<double> pipeTorRadius  = cal.caloInfo().getVDouble("pipeTorRadius");
    const std::vector<double> largeTorPhi    = cal.caloInfo().getVDouble("largeTorPhi");
    const std::vector<double> smallTorPhi    = cal.caloInfo().getVDouble("smallTorPhi");
    const std::vector<double> yposition      = cal.caloInfo().getVDouble("yposition");
    const std::vector<double> straightEndPhi = cal.caloInfo().getVDouble("straightEndPhi");
    const double radSmTor                    = cal.caloInfo().getDouble ("radSmTor");
    const double xsmall                      = cal.caloInfo().getDouble ("xsmall");
    const double xdistance                   = cal.caloInfo().getDouble ("xdistance");

    const double frontPanelHalfThick         = (2.0*FPCarbonDZ+2.0*FPFoamDZ-pipeRadius+FPCoolPipeRadius)/2.0;
    const double ZposCarbon2                 = frontPanelHalfThick-FPCarbonDZ;
    const double ZposFoam                    = ZposCarbon2-FPCarbonDZ-FPFoamDZ;
    const double ZposCarbon1                 = ZposFoam-FPFoamDZ-FPCarbonDZ;
    const double ZposPipe                    = frontPanelHalfThick-2*FPCarbonDZ-2*FPFoamDZ+pipeRadius;
    const double pipeTotalSep                = pipeTorRadius.back()-pipeTorRadius.front();
    const double angMax                      = CLHEP::pi/2.0-std::asin((pipeInitSeparation+pipeTotalSep)/FPOuterRadius)-0.1;

    G4Material* vacuumMaterial               = materialFinder.get("calorimeter.vacuumMaterial");
    G4Material* FPFoamMaterial               = materialFinder.get("calorimeter.FPFoamMaterial");
    G4Material* FPCarbonMaterial             = materialFinder.get("calorimeter.FPCarbonMaterial");
    G4Material* pipeMaterial                 = materialFinder.get("calorimeter.pipeMaterial");


    G4RotationMatrix* rotFPPipe = reg.add(new G4RotationMatrix());
    rotFPPipe->rotateY(CLHEP::pi);

    G4RotationMatrix* rotPipe1 = nullptr; //using the identity matrix
    G4RotationMatrix* rotPipe2 = reg.add(new G4RotationMatrix());
    rotPipe2->rotateZ(1.0 * CLHEP::pi);
    G4RotationMatrix* rotPipe3 = reg.add(new G4RotationMatrix());
    rotPipe3->rotateY(1.0 * CLHEP::pi);
    G4RotationMatrix* rotPipe4 = reg.add(new G4RotationMatrix());
    rotPipe4->rotateX(1.0 * CLHEP::pi);

    // front plate volume with carbon fiber, foam and cooling pipes
    VolumeInfo frontPlate    ("CaloFP_"+std::to_string(idisk));
    VolumeInfo frontCFPanel  ("CaloFPfrontCF_"+std::to_string(idisk));
    VolumeInfo frontPanelFoam("CaloFPfoam_"+std::to_string(idisk));
    VolumeInfo backCFPanel   ("CaloFPbackCF_"+std::to_string(idisk));
    VolumeInfo coolFP        ("CaloFPcool_"+std::to_string(idisk));

    frontPlate.solid        = new G4Tubs(frontPlate.name,     FPInnerRadius,FPCoolPipeTorRadius+FPCoolPipeRadius,frontPanelHalfThick,0,CLHEP::twopi);
    frontCFPanel.solid      = new G4Tubs(frontCFPanel.name,   FPInnerRadius, FPOuterRadius, FPCarbonDZ, 0, CLHEP::twopi);
    frontPanelFoam.solid    = new G4Tubs(frontPanelFoam.name, FPInnerRadius, FPOuterRadius, FPFoamDZ,   0, CLHEP::twopi);
    backCFPanel.solid       = new G4Tubs(backCFPanel.name,    FPInnerRadius, FPOuterRadius, FPCarbonDZ, 0, CLHEP::twopi);
    coolFP.solid            = new G4Torus(coolFP.name,        FPCoolPipeRadiusIn, FPCoolPipeRadius, FPCoolPipeTorRadius, angMax, CLHEP::twopi-2.0*angMax);

    frontPlate.logical      = caloLogical(frontPlate,     vacuumMaterial,               0, G4Color::Black(), 0, 0);
    frontCFPanel.logical    = caloLogical(frontCFPanel,   FPCarbonMaterial, isPipeVisible, G4Color::Gray(),  isPipeSolid, forceEdge);
    frontPanelFoam.logical  = caloLogical(frontPanelFoam, FPFoamMaterial,   isPipeVisible, G4Color::Brown(), isPipeSolid, forceEdge);
    backCFPanel.logical     = caloLogical(backCFPanel,    FPCarbonMaterial, isPipeVisible, G4Color::Gray(),  isPipeSolid, forceEdge);
    coolFP.logical          = caloLogical(coolFP,         pipeMaterial,     isPipeVisible, G4Color::Red(),   isPipeSolid, forceEdge);

    frontCFPanel.physical   = caloPlacement(frontCFPanel,  frontPlate, 0, G4ThreeVector(0,0,ZposCarbon1), false, 0, config, doSurfaceCheck, verbosity);
    frontPanelFoam.physical = caloPlacement(frontPanelFoam,frontPlate, 0, G4ThreeVector(0,0,ZposFoam),    false, 0, config, doSurfaceCheck, verbosity);
    backCFPanel.physical    = caloPlacement(backCFPanel,   frontPlate, 0, G4ThreeVector(0,0,ZposCarbon2), false, 0, config, doSurfaceCheck, verbosity);
    coolFP.physical         = caloPlacement(coolFP,        frontPlate, 0, G4ThreeVector(0,0,ZposPipe),    false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(frontCFPanel);
    helper.addVolInfo(frontPanelFoam);
    helper.addVolInfo(backCFPanel);
    helper.addVolInfo(coolFP);


    // pipes inside foam
    for (int ipipe=0; ipipe<nPipes; ++ipipe) {
      double angle  = largeTorPhi[ipipe] * CLHEP::degree / 2;
      double sAngle = smallTorPhi[ipipe] * CLHEP::degree;
      double sxPos  = xsmall + xdistance * ipipe;
      double syPos  = yposition[ipipe];
      double z      = -FPFoamDZ+pipeRadius;

      // calculate the parameters of the straight pipes
      // minus 2.0 mm from the manifold pipe's most inner radius to avoid the overlap.
      double xEnd         = (FPCoolPipeTorRadius -FPCoolPipeRadius) * sin(straightEndPhi[ipipe] * CLHEP::degree);
      double yEnd         = (FPCoolPipeTorRadius -FPCoolPipeRadius) * cos(straightEndPhi[ipipe] * CLHEP::degree);
      double xStart       = sxPos - radSmTor * cos(straightEndPhi[ipipe] * CLHEP::degree);
      double yStart       = syPos + radSmTor * sin(straightEndPhi[ipipe] * CLHEP::degree);
      double sLength      = sqrt((xEnd - xStart)*(xEnd - xStart) + (yEnd - yStart)*(yEnd - yStart));
      double zRotateAngle = std::atan((yEnd - yStart) / (xEnd - xStart));
      double xCenter      = 0.5 * (xEnd + xStart);
      double yCenter      = 0.5 * (yEnd + yStart);

      // rotation coordinate
      G4RotationMatrix* sPipeRotate1 = reg.add(new G4RotationMatrix());
      sPipeRotate1 -> rotateX(0.5 * CLHEP::pi);
      sPipeRotate1 -> rotateY((-0.5*CLHEP::pi + zRotateAngle));
      G4RotationMatrix* sPipeRotate2 = reg.add(new G4RotationMatrix());
      sPipeRotate2 -> rotateX(0.5 * CLHEP::pi);
      sPipeRotate2 -> rotateY((-0.5*CLHEP::pi - zRotateAngle));
      G4RotationMatrix* sPipeRotate3 = reg.add(new G4RotationMatrix());
      sPipeRotate3 -> rotateX(0.5 * CLHEP::pi);
      sPipeRotate3 -> rotateY((-1.5*CLHEP::pi + zRotateAngle));
      G4RotationMatrix* sPipeRotate4 = reg.add(new G4RotationMatrix());
      sPipeRotate4 -> rotateX(0.5 * CLHEP::pi);
      sPipeRotate4 -> rotateY((-1.5*CLHEP::pi - zRotateAngle));

      VolumeInfo pipe1("CaloFP_pipe1_"+std::to_string(ipipe+10*idisk));
      VolumeInfo pipe2("CaloFP_pipe2_"+std::to_string(ipipe+10*idisk));
      VolumeInfo pipe3("CaloFP_pipe3_"+std::to_string(ipipe+10*idisk));

      pipe1.solid    = new G4Torus(pipe1.name, pipeRadius-pipeThickness, pipeRadius, pipeTorRadius[ipipe], -angle*0.9975, 2*angle*0.9975);
      pipe2.solid    = new G4Torus(pipe2.name, pipeRadius-pipeThickness, pipeRadius, radSmTor, angle-sAngle+CLHEP::pi, sAngle);
      pipe3.solid    = new G4Tubs (pipe3.name, pipeRadius-pipeThickness, pipeRadius, 0.5 * (sLength - 4.), 0, CLHEP::twopi);

      pipe1.logical  = caloLogical(pipe1, pipeMaterial, isPipeVisible, G4Color::Red(), isPipeSolid, forceEdge);
      pipe2.logical  = caloLogical(pipe2, pipeMaterial, isPipeVisible, G4Color::Red(), isPipeSolid, forceEdge);
      pipe3.logical  = caloLogical(pipe3, pipeMaterial, isPipeVisible, G4Color::Red(), isPipeSolid, forceEdge);

      // large bending torus
      caloPlacement(pipe1, frontPanelFoam, rotPipe1, G4ThreeVector(0,0,z), true, 0, config, doSurfaceCheck, verbosity);
      caloPlacement(pipe1, frontPanelFoam, rotPipe2, G4ThreeVector(0,0,z), true, 1, config, doSurfaceCheck, verbosity);

      // small bending torus
      caloPlacement(pipe2, frontPanelFoam, rotPipe1, G4ThreeVector(sxPos,  syPos, z), true, 0, config, doSurfaceCheck, verbosity);
      caloPlacement(pipe2, frontPanelFoam, rotPipe2, G4ThreeVector(-sxPos,-syPos, z), true, 1, config, doSurfaceCheck, verbosity);
      caloPlacement(pipe2, frontPanelFoam, rotPipe3, G4ThreeVector(-sxPos, syPos, z), true, 2, config, doSurfaceCheck, verbosity);
      caloPlacement(pipe2, frontPanelFoam, rotPipe4, G4ThreeVector(sxPos, -syPos, z), true, 3, config, doSurfaceCheck, verbosity);

      // straight pipes
      caloPlacement(pipe3, frontPanelFoam, sPipeRotate1, G4ThreeVector(xCenter,  yCenter, z), true, 0, config, doSurfaceCheck, verbosity);
      caloPlacement(pipe3, frontPanelFoam, sPipeRotate2, G4ThreeVector(-xCenter, yCenter, z), true, 1, config, doSurfaceCheck, verbosity);
      caloPlacement(pipe3, frontPanelFoam, sPipeRotate3, G4ThreeVector(-xCenter,-yCenter, z), true, 2, config, doSurfaceCheck, verbosity);
      caloPlacement(pipe3, frontPanelFoam, sPipeRotate4, G4ThreeVector(xCenter, -yCenter, z), true, 3, config, doSurfaceCheck, verbosity);

      helper.addVolInfo(pipe1);
      helper.addVolInfo(pipe2);
      helper.addVolInfo(pipe3);
    }

    return frontPlate;
  }



  //--------------------------------------------------------------------------------------------------------------------------------
  //construct central part of the disk with crystals
  VolumeInfo caloBuildCase(const SimpleConfig& config, unsigned idisk)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry(config, "calorimeterCase",    "calorimeter.case");
    geomOptions->loadEntry(config, "calorimeterCrystal", "calorimeter.crystal");

    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper                  = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry& reg                 = helper.antiLeakRegistry();
    const DiskCalorimeter& cal            = *(GeomHandle<DiskCalorimeter>());

    const bool isDiskVisible              = geomOptions->isVisible          ("calorimeterCase");
    const bool isDiskSolid                = geomOptions->isSolid            ("calorimeterCase");
    const bool isCrystalVisible           = geomOptions->isVisible          ("calorimeterCrystal");
    const bool isCrystalSolid             = geomOptions->isSolid            ("calorimeterCrystal");
    const bool forceEdge                  = geomOptions->forceAuxEdgeVisible("calorimeterCase");
    const bool doSurfaceCheck             = geomOptions->doSurfaceCheck     ("calorimeterCase");
    const int  verbosity                  = config.getInt                   ("calorimeter.verbosityLevel",1);

    G4Material* vacuumMaterial            = materialFinder.get("calorimeter.vacuumMaterial");
    G4Material* shimMaterial              = materialFinder.get("calorimeter.vacuumMaterial");
    G4Material* crystalMaterial           = materialFinder.get("calorimeter.crystalMaterial");
    G4Material* crystalCapMaterial        = materialFinder.get("calorimeter.crystalCapMaterial");
    G4Material* wrapMaterial              = materialFinder.get("calorimeter.wrapperMaterial");
    G4Material* caphriCrystalMaterial     = materialFinder.get("calorimeter.caphriCrystalMaterial");
    G4Material* innerAlRingMaterial       = materialFinder.get("calorimeter.innerAlRingMaterial");
    G4Material* innerCFRingMaterial       = materialFinder.get("calorimeter.innerCFRingMaterial");
    G4Material* innerStepMaterial         = materialFinder.get("calorimeter.innerStepMaterial");
    G4Material* outerRingMaterial         = materialFinder.get("calorimeter.outerRingMaterial");
    G4Material* coolPipeMaterial          = materialFinder.get("calorimeter.coolPipeMaterial");

    const double vdThickness              = cal.caloInfo().getDouble("vdThickness");
    const double crystalDXY               = cal.caloInfo().getDouble("crystalXYLength")/2.0;
    const double crystalDZ                = cal.caloInfo().getDouble("crystalZLength")/2.0;
    const double crystalCapDZ             = cal.caloInfo().getDouble("crystalCapZLength")/2.0;
    const double wrapperHalfThick         = cal.caloInfo().getDouble("wrapperThickness")/2.0;
    const double wrapperDXY               = crystalDXY + 2*wrapperHalfThick;
    const double wrapperDZ                = crystalDZ+crystalCapDZ;
    const std::vector<int> caphriCystalId = cal.caloInfo().getVInt("caphriCrystalId");

    const double diskCaseDZ               = wrapperDZ;
    const double diskInAlRingRIn          = cal.caloInfo().getDouble("diskInAlRingRIn");
    const double diskInAlRingDZ           = cal.caloInfo().getDouble("diskInAlRingZLength")/2.0;
    const double diskInCFRingRIn          = cal.caloInfo().getDouble("diskInCFRingRIn");
    const double diskInCFRingROut         = cal.caloInfo().getDouble("diskInCFRingROut");
    const double diskCaseRingROut         = cal.caloInfo().getDouble("diskCaseRingROut");
    const double diskOutRailROut          = cal.caloInfo().getDouble("diskOutRailROut") - vdThickness;
    const double diskOutRailDZ            = cal.caloInfo().getDouble("diskOutRailZLength")/2.0;

    const double FPCoolPipeRadius         = cal.caloInfo().getDouble("FPCoolPipeRadius");
    const double FPCoolPipeThickness      = cal.caloInfo().getDouble("FPCoolPipeThickness");
    const double coolPipeZpos             = (diskCaseDZ - 2*FPCoolPipeRadius - 2.0*diskOutRailDZ)/2.0 + FPCoolPipeRadius;

    const std::vector<double> stepsInX    = cal.caloInfo().getVDouble("stepsInsideX");
    const std::vector<double> stepsInY    = cal.caloInfo().getVDouble("stepsInsideY");
    const std::vector<double> stepsOutX   = cal.caloInfo().getVDouble("stepsOutsideX");
    const std::vector<double> stepsOutY   = cal.caloInfo().getVDouble("stepsOutsideY");
    const double diskStepThickness        = cal.caloInfo().getDouble ("diskStepThickness");


    //------------------------------------------------------------
    // Build disk casing, inner carbon fiber ring, inner Al rings and outer rails
    // Crystals, inner and outer steps are placed into the casing.
    // The case ring is made out of shim material to fill holes between crystals and walls
    //

    VolumeInfo fullCase    ("CaloFullCase_"+std::to_string(idisk));
    VolumeInfo innerAlRingF("CaloInnerAlRingF_"+std::to_string(idisk));
    VolumeInfo innerAlRingM("CaloInnerAlRingM_"+std::to_string(idisk));
    VolumeInfo innerAlRingB("CaloInnerAlRingB_"+std::to_string(idisk));
    VolumeInfo innerCFRing ("CaloInnerCFRing_"+std::to_string(idisk));
    VolumeInfo caseRing    ("CaloCaseRing_"+std::to_string(idisk));
    VolumeInfo outerRailF  ("CaloOuterRailF_"+std::to_string(idisk));
    VolumeInfo outerRailB  ("CaloOuterRailB_"+std::to_string(idisk));

    fullCase.solid        = new G4Tubs(fullCase.name,     diskInAlRingRIn,  diskOutRailROut,  diskCaseDZ,     0, CLHEP::twopi);
    innerAlRingF.solid    = new G4Tubs(innerAlRingF.name, diskInAlRingRIn,  diskInCFRingRIn,  diskInAlRingDZ, 0, CLHEP::twopi);
    innerAlRingM.solid    = new G4Tubs(innerAlRingM.name, diskInAlRingRIn,  diskInCFRingRIn,  diskInAlRingDZ, 0, CLHEP::twopi);
    innerAlRingB.solid    = new G4Tubs(innerAlRingB.name, diskInAlRingRIn,  diskInCFRingRIn,  diskInAlRingDZ, 0, CLHEP::twopi);
    innerCFRing.solid     = new G4Tubs(innerCFRing.name,  diskInCFRingRIn,  diskInCFRingROut, diskCaseDZ,     0, CLHEP::twopi);
    caseRing.solid        = new G4Tubs(caseRing.name,     diskInCFRingROut, diskCaseRingROut, diskCaseDZ,     0, CLHEP::twopi);
    outerRailF.solid      = new G4Tubs(outerRailF.name,   diskCaseRingROut, diskOutRailROut,  diskOutRailDZ,  0, CLHEP::twopi);
    outerRailB.solid      = new G4Tubs(outerRailB.name,   diskCaseRingROut, diskOutRailROut,  diskOutRailDZ,  0, CLHEP::twopi);

    fullCase.logical      = caloLogical(fullCase,     vacuumMaterial,     0, G4Color::Black(), 0, 0);
    innerAlRingF.logical  = caloLogical(innerAlRingF, innerAlRingMaterial,isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);
    innerAlRingM.logical  = caloLogical(innerAlRingM, innerAlRingMaterial,isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);
    innerAlRingB.logical  = caloLogical(innerAlRingB, innerAlRingMaterial,isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);
    innerCFRing.logical   = caloLogical(innerCFRing,  innerCFRingMaterial,isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);
    caseRing.logical      = caloLogical(caseRing,     shimMaterial,       isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);
    outerRailF.logical    = caloLogical(outerRailF,   outerRingMaterial,  isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);
    outerRailB.logical    = caloLogical(outerRailB,   outerRingMaterial,  isDiskVisible, G4Color::Yellow(), isDiskSolid, forceEdge);

    innerAlRingF.physical = caloPlacement(innerAlRingF, fullCase, 0, G4ThreeVector(0,0, diskCaseDZ-diskInAlRingDZ), false, 0, config, doSurfaceCheck, verbosity);
    innerAlRingM.physical = caloPlacement(innerAlRingM, fullCase, 0, G4ThreeVector(0,0,                         0), false, 0, config, doSurfaceCheck, verbosity);
    innerAlRingB.physical = caloPlacement(innerAlRingB, fullCase, 0, G4ThreeVector(0,0,-diskCaseDZ+diskInAlRingDZ), false, 0, config, doSurfaceCheck, verbosity);
    innerCFRing.physical  = caloPlacement(innerCFRing,  fullCase, 0, G4ThreeVector(0,0,0),                          false, 0, config, doSurfaceCheck, verbosity);
    caseRing.physical     = caloPlacement(caseRing,     fullCase, 0, G4ThreeVector(0,0,0),                          false, 0, config, doSurfaceCheck, verbosity);
    outerRailF.physical   = caloPlacement(outerRailF,   fullCase, 0, G4ThreeVector(0,0, diskCaseDZ-diskOutRailDZ),  false, 0, config, doSurfaceCheck, verbosity);
    outerRailB.physical   = caloPlacement(outerRailB,   fullCase, 0, G4ThreeVector(0,0,-diskCaseDZ+diskOutRailDZ),  false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(innerAlRingB);
    helper.addVolInfo(innerAlRingM);
    helper.addVolInfo(innerAlRingF);
    helper.addVolInfo(innerCFRing);
    helper.addVolInfo(caseRing);
    helper.addVolInfo(outerRailF);
    helper.addVolInfo(outerRailB);


    //------------------------------------------------------------
    // Two cooling pipes running between the rails. Make them similar to front plate cooling pipes.
    //
    G4RotationMatrix* rotPipe = reg.add(new G4RotationMatrix());
    rotPipe->rotateZ(0.1*CLHEP::pi);

    VolumeInfo coolPipeF("CaloCasePipeF_"+std::to_string(idisk));
    VolumeInfo coolPipeB("CaloCasePipeB_"+std::to_string(idisk));

    coolPipeF.solid    = new G4Torus(coolPipeF.name, FPCoolPipeRadius-FPCoolPipeThickness, FPCoolPipeRadius, diskCaseRingROut+FPCoolPipeRadius,0,1.2*CLHEP::pi);
    coolPipeB.solid    = new G4Torus(coolPipeB.name, FPCoolPipeRadius-FPCoolPipeThickness, FPCoolPipeRadius, diskCaseRingROut+FPCoolPipeRadius,0,1.2*CLHEP::pi);

    coolPipeF.logical  = caloLogical(coolPipeF, coolPipeMaterial, isDiskVisible, G4Color::Blue(), isDiskSolid, forceEdge);
    coolPipeB.logical  = caloLogical(coolPipeB, coolPipeMaterial, isDiskVisible, G4Color::Blue(), isDiskSolid, forceEdge);

    coolPipeF.physical = caloPlacement(coolPipeF, fullCase, rotPipe, G4ThreeVector(0,0, coolPipeZpos), false, 0, config, doSurfaceCheck, verbosity);;
    coolPipeB.physical = caloPlacement(coolPipeB, fullCase, rotPipe, G4ThreeVector(0,0,-coolPipeZpos), false, 0, config, doSurfaceCheck, verbosity);;

    helper.addVolInfo(coolPipeF);
    helper.addVolInfo(coolPipeB);


    //------------------------------------------------------------------
    // Inner and outer carbon fiber steps
    //
    std::vector<G4TwoVector> polygonOut = caloExtrudedVertices(stepsOutX,stepsOutY, 0.01);
    std::vector<G4TwoVector> polygonIn  = caloExtrudedVertices(stepsInX, stepsInY, -0.01);
    std::vector<G4TwoVector> polygonIn2 = caloExtrudedVertices(stepsInX, stepsInY, -diskStepThickness);

    VolumeInfo outerDiskSteps("CaloOuterDiskSteps_"+std::to_string(idisk));
    VolumeInfo innerDiskSteps("CaloInnerDiskSteps_"+std::to_string(idisk));
    VolumeInfo innerStepsHole("CaloInnerStepsHole_"+std::to_string(idisk));

    G4Tubs* innerStepHole        = new G4Tubs("caloInnerStepHole", 0, diskInCFRingROut, diskCaseDZ, 0,CLHEP::twopi);
    G4ExtrudedSolid* outerSteps  = new G4ExtrudedSolid("caloOuterSteps",  polygonOut, diskCaseDZ, 0,1,0,1);
    G4ExtrudedSolid* innerStepsO = new G4ExtrudedSolid("caloInnerStepsO", polygonIn,  diskCaseDZ, 0,1,0,1);
    G4ExtrudedSolid* innerStepsI = new G4ExtrudedSolid("caloInnerStepsI", polygonIn2, diskCaseDZ, 0,1,0,1);

    outerDiskSteps.solid         = new G4SubtractionSolid(outerDiskSteps.name, caseRing.solid, outerSteps,    0, G4ThreeVector(0,0,0));
    innerDiskSteps.solid         = new G4SubtractionSolid(innerDiskSteps.name, innerStepsO,    innerStepsI,   0, G4ThreeVector(0,0,0));
    innerStepsHole.solid         = new G4SubtractionSolid(innerStepsHole.name, innerStepsI,    innerStepHole, 0, G4ThreeVector(0,0,0));

    outerDiskSteps.logical       = caloLogical(outerDiskSteps, outerRingMaterial, isDiskVisible, G4Color::Red(), isDiskSolid, forceEdge);
    innerDiskSteps.logical       = caloLogical(innerDiskSteps, innerStepMaterial, isDiskVisible, G4Color::Red(), isDiskSolid, forceEdge);
    innerStepsHole.logical       = caloLogical(innerStepsHole, innerStepMaterial, isDiskVisible, G4Color::Red(), isCrystalSolid, forceEdge);

    outerDiskSteps.physical      = caloPlacement(outerDiskSteps, caseRing, 0, G4ThreeVector(0,0,0), false, 0, config, doSurfaceCheck, verbosity);
    innerDiskSteps.physical      = caloPlacement(innerDiskSteps, caseRing, 0, G4ThreeVector(0,0,0), false, 0, config, doSurfaceCheck, verbosity);
    innerStepsHole.physical      = caloPlacement(innerStepsHole, caseRing, 0, G4ThreeVector(0,0,0), false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(outerDiskSteps);
    helper.addVolInfo(innerDiskSteps);
    helper.addVolInfo(innerStepsHole);


    //-----------------------------------------
    // Wrapped crystal with reflector cap at the front and place them in the disk - for CsI and LYSO crystals
    // Create a single CsI / LYso crystal and place it multiple times to limit overhead -> wrapper.physical field and position are meaningless
    // Move the VolumeInfo creation inside the loop if there is a need to use individual crystal dimensions (not recommended)
    //
    VolumeInfo crystalCsI  ("CaloCrystalCsI_"+std::to_string(idisk));
    VolumeInfo crystalLYSO ("CaloCrystalLYSO_"+std::to_string(idisk));
    VolumeInfo wrapperCsI  ("CaloWrapperCsI"+std::to_string(idisk));
    VolumeInfo wrapperLYSO ("CaloWrapperLYSO"+std::to_string(idisk));
    VolumeInfo frontCapCsI ("CaloFrontCapCsI"+std::to_string(idisk));
    VolumeInfo frontCapLYSO("CaloFrontCapLYSO"+std::to_string(idisk));

    crystalCsI.solid      = new G4Box(crystalCsI.name,  crystalDXY, crystalDXY, crystalDZ);
    crystalLYSO.solid     = new G4Box(crystalLYSO.name, crystalDXY, crystalDXY, crystalDZ);
    frontCapCsI.solid     = new G4Box(frontCapCsI.name, crystalDXY, crystalDXY, crystalCapDZ);
    frontCapLYSO.solid    = new G4Box(frontCapLYSO.name,crystalDXY, crystalDXY, crystalCapDZ);
    wrapperCsI.solid      = new G4Box(wrapperCsI.name,  wrapperDXY, wrapperDXY, wrapperDZ);
    wrapperLYSO.solid     = new G4Box(wrapperLYSO.name, wrapperDXY, wrapperDXY, wrapperDZ);

    crystalCsI.logical    = caloLogical(crystalCsI,   crystalMaterial,       isCrystalVisible, G4Color::Cyan(),  isCrystalSolid, forceEdge);
    crystalLYSO.logical   = caloLogical(crystalLYSO,  caphriCrystalMaterial, isCrystalVisible, G4Color::Green(), isCrystalSolid, forceEdge);
    frontCapCsI.logical   = caloLogical(frontCapCsI,  crystalCapMaterial,    isCrystalVisible, G4Color::Cyan(),  isCrystalSolid, forceEdge);
    frontCapLYSO.logical  = caloLogical(frontCapLYSO, crystalCapMaterial,    isCrystalVisible, G4Color::Green(), isCrystalSolid, forceEdge);
    wrapperCsI.logical    = caloLogical(wrapperCsI,   wrapMaterial,          isCrystalVisible, G4Color::Cyan(),  isCrystalSolid, forceEdge);
    wrapperLYSO.logical   = caloLogical(wrapperLYSO,  wrapMaterial,          isCrystalVisible, G4Color::Green(), isCrystalSolid, forceEdge);

    auto crystalPos       = G4ThreeVector(0,0,crystalCapDZ);
    auto frontCapPos      = G4ThreeVector(0,0,crystalCapDZ-wrapperDZ);

    crystalCsI.physical   = caloPlacement(crystalCsI,   wrapperCsI,  0, crystalPos,  false, 0, config, doSurfaceCheck, verbosity);
    crystalLYSO.physical  = caloPlacement(crystalLYSO,  wrapperLYSO, 0, crystalPos,  false, 0, config, doSurfaceCheck, verbosity);
    frontCapCsI.physical  = caloPlacement(frontCapCsI,  wrapperCsI,  0, frontCapPos, false, 0, config, doSurfaceCheck, verbosity);
    frontCapLYSO.physical = caloPlacement(frontCapLYSO, wrapperLYSO, 0, frontCapPos, false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(crystalCsI);
    helper.addVolInfo(frontCapCsI);
    helper.addVolInfo(crystalLYSO);
    helper.addVolInfo(frontCapLYSO);
    helper.addVolInfo(wrapperCsI);
    helper.addVolInfo(wrapperLYSO);

    int nTotCrystal(0);
    for (unsigned i=0;i<idisk;++i) nTotCrystal+=cal.disk(idisk).nCrystals();
    for (size_t ic=0; ic<cal.disk(idisk).nCrystals(); ++ic) {

        G4int id = nTotCrystal+ic;
        CLHEP::Hep3Vector position = cal.disk(idisk).crystal(ic).localPosition();
        position.setZ(diskCaseDZ-wrapperDZ);

        bool isCaphri = std::find(caphriCystalId.begin(),caphriCystalId.end(),id) != caphriCystalId.end();
        if (isCaphri) caloPlacement(wrapperLYSO, caseRing, 0, position, true, id, config, doSurfaceCheck, verbosity);
        else          caloPlacement(wrapperCsI,  caseRing, 0, position, true, id, config, doSurfaceCheck, verbosity);
    }

    return fullCase;
  }


  //--------------------------------------------------------------------------------------------------------------------------------
  // build full backplate - yes this was annoying
  VolumeInfo caloBuildBackPlate(const SimpleConfig& config, unsigned idisk)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry( config, "calorimeterRO", "calorimeter.readout" );

    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper             = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry& reg            = helper.antiLeakRegistry();
    const DiskCalorimeter& cal       = *(GeomHandle<DiskCalorimeter>());

    const bool isROVisible           = geomOptions->isVisible          ("calorimeterRO");
    const bool isROSolid             = geomOptions->isSolid            ("calorimeterRO");
    const bool forceEdge             = geomOptions->forceAuxEdgeVisible("calorimeterRO");
    const bool doSurfaceCheck        = geomOptions->doSurfaceCheck     ("calorimeterRO");
    const int  verbosity             = config.getInt                   ("calorimeter.verbosityLevel",1);

    G4Material* vacuumMaterial       = materialFinder.get("calorimeter.vacuumMaterial");
    G4Material* ROMaterial           = materialFinder.get("calorimeter.readoutMaterial");
    G4Material* FEEMaterial          = materialFinder.get("calorimeter.FEEMaterial");
    G4Material* FEEBoxMaterial       = materialFinder.get("calorimeter.FEEBoxMaterial");
    G4Material* backPlateMaterial    = materialFinder.get("calorimeter.BackPlateMaterial");
    G4Material* pipeMaterial         = materialFinder.get("calorimeter.coolPipeMaterial");
    G4Material* stripMaterial        = materialFinder.get("calorimeter.BPStripMaterial");

    const double BPInnerRadius       = cal.caloInfo().getDouble("FPInnerRadius"); //same as front plate
    const double BPOuterRadius       = cal.caloInfo().getDouble("BPOuterRadius");
    const double diskCrystalRIn      = cal.caloInfo().getDouble("diskCrystalRIn");
    const double diskCrystalROut     = cal.caloInfo().getDouble("diskCrystalROut");

    const double crystalDXY          = cal.caloInfo().getDouble("crystalXYLength")/2.0;
    const double wrapperDXY          = crystalDXY + cal.caloInfo().getDouble("wrapperThickness");
    const double RODX                = cal.caloInfo().getDouble("readoutXLength")/2.0;
    const double RODY                = cal.caloInfo().getDouble("readoutYLength")/2.0;
    const double RODZ                = cal.caloInfo().getDouble("readoutZLength")/2.0;
    const double holeDX              = cal.caloInfo().getDouble("BPHoleXLength")/2.0;
    const double holeDY              = cal.caloInfo().getDouble("BPHoleYLength")/2.0;
    const double holeDZ              = cal.caloInfo().getDouble("BPHoleZLength")/2.0;
    const double stripDY             = wrapperDXY-holeDY-1.0;
    const double stripDZ             = cal.caloInfo().getDouble("BPStripThickness")/2.0;
    const double FEEDX               = cal.caloInfo().getDouble("FEEXLength")/2.0;
    const double FEEDY               = cal.caloInfo().getDouble("FEEYLength")/2.0;
    const double FEEDZ               = cal.caloInfo().getDouble("FEEZLength")/2.0;

    const double FEEBoxThickness     = cal.caloInfo().getDouble("FEEBoxThickness");
    const double FEEBoxDX            = holeDX + 2*FEEBoxThickness;
    const double FEEBoxDY            = FEEDY + 2*FEEBoxThickness;
    const double FEEBoxDZ            = FEEDZ + 2*FEEBoxThickness;

    const double BPPipeRadiusHigh    = cal.caloInfo().getDouble("BPPipeRadiusHigh");
    const double BPPipeRadiusLow     = cal.caloInfo().getDouble("BPPipeRadiusLow");
    const double BPPipeThickness     = cal.caloInfo().getDouble("BPPipeThickness");
    const double BPPipeTorRadiusHigh = BPOuterRadius - BPPipeRadiusHigh;
    const double BPPipeTorRadiusLow  = BPOuterRadius - 3.0*BPPipeRadiusHigh;
    const double BPPipeDZOffset      = cal.caloInfo().getDouble("BPPipeZOffset")/2.0;
    const double BPFEEDZ             = FEEBoxDZ + BPPipeDZOffset + BPPipeRadiusHigh;


    //----------------------
    // Full back plane, back plane and back plane FEE
    //----
    VolumeInfo fullBackPlate("CaloFullBackPlate_"+std::to_string(idisk));
    VolumeInfo backPlateFEE ("CaloBackPlateFEE_"+std::to_string(idisk));
    VolumeInfo backPlate    ("CaloBackPlate_"+std::to_string(idisk));

    fullBackPlate.solid   = new G4Tubs(fullBackPlate.name,BPInnerRadius,BPOuterRadius,BPFEEDZ+holeDZ,0,CLHEP::twopi);
    backPlateFEE.solid    = new G4Tubs(backPlateFEE.name, BPInnerRadius,BPOuterRadius,BPFEEDZ,       0,CLHEP::twopi);
    backPlate.solid       = new G4Tubs(backPlate.name,    BPInnerRadius,BPOuterRadius,holeDZ,        0,CLHEP::twopi);

    fullBackPlate.logical = caloLogical(fullBackPlate, vacuumMaterial,    0, G4Color::Black(),0,0);
    backPlateFEE.logical  = caloLogical(backPlateFEE,  vacuumMaterial,    0, G4Color::Black(),0,0);
    backPlate.logical     = caloLogical(backPlate,     backPlateMaterial, isROVisible,G4Color::Green(),isROSolid,forceEdge);

    backPlateFEE.physical = caloPlacement(backPlateFEE, fullBackPlate, 0, G4ThreeVector(0,0,holeDZ),   false, 0, config, doSurfaceCheck, verbosity);
    backPlate.physical    = caloPlacement(backPlate,    fullBackPlate, 0, G4ThreeVector(0,0,-BPFEEDZ), false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(backPlateFEE);
    helper.addVolInfo(backPlate);


    //-------------------------------------------------------------------------------------------
    // Build hole in back Plane, place SiPM at bottom of hole
    //----
    VolumeInfo holeBack ("CaloHoleBack_"+std::to_string(idisk));
    VolumeInfo crystalRO("CaloCrystalRO_"+std::to_string(idisk));

    holeBack.solid      = new G4Box(holeBack.name,  holeDX,holeDY,holeDZ);
    crystalRO.solid     = new G4Box(crystalRO.name, RODX,RODY,RODZ);

    holeBack.logical    = caloLogical(holeBack,  vacuumMaterial, 0,G4Color::Green(),0,0);
    crystalRO.logical   = caloLogical(crystalRO, ROMaterial,     isROVisible, G4Color::Yellow(),isROSolid, forceEdge);

    caloPlacement(crystalRO, holeBack, 0, G4ThreeVector( RODX, 0, -holeDZ+RODZ), true, 0, config, doSurfaceCheck, verbosity);
    caloPlacement(crystalRO, holeBack, 0, G4ThreeVector(-RODX, 0, -holeDZ+RODZ), true, 1, config, doSurfaceCheck, verbosity);

    int nTotCrystal(0);
    for (unsigned i=0;i<idisk;++i) nTotCrystal+=cal.disk(idisk).nCrystals();

    for(unsigned ic=0; ic <cal.disk(idisk).nCrystals(); ++ic)
    {
      G4int id = nTotCrystal+ic;
      CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).idealLocalPosition();
      unitPosition.setZ(0.0);

      caloPlacement(holeBack, backPlate, 0, unitPosition, true, id, config, doSurfaceCheck, verbosity);
    }
    helper.addVolInfo(holeBack);
    helper.addVolInfo(crystalRO);


    //----------------------
    // Build FEE in copper box
    //----
    VolumeInfo FEEBox  ("CaloFEEBox_"+std::to_string(idisk));
    VolumeInfo FEEBoxIn("CaloFEEBoxIn_"+std::to_string(idisk));
    VolumeInfo FEECard ("CaloFEECard_"+std::to_string(idisk));

    FEEBox.solid      = new G4Box(FEEBox.name,   FEEBoxDX,FEEBoxDY,FEEBoxDZ);
    FEEBoxIn.solid    = new G4Box(FEEBoxIn.name, FEEBoxDX-FEEBoxThickness,FEEBoxDY-FEEBoxThickness,FEEBoxDZ-0.5*FEEBoxThickness);
    FEECard.solid     = new G4Box(FEECard.name,  FEEDX,FEEDY,FEEDZ);

    FEEBox.logical    = caloLogical(FEEBox,   FEEBoxMaterial, isROVisible, G4Color::Yellow(),isROSolid, forceEdge);
    FEEBoxIn.logical  = caloLogical(FEEBoxIn, vacuumMaterial, 0, G4Color::Black(),0, 0);
    FEECard.logical   = caloLogical(FEECard,  FEEMaterial,    isROVisible, G4Color::Yellow(),isROSolid, forceEdge);


    // cards touch top of the box
    double dFEEsize   = FEEBoxDZ-0.5*FEEBoxThickness - FEEDZ;
    FEEBoxIn.physical = caloPlacement(FEEBoxIn, FEEBox, 0, G4ThreeVector(0, 0, -0.5*FEEBoxThickness), false, 0, config, doSurfaceCheck, verbosity);
    caloPlacement(FEECard, FEEBoxIn, 0, G4ThreeVector( RODX, 0, -dFEEsize), true, 0, config, doSurfaceCheck, verbosity);
    caloPlacement(FEECard, FEEBoxIn, 0, G4ThreeVector(-RODX, 0, -dFEEsize), true, 1, config, doSurfaceCheck, verbosity);

    for (unsigned ic=0; ic <cal.disk(idisk).nCrystals(); ++ic) {
      G4int id = nTotCrystal+ic;
      CLHEP::Hep3Vector unitPosition = cal.disk(idisk).crystal(ic).idealLocalPosition();
      unitPosition.setZ(-BPFEEDZ + FEEBoxDZ);

      caloPlacement(FEEBox, backPlateFEE, 0, unitPosition, true, id, config, doSurfaceCheck, verbosity);
    }

    helper.addVolInfo(FEEBox);
    helper.addVolInfo(FEEBoxIn);
    helper.addVolInfo(FEECard);


    //----------------------
    // Add cooling strips at the back of the backPlate.
    // Approximate model used to avoid recording the length of each strip - good enough in our case
    //----
    double yminStrip(diskCrystalROut);
    for(unsigned ic=0; ic <cal.disk(idisk).nCrystals(); ++ic) {
      CLHEP::Hep3Vector position = cal.disk(idisk).crystal(ic).idealLocalPosition();
      yminStrip = std::min(yminStrip,position.y()-wrapperDXY);
    }

    int istrip(0);
    double yStrip(yminStrip);
    while (yStrip < 0){
      double zpos = holeDZ-stripDZ;
      double ymin = (yStrip<0) ? yStrip+stripDY : yStrip-stripDY;
      double xmin = (abs(yStrip) < diskCrystalRIn) ? sqrt(diskCrystalRIn*diskCrystalRIn-ymin*ymin) : 0;
      double xmax = sqrt(diskCrystalROut*diskCrystalROut-ymin*ymin);

      if (xmin > 0)
      {
        double stripDX = (xmax-xmin)/2.0;
        double xpos    = xmin+stripDX;

        VolumeInfo strip("Calostrip_"+std::to_string(idisk*100+istrip));
        strip.solid      = new G4Box(strip.name,stripDX,stripDY,stripDZ);
        strip.logical    = caloLogical(strip, stripMaterial, isROVisible, G4Color::Red(),isROSolid, forceEdge);

        caloPlacement(strip, backPlate, 0, G4ThreeVector( xpos,  yStrip, zpos), true, 0, config, doSurfaceCheck, verbosity);
        caloPlacement(strip, backPlate, 0, G4ThreeVector(-xpos,  yStrip, zpos), true, 1, config, doSurfaceCheck, verbosity);
        caloPlacement(strip, backPlate, 0, G4ThreeVector( xpos, -yStrip, zpos), true, 2, config, doSurfaceCheck, verbosity);
        caloPlacement(strip, backPlate, 0, G4ThreeVector(-xpos, -yStrip, zpos), true, 3, config, doSurfaceCheck, verbosity);

        helper.addVolInfo(strip);
      }
      else
      {
        double stripDX = xmax-xmin;
        VolumeInfo strip("Calostrip_"+std::to_string(idisk*100+istrip));
        strip.solid      = new G4Box(strip.name,stripDX,stripDY,stripDZ);
        strip.logical    = caloLogical(strip,stripMaterial, isROVisible, G4Color::Red(),isROSolid, forceEdge);

        caloPlacement(strip, backPlate, 0, G4ThreeVector(0, yStrip, zpos), true, 0, config, doSurfaceCheck, verbosity);
        caloPlacement(strip, backPlate, 0, G4ThreeVector(0, -yStrip, zpos), true, 1, config, doSurfaceCheck, verbosity);

        helper.addVolInfo(strip);
      }
      yStrip += 2*wrapperDXY;
      ++istrip;
    }

    //----------------------
    // Build last cooling pipes
    // Trick: pipe in xy plane going from -alpha to beta = pipe going from pi-beta to pi+alpha, then rotated by pi around y axis
    //----
    VolumeInfo BPPipe1("CaloBPPipe1_"+std::to_string(idisk));
    VolumeInfo BPPipe2("CaloBPPipe2_"+std::to_string(idisk));
    VolumeInfo BPPipe3("CaloBPPipe3_"+std::to_string(idisk));
    VolumeInfo BPPipe4("CaloBPPipe4_"+std::to_string(idisk));

    BPPipe1.solid    = new G4Torus(BPPipe1.name,BPPipeRadiusHigh-BPPipeThickness, BPPipeRadiusHigh, BPPipeTorRadiusHigh, 0.5*CLHEP::pi, 0.9*CLHEP::pi);
    BPPipe2.solid    = new G4Torus(BPPipe2.name,BPPipeRadiusLow-BPPipeThickness,  BPPipeRadiusLow,  BPPipeTorRadiusHigh, 0.5*CLHEP::pi, 0.7*CLHEP::pi);
    BPPipe3.solid    = new G4Torus(BPPipe3.name,BPPipeRadiusHigh-BPPipeThickness, BPPipeRadiusHigh, BPPipeTorRadiusLow,  0.5*CLHEP::pi, 0.9*CLHEP::pi);
    BPPipe4.solid    = new G4Torus(BPPipe4.name,BPPipeRadiusLow-BPPipeThickness,  BPPipeRadiusLow,  BPPipeTorRadiusLow,  0.5*CLHEP::pi, 0.7*CLHEP::pi);

    BPPipe1.logical  = caloLogical(BPPipe1,pipeMaterial, isROVisible, G4Color::Blue(),isROSolid, forceEdge);
    BPPipe2.logical  = caloLogical(BPPipe2,pipeMaterial, isROVisible, G4Color::Blue(),isROSolid, forceEdge);
    BPPipe3.logical  = caloLogical(BPPipe3,pipeMaterial, isROVisible, G4Color::Red(),isROSolid, forceEdge);
    BPPipe4.logical  = caloLogical(BPPipe4,pipeMaterial, isROVisible, G4Color::Red(),isROSolid, forceEdge);

    G4RotationMatrix* rotY = reg.add(new G4RotationMatrix(CLHEP::HepRotation::IDENTITY));
    rotY->rotateY(CLHEP::pi);
    BPPipe1.physical = caloPlacement(BPPipe1, backPlateFEE, rotY, G4ThreeVector(0,0,BPFEEDZ-BPPipeRadiusHigh), false, 0, config, doSurfaceCheck, verbosity);
    BPPipe2.physical = caloPlacement(BPPipe2, backPlateFEE,    0, G4ThreeVector(0,0,BPFEEDZ-BPPipeRadiusHigh), false, 0, config, doSurfaceCheck, verbosity);
    BPPipe3.physical = caloPlacement(BPPipe3, backPlateFEE,    0, G4ThreeVector(0,0,BPFEEDZ-BPPipeRadiusHigh), false, 0, config, doSurfaceCheck, verbosity);
    BPPipe4.physical = caloPlacement(BPPipe4, backPlateFEE, rotY, G4ThreeVector(0,0,BPFEEDZ-BPPipeRadiusHigh), false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(BPPipe1);
    helper.addVolInfo(BPPipe2);
    helper.addVolInfo(BPPipe3);
    helper.addVolInfo(BPPipe4);

    return fullBackPlate;
  }




  //--------------------------------------------------------------------------------------------------------------------------------
  // build full FEB - some parameters are common with the detector solenoid and taken from the corresponding class
  VolumeInfo caloBuildFEB(const SimpleConfig& config, unsigned idisk)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry(config, "calorimeterCrate", "calorimeter.crate" );

    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper           = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry& reg          = helper.antiLeakRegistry();
    const DiskCalorimeter& cal     = *(GeomHandle<DiskCalorimeter>());
    const DetectorSolenoid& ds     = *(GeomHandle<DetectorSolenoid>());

    G4Material* vacuumMaterial     = materialFinder.get("calorimeter.vacuumMaterial");
    G4Material* cableMaterial      = findMaterialOrThrow(ds.calCableRunMaterial());
    G4Material* cableCoreMaterial  = findMaterialOrThrow(ds.materialCableRunCalCore());

    const bool isCrateVisible      = geomOptions->isVisible("calorimeterCrate");
    const bool isCrateSolid        = geomOptions->isSolid("calorimeterCrate");
    const bool forceEdge           = geomOptions->forceAuxEdgeVisible("calorimeterCrate");
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("calorimeterCrate");
    const int  verbosity           = config.getInt("calorimeter.verbosityLevel",1);

    const int    nCrates           = cal.caloInfo().getInt("nCrates");
    const double vdThickness       = cal.caloInfo().getDouble("vdThickness");
    const double caloDiskRadiusOut = cal.caloInfo().getDouble("caloDiskRadiusOut");
    const auto   cratePhiAngle     = cal.caloInfo().getVDouble("cratePhiAngles");

    // get me a lil' crate
    VolumeInfo crateBox = caloBuildCrate(config,idisk);
    G4Box* crate = dynamic_cast<G4Box*>(crateBox.solid);

    auto  FEBPhiMinMax             = calcFEBPhiRange(cal);
    G4double crateXHalfLength      = crate->GetXHalfLength();
    G4double crateYHalfLength      = crate->GetYHalfLength();
    G4double crateZHalfLength      = crate->GetZHalfLength();
    G4double crateMaxRadius        = caloDiskRadiusOut+ vdThickness+2*crateYHalfLength;
    G4double crateRadius           = sqrt(crateMaxRadius*crateMaxRadius+crateXHalfLength*crateXHalfLength)-caloDiskRadiusOut;
    G4double FEBRadIn              = caloDiskRadiusOut;
    G4double cratePosY             = FEBRadIn + crateYHalfLength + vdThickness;
    G4double FEBRadOut             = FEBRadIn + crateRadius + 2.0*vdThickness;
    G4double FEBZHalfLength        = crateZHalfLength + vdThickness;

    VolumeInfo fullFEB("CaloFEB_"+std::to_string(idisk));
    fullFEB.solid   = new G4Tubs(fullFEB.name,FEBRadIn, FEBRadOut,FEBZHalfLength,FEBPhiMinMax[0],FEBPhiMinMax[1]-FEBPhiMinMax[0]);
    fullFEB.logical = caloLogical(fullFEB, vacuumMaterial, 0, G4Color::Black(), 0, 0);


    for (G4int icrt=0;icrt < nCrates; ++icrt) {
       double phiCrate = cratePhiAngle[icrt]*CLHEP::degree;

       G4RotationMatrix* rotCrate = reg.add(new G4RotationMatrix());
       rotCrate->rotateZ(CLHEP::pi/2-phiCrate);

       G4ThreeVector posCrate = G4ThreeVector(cratePosY*std::cos(phiCrate),cratePosY*std::sin(phiCrate),0.0);
       caloPlacement(crateBox, fullFEB, rotCrate,posCrate,true,icrt,config,doSurfaceCheck,verbosity);
    }
    helper.addVolInfo(crateBox);

    if (ds.hasCableRunCal()) {
       G4double crRin  = ds.upRInCableRunCal();
       G4double crRout = ds.upROutCableRunCal();
       G4double dPhi   = ds.dPhiCableRunCal()*CLHEP::degree;
       G4double phi0   = CLHEP::pi/2-dPhi/2.0;

       VolumeInfo cableRunCal("CaloCableRunCal_"+std::to_string(idisk));
       cableRunCal.solid    = new G4Tubs(cableRunCal.name, crRin, crRout, crateZHalfLength - vdThickness, phi0, dPhi);
       cableRunCal.logical  = caloLogical(cableRunCal,   cableMaterial, isCrateVisible,G4Color::Magenta(),isCrateSolid,forceEdge);
       cableRunCal.physical = caloPlacement(cableRunCal, fullFEB, 0,G4ThreeVector(0,0,0),false,0,config,doSurfaceCheck,verbosity);
       helper.addVolInfo(cableRunCal);

       if (ds.cableRunVersion()>2) {
         TubsParams caloCableRunCalParams(crRin, crRout, crateZHalfLength, phi0, dPhi);
         TubsParams cableRunCorePars    = calculateTubeCoreParams(caloCableRunCalParams,ds.rCableRunCalCoreFract(),
                                                                  ds.rdCableRunCalCoreFract(),ds.dPhiCableRunCalCoreFract());
         const auto& pars               = cableRunCorePars.data();
         VolumeInfo cableRunCalCore("CaloCableRunCalCore_"+std::to_string(idisk));
         cableRunCalCore.solid    = new G4Tubs(cableRunCalCore.name,pars[0],pars[1],0.999*pars[2],pars[3],pars[4]);
         cableRunCalCore.logical  = caloLogical(cableRunCalCore, cableCoreMaterial, isCrateVisible, G4Color::Grey(),isCrateSolid,forceEdge);
         cableRunCalCore.physical = caloPlacement(cableRunCalCore, cableRunCal, 0,G4ThreeVector(0,0,0),false,0,config,doSurfaceCheck,verbosity);
         helper.addVolInfo(cableRunCalCore);
       }
     }

     double delta = dynamic_cast<G4Tubs*>(fullFEB.solid)->GetZHalfLength()*2 - cal.disk(idisk).geomInfo().FEBZLength();

     if (std::abs(delta) > 1e-3)  G4cout << __func__ <<"PANIC..... geometry description in Geant4 and DiskMaker do NOT match  - FEB size: "
                                         <<dynamic_cast<G4Tubs*>(fullFEB.solid)->GetZHalfLength()*2<<" vs "
                                         <<cal.disk(idisk).geomInfo().FEBZLength()<<G4endl;

     if (verbosity) G4cout << __func__ <<" Compare FEB size  Geant4 / CaloInfo "<<dynamic_cast<G4Tubs*>(fullFEB.solid)->GetZHalfLength()*2
                                <<" / "<<cal.disk(idisk).geomInfo().FEBZLength()<<G4endl;

     return fullFEB;
  }


  //--------------------------------------------------------------------------------------------------------------------------------
  // build crate
  VolumeInfo caloBuildCrate(const SimpleConfig& config, unsigned idisk)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry(config, "calorimeterCrate",      "calorimeter.crate" );
    geomOptions->loadEntry(config, "calorimeterCrateBoard", "calorimeter.crateBoard" );

    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper             = *(art::ServiceHandle<Mu2eG4Helper>());
    const DiskCalorimeter& cal       = *(GeomHandle<DiskCalorimeter>());

    const bool isCrateVisible        = geomOptions->isVisible          ("calorimeterCrate");
    const bool isBoardVisible        = geomOptions->isVisible          ("calorimeterCrateBoard");
    const bool isCrateSolid          = geomOptions->isSolid            ("calorimeterCrate");
    const bool isBoardSolid          = geomOptions->isSolid            ("calorimeterCrateBoard");
    const bool forceEdge             = geomOptions->forceAuxEdgeVisible("calorimeterCrate");
    const bool doSurfaceCheck        = geomOptions->doSurfaceCheck     ("calorimeterCrate");
    const int  verbosity             = config.getInt                   ("calorimeter.verbosityLevel",1);

    G4Material* vacuumMaterial       = materialFinder.get("calorimeter.vacuumMaterial");
    G4Material* crateMaterial        = materialFinder.get("calorimeter.crateMaterial");
    G4Material* crateShieldMaterial  = materialFinder.get("calorimeter.shieldMaterial");
    G4Material* radiatorMaterial     = materialFinder.get("calorimeter.radiatorMaterial");
    G4Material* activeStripMaterial  = materialFinder.get("calorimeter.activeStripMaterial");
    G4Material* passiveStripMaterial = materialFinder.get("calorimeter.passiveStripMaterial");

    const int nBoards                = cal.caloInfo().getInt   ("nBoards");
    const double crateDX             = cal.caloInfo().getDouble("crateXLength")/2.0;
    const double crateDY             = cal.caloInfo().getDouble("crateYLength")/2.0;
    const double crateDZ             = cal.caloInfo().getDouble("crateZLength")/2.0;
    const double crateFShieldDisp    = cal.caloInfo().getDouble("crateFShieldDeltaZ");
    const double crateFShieldThick   = cal.caloInfo().getDouble("crateFShieldThickness");
    const double crateBottomThick    = cal.caloInfo().getDouble("crateBShieldThickness");
    const double crateBottomLength   = cal.caloInfo().getDouble("crateBShieldLength");
    const double crateTopThick       = cal.caloInfo().getDouble("crateTThickness");
    const double crateSideThick      = cal.caloInfo().getDouble("crateSThickness");
    const double crateFShieldDY      = cal.caloInfo().getDouble("crateFShieldYLength")/2.0;

    const double radiatorDY          = cal.caloInfo().getDouble("radiatorThickness")/2.0;
    const double radiatorDZ          = cal.caloInfo().getDouble("radiatorZLength")/2.0;
    const double activeStripDY       = cal.caloInfo().getDouble("activeStripThickness")/2.0;
    const double passiveStripDY      = cal.caloInfo().getDouble("passiveStripThickness")/2.0;
    const double crateFullDZ         = crateDZ+crateFShieldDisp/2.0+crateFShieldThick/2.0;
    const double crateFullDY         = crateDY+crateBottomThick/2.0;
    const double crateBoxInDY        = crateDY-crateTopThick/2.0;
    const double boardDX             = crateDX-crateSideThick;
    const double boardDY             = radiatorDY+activeStripDY+passiveStripDY;
    const double boardDZ             = crateDZ;
    const double radiatorPosY        = -boardDY + radiatorDY;
    const double activeStripPosY     = radiatorPosY+radiatorDY+activeStripDY;
    const double passiveStripPosY    = activeStripPosY + activeStripDY + passiveStripDY;
    const double deltaY              = 2.0*crateBoxInDY/float(nBoards+1);

    // ---------------------------------------
    // define the crate box and shield pieces
    //
    VolumeInfo ccrateFullBox   ("ccrateFullBox_"+std::to_string(idisk));
    VolumeInfo ccrateBoxTop    ("ccrateBoxTop_"+std::to_string(idisk));
    VolumeInfo ccrateBoxIn     ("ccrateBoxIn_"+std::to_string(idisk));
    VolumeInfo ccrateBoxBottom ("ccrateBoxBottom_"+std::to_string(idisk));
    VolumeInfo ccrateBoxShieldF("ccrateBoxShieldF_"+std::to_string(idisk));

    ccrateFullBox.solid       = new G4Box(ccrateFullBox.name,    crateDX,crateFullDY,crateFullDZ);
    ccrateBoxTop.solid        = new G4Box(ccrateBoxTop.name,     crateDX,crateDY,crateDZ);
    ccrateBoxIn.solid         = new G4Box(ccrateBoxIn.name,      crateDX-crateSideThick,crateBoxInDY,crateDZ);
    ccrateBoxBottom.solid     = new G4Box(ccrateBoxBottom.name,  crateDX,crateBottomThick/2.0,crateBottomLength/2.);
    ccrateBoxShieldF.solid    = new G4Box(ccrateBoxShieldF.name, crateDX,crateFShieldDY,crateFShieldThick/2.0);

    ccrateFullBox.logical     = caloLogical(ccrateFullBox,    vacuumMaterial,      0,              G4Color::Black(),  0, 0);
    ccrateBoxTop.logical      = caloLogical(ccrateBoxTop,     crateMaterial,       isCrateVisible, G4Color::Blue(),   isCrateSolid, forceEdge);
    ccrateBoxIn.logical       = caloLogical(ccrateBoxIn,      vacuumMaterial,      0,              G4Color::Black(),  0, 0);
    ccrateBoxBottom.logical   = caloLogical(ccrateBoxBottom,  crateShieldMaterial, isCrateVisible, G4Color::Yellow(), isCrateSolid, forceEdge);
    ccrateBoxShieldF.logical  = caloLogical(ccrateBoxShieldF, crateShieldMaterial, isCrateVisible, G4Color::Yellow(), isCrateSolid, forceEdge);

    auto ccBoxTopPos          = G4ThreeVector(0.0,crateBottomThick/2.0,crateFullDZ-crateDZ);
    auto ccBoxInPos           = G4ThreeVector(0,-crateTopThick/2.0,0);
    auto ccBoxBottomPos       = G4ThreeVector(0,-crateFullDY + crateBottomThick/2.0,-crateFullDZ + crateFShieldThick + crateBottomLength/2.0);
    auto ccBoxShieldFPos      = G4ThreeVector(0,-crateFullDY + crateFShieldDY,-crateFullDZ + crateFShieldThick/2.0);

    ccrateBoxTop.physical     = caloPlacement(ccrateBoxTop,     ccrateFullBox, 0, ccBoxTopPos,     false, 0, config, doSurfaceCheck, verbosity);
    ccrateBoxIn.physical      = caloPlacement(ccrateBoxIn,      ccrateBoxTop,  0, ccBoxInPos,      false, 0, config, doSurfaceCheck, verbosity);
    ccrateBoxBottom.physical  = caloPlacement(ccrateBoxBottom,  ccrateFullBox, 0, ccBoxBottomPos,  false, 0, config, doSurfaceCheck, verbosity);
    ccrateBoxShieldF.physical = caloPlacement(ccrateBoxShieldF, ccrateFullBox, 0, ccBoxShieldFPos, false, 0, config, doSurfaceCheck, verbosity);

    helper.addVolInfo(ccrateBoxTop);
    helper.addVolInfo(ccrateBoxIn);
    helper.addVolInfo(ccrateBoxBottom);
    helper.addVolInfo(ccrateBoxShieldF);

    // ---------------
    // add the boards
    //
    VolumeInfo ccrateBoard     ("ccrateBoard_"+std::to_string(idisk));
    VolumeInfo ccrateRadiator  ("ccrateRadiator_"+std::to_string(idisk));
    VolumeInfo ccrateActiveStr ("ccrateActiveStrip_"+std::to_string(idisk));
    VolumeInfo ccratePassiveStr("ccratePassiveStrip_"+std::to_string(idisk));

    ccrateBoard.solid         = new G4Box(ccrateBoard.name,      boardDX, boardDY,        boardDZ);
    ccrateRadiator.solid      = new G4Box(ccrateRadiator.name,   boardDX, radiatorDY,     radiatorDZ);
    ccrateActiveStr.solid     = new G4Box(ccrateActiveStr.name,  boardDX, activeStripDY,  boardDZ);
    ccratePassiveStr.solid    = new G4Box(ccratePassiveStr.name, boardDX, passiveStripDY, boardDZ);

    ccrateBoard.logical       = caloLogical(ccrateBoard,      vacuumMaterial,       0,              G4Color::Black(),0, 0);
    ccrateRadiator.logical    = caloLogical(ccrateRadiator,   radiatorMaterial,     isBoardVisible, G4Color::Green(), isBoardSolid, forceEdge);
    ccrateActiveStr.logical   = caloLogical(ccrateActiveStr,  activeStripMaterial,  isBoardVisible, G4Color::Red(),   isBoardSolid, forceEdge);
    ccratePassiveStr.logical  = caloLogical(ccratePassiveStr, passiveStripMaterial, isBoardVisible, G4Color::Green(), isBoardSolid, forceEdge);

    auto ccRadiatorPos        = G4ThreeVector(0.0,radiatorPosY,-boardDZ+radiatorDZ);
    auto ccActivePos          = G4ThreeVector(0.0,activeStripPosY,0.0);
    auto ccPassivePos         = G4ThreeVector(0.0,passiveStripPosY,0.0);
    ccrateRadiator.physical   = caloPlacement(ccrateRadiator,   ccrateBoard, 0, ccRadiatorPos, false, 0, config, doSurfaceCheck, verbosity);
    ccrateActiveStr.physical  = caloPlacement(ccrateActiveStr,  ccrateBoard, 0, ccActivePos,   false, 0, config, doSurfaceCheck, verbosity);
    ccratePassiveStr.physical = caloPlacement(ccratePassiveStr, ccrateBoard, 0, ccPassivePos,  false, 0, config, doSurfaceCheck, verbosity);

    for (G4int ibrd=0; ibrd < nBoards; ++ibrd) {
      caloPlacement(ccrateBoard, ccrateBoxIn, 0, G4ThreeVector(0.0,(ibrd+1)*deltaY-crateBoxInDY,0), true, ibrd, config, doSurfaceCheck, verbosity);
    }

    helper.addVolInfo(ccrateRadiator);
    helper.addVolInfo(ccrateActiveStr);
    helper.addVolInfo(ccratePassiveStr);
    helper.addVolInfo(ccrateBoard);

    return ccrateFullBox;
  }



  //--------------------------------------------------------------------------------------------------------------------------------
  // build calo cable run - some parameters are common with the detector solenoid and taken from the corresponding class
  VolumeInfo caloBuildCable(const SimpleConfig& config, unsigned idisk, const VolumeInfo& FEBvol)
  {
    const auto geomOptions(art::ServiceHandle<GeometryService>()->geomOptions());
    geomOptions->loadEntry(config, "calorimeterCrate", "calorimeter.crate" );

    MaterialFinder materialFinder(config);

    Mu2eG4Helper& helper          = *(art::ServiceHandle<Mu2eG4Helper>());
    const DiskCalorimeter& cal    = *(GeomHandle<DiskCalorimeter>());
    const DetectorSolenoid& ds    = *(GeomHandle<DetectorSolenoid>());

    const bool isCrateVisible     = geomOptions->isVisible("calorimeterCrate");
    const bool isCrateSolid       = geomOptions->isSolid("calorimeterCrate");
    const bool forceEdge          = geomOptions->forceAuxEdgeVisible("calorimeterCrate");
    const bool doSurfaceCheck     = geomOptions->doSurfaceCheck("calorimeterCrate");
    const int  verbosity          = config.getInt("calorimeter.verbosityLevel",1);

    G4Material* cableMaterial     = findMaterialOrThrow(ds.calCableRunMaterial());
    G4Material* cableCoreMaterial = findMaterialOrThrow(ds.materialCableRunCalCore());

    const double mother_z0        = cal.caloInfo().getDouble("caloMotherZ0");
    const double mother_z1        = cal.caloInfo().getDouble("caloMotherZ1");
    const double mother_zlength   = mother_z1-mother_z0;
    const double crRin            = ds.upRInCableRunCal();
    const double crRout           = ds.upROutCableRunCal();
    const double dPhi             = ds.dPhiCableRunCal()*CLHEP::degree;
    const double phi0             = CLHEP::pi/2-dPhi/2.0;


    //length of the cable run = distance between disk for first gaps, ditance betweenm end of feb and mother volume for last gap
    G4Tubs*  FEB             = dynamic_cast<G4Tubs*>(FEBvol.solid);
    G4double distFEBtoNext   = (idisk+1 < cal.nDisks()) ? cal.disk(idisk+1).geomInfo().origin().z() - cal.disk(idisk).geomInfo().origin().z()
                                                       : mother_zlength/2.0 - FEBvol.physical->GetTranslation().z()-FEB->GetZHalfLength();
    G4double cableHalfLength = (idisk+1 < cal.nDisks()) ? 0.5*distFEBtoNext - FEB->GetZHalfLength() : 0.5*distFEBtoNext;

    VolumeInfo cableRunGap("CaloCableRunGap_"+std::to_string(idisk));
    cableRunGap.solid   = new G4Tubs(cableRunGap.name,crRin, crRout, cableHalfLength, phi0, dPhi);
    cableRunGap.logical = caloLogical(cableRunGap, cableMaterial, isCrateVisible,G4Color::Magenta(),isCrateSolid,forceEdge);

    if (ds.cableRunVersion()>2) {
      TubsParams caloCableRunCalParams(crRin, crRout, cableHalfLength, phi0, dPhi);
      TubsParams cableRunCorePars = calculateTubeCoreParams(caloCableRunCalParams,ds.rCableRunCalCoreFract(),
                                                            ds.rdCableRunCalCoreFract(),ds.dPhiCableRunCalCoreFract());
      const auto pars = cableRunCorePars.data();
      VolumeInfo cableRunGapCore("CaloCableRunGapCore_"+std::to_string(idisk));
      cableRunGapCore.solid    = new G4Tubs(cableRunGapCore.name,pars[0],pars[1],0.999*pars[2],pars[3],pars[4]);
      cableRunGapCore.logical  = caloLogical(cableRunGapCore, cableCoreMaterial, isCrateVisible, G4Color::Grey(),isCrateSolid,forceEdge);
      cableRunGapCore.physical = caloPlacement(cableRunGapCore, cableRunGap, 0,G4ThreeVector(0,0,0),false,0,config,doSurfaceCheck,verbosity);
      helper.addVolInfo(cableRunGapCore);
    }

    return cableRunGap;
  }




  //--------------------------------------------------------------------------------------------------------------------------------
  // make logical volume and attach visible attribute to it
  G4LogicalVolume* caloLogical(const VolumeInfo& volume, G4Material* mat, bool isVisible, const G4Color& color, bool isSolid, bool forceEdge)
  {
     G4LogicalVolume* logical = new G4LogicalVolume(volume.solid, mat, volume.name);

     AntiLeakRegistry& reg    = (art::ServiceHandle<Mu2eG4Helper>())->antiLeakRegistry();
     G4VisAttributes* visAtt  = reg.add(G4VisAttributes(true, color));

     if (!isVisible) {
       logical->SetVisAttributes(G4VisAttributes::GetInvisible());
     } else {
      visAtt->SetForceSolid(isSolid);
      visAtt->SetForceAuxEdgeVisible(forceEdge);
      logical->SetVisAttributes(visAtt);
     }
     return logical;
  }

  //--------------------------------------------------------------------------------------------------------------------------------
  // place volumeInfo and calculate their position
  G4PVPlacement* caloPlacement(VolumeInfo& volume, const VolumeInfo& parent, G4RotationMatrix* rot, const G4ThreeVector& position,
                               bool pMany, int copyNo, const SimpleConfig& config, bool doSurfaceCheck, int verbosity)
  {
     if (volume.logical == nullptr) return nullptr;
     G4PVPlacement* physical = new G4PVPlacement(rot, position, volume.logical, volume.name, parent.logical, pMany, copyNo, false);
     doSurfaceCheck && checkForOverlaps(physical, config, verbosity>0);

     if (!pMany){
       volume.centerInParent = position;
       caloVolInfG4[volume.name] = parent.name;
       //volume.centerInWorld  = parent.centerInWorld + position;
     }
     return physical;
  }

  //--------------------------------------------------------------------------------------------------------------------------------
  // utility to make list of vertexes for extruded polygon solid
  std::vector<G4TwoVector> caloExtrudedVertices(const std::vector<double>& stepsX, const std::vector<double>& stepsY, double delta)
  {
    std::vector<double> stepsNewX,stepsNewY;
    for (const auto& val : stepsX) stepsNewX.push_back(val>0 ? val+delta : val-delta);
    for (const auto& val : stepsY) stepsNewY.push_back(val>0 ? val+delta : val-delta);

    std::vector<G4TwoVector> polygon,polygon2;
    for (size_t i=0;i<stepsNewX.size();i+=2)
    {
       if (i==0 || abs(stepsNewX[i]-polygon.back().x()) > 1e-3 ){
         polygon.push_back(G4TwoVector(stepsNewX[i],stepsNewY[i]));
         polygon.push_back(G4TwoVector(stepsNewX[i],stepsNewY[i+1]));
       } else {
         polygon.back() = G4TwoVector(stepsNewX[i],stepsNewY[i+1]);
       }
       if (i==0 || abs(stepsNewX[i+1]-polygon2.back().x()) > 1e-3 ){
         polygon2.push_back(G4TwoVector(stepsNewX[i+1],stepsNewY[i]));
         polygon2.push_back(G4TwoVector(stepsNewX[i+1],stepsNewY[i+1]));
       } else {
         polygon2.back() = G4TwoVector(stepsNewX[i+1],stepsNewY[i+1]);
       }
    }
    std::reverse_copy(polygon2.begin(), polygon2.end(), std::back_inserter(polygon));
    return polygon;
  }

  //--------------------------------------------------------------------------------------------------------------------------------
  // utility to make list of vertexes for extruded polygon solid
  std::vector<G4double> calcFEBPhiRange(const DiskCalorimeter& cal)
  {
    G4double crateXLength  = cal.caloInfo().getDouble("crateXLength");
    G4double crateRin      = cal.caloInfo().getDouble("caloDiskRadiusOut");
    auto cratePhis         = cal.caloInfo().getVDouble("cratePhiAngles");

    for (auto& val : cratePhis) val *= CLHEP::degree;
    G4double cratePhiSpan  = atan(crateXLength/2.0/crateRin);
    G4double phiMin        = (*std::min_element(cratePhis.begin(),cratePhis.end()))-cratePhiSpan-0.02;
    G4double phiMax        = (*std::max_element(cratePhis.begin(),cratePhis.end()))+cratePhiSpan+0.02;

    std::vector<G4double> range{phiMin,phiMax};
    return range;
  }

  //--------------------------------------------------------------------------------------------------------------------------------
  // utility to traverse G4LogicalVolumes using recursion
  const G4LogicalVolume* findCaloSolid(const G4LogicalVolume* volume, const G4String& objectName, std::vector<const G4LogicalVolume*>& nodes)
  {
     if (volume->GetSolid()->GetName() == objectName) return volume;

     std::queue<G4LogicalVolume*> toProcess;
     for (size_t i=0;i<volume->GetNoDaughters();++i) toProcess.push(volume->GetDaughter(i)->GetLogicalVolume());

     while (!toProcess.empty()){
       const G4LogicalVolume* result = findCaloSolid(toProcess.front(),objectName,nodes);
       if (result != nullptr){
          nodes.push_back(toProcess.front());
          return result;
       }
       toProcess.pop();
     }
     return nullptr;
  }


  //--------------------------------------------------------------------------------------------------------------------------------
  // utility to traverse G4LogicalVolumes using recursion
  void browseCaloSolids(const G4LogicalVolume* volume)
  {

     if (G4Box* obj = dynamic_cast<G4Box*>(volume->GetSolid()))
        std::cout<<obj->GetName()<<",XXXX,G4Box,dX="<<2*obj->GetXHalfLength()<<",dY="<<2*obj->GetYHalfLength()
                 <<",dZ="<<2*obj->GetZHalfLength()<<","<<volume->GetMaterial()->GetName()<<std::endl;

     else if (G4Tubs* obj = dynamic_cast<G4Tubs*>(volume->GetSolid()))
        std::cout<<obj->GetName()<<",XXXX,G4Tubs,Rin="<<obj->GetInnerRadius()<<",Rout="<<obj->GetOuterRadius()
                 <<",dZ="<<2*obj->GetZHalfLength()<<",dPhi="<<obj->GetDeltaPhiAngle()<<","<<volume->GetMaterial()->GetName()<<std::endl;

     else if (G4Torus* obj = dynamic_cast<G4Torus*>(volume->GetSolid()))
        std::cout<<obj->GetName()<<",XXXX,G4Torus,Rmin="<<obj->GetRmin()<<",Rmax="<<obj->GetRmax()<<",Rtot="
                 <<obj->GetRtor()<<",dPhi="<<obj->GetDPhi()<<","<<volume->GetMaterial()->GetName()<<std::endl;

     else if (volume != nullptr && volume->GetSolid() != nullptr)
        std::cout<<volume->GetSolid()->GetName()<<",XXXX ,OTHER,"<<volume->GetMaterial()->GetName()<<std::endl;

     std::queue<G4LogicalVolume*> toProcess;
     for (size_t i=0;i<volume->GetNoDaughters();++i){
       if (volume->GetDaughter(i)->GetCopyNo()==0) toProcess.push(volume->GetDaughter(i)->GetLogicalVolume());
     }

     while (!toProcess.empty()){
       browseCaloSolids(toProcess.front());
       toProcess.pop();
     }
     return;
  }

}
