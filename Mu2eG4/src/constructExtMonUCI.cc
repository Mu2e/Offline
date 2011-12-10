#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"

#include "GeometryService/inc/GeomHandle.hh"
//#include "GeometryService/inc/Mu2eBuilding.hh"

#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/constructExtMonUCI.hh"

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"

#define YZYDEBUG(stuff) std::cerr<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;

namespace mu2e {
  void constructExtMonUCI(const VolumeInfo& parent, const SimpleConfig& config) {
    YZYDEBUG("start");

    int const verbosity = config.getInt("extmon_uci.verbosity",0);

    const string extmonBaseName = "ExtMonUCI";

    ExtMonUCI::ExtMon const & det = *(GeomHandle<ExtMonUCI::ExtMon>());    
    G4ThreeVector hallOriginInMu2e = parent.centerInMu2e();
    G4ThreeVector detOriginLocal = det.origin() - hallOriginInMu2e;

    MaterialFinder materialFinder(config);
//    G4VSensitiveDetector* emSD = G4SDManager::GetSDMpointer()->
//      FindSensitiveDetector(SensitiveDetectorName::ExtMonUCI());

    bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV       = true;

    bool envelopeVisible = config.getBool("extmon_uci.envelopeVisible",false);
    bool envelopeSolid   = config.getBool("extmon_uci.envelopeSolid",false);

    if (verbosity >= 1)
    {
      YZYDEBUG("ExtMonUCI hallOriginInMu2e " << det.hallOriginInMu2e());
      YZYDEBUG("ExtMonUCI envelope originLocal " << det.originLocal());  
    }

    VolumeInfo envelopeInfo = nestBox( "ExtMonUCI", 
                                       det.envelopeParams(),
                                       materialFinder.get("extmon_uci.envelopeMaterialName"),
                                       &det.rotation(),
                                       detOriginLocal,
                                       parent, 0,
                                       envelopeVisible,
                                       G4Colour::Green(), envelopeSolid,
                                       forceAuxEdgeVisible, placePV, doSurfaceCheck
                                     );

    //
    // Building Collimators
    //
    bool colVisible = config.getBool("extmon_uci.colVisible", true);
    bool colSolid   = config.getBool("extmon_uci.colSolid", false);
    G4Material* colMaterial = materialFinder.get("extmon_uci.colMaterialName");

    const int kNCol = det.nCols();

    VolumeInfo colInfo[kNCol];
    VolumeInfo colShieldInfo[kNCol];
    for (int iCol = 0; iCol < kNCol; iCol++)
    {
      // collimators
      colInfo[iCol].name = det.col(iCol)->name(extmonBaseName);
      if (verbosity >= 2) YZYDEBUG("ExtMonUCI colInfo " << iCol << " name " << colInfo[iCol].name);

      const std::vector<double> colOuterParams = det.col(iCol)->paramsOuter();
      const std::vector<double> colInnerParams = det.col(iCol)->paramsInner();

      if (verbosity >= 2) 
      {
        YZYDEBUG("colOuterParams " << colOuterParams[0] << "  " << colOuterParams[1] << "  " << colOuterParams[2]);
        YZYDEBUG("colInnerParams " << colInnerParams[0] << "  " << colInnerParams[1] << "  " << colInnerParams[2]);
      }

      G4Box *colOuter = new G4Box(colInfo[iCol].name+"_Outer", colOuterParams[0], colOuterParams[1], colOuterParams[2]);
      G4Box *colInner = new G4Box(colInfo[iCol].name+"_Inner", colInnerParams[0], colInnerParams[1], colInnerParams[2]);
      colInfo[iCol].solid = new G4SubtractionSolid(colInfo[iCol].name,
                                                   colOuter,
                                                   colInner,
                                                   0,
                                                   G4ThreeVector(0, 0, 0));
      G4RotationMatrix *colRot = new G4RotationMatrix();
      *colRot = det.col(iCol)->rotation();
      G4ThreeVector colOriginLocal = det.col(iCol)->originLocal();

      int nestVerbosity = 0;
      if (verbosity >= 2) nestVerbosity = 1;
      finishNesting(colInfo[iCol],
                    colMaterial,
                    colRot, 
                    colOriginLocal,
                    envelopeInfo.logical,
                    0,
                    colVisible,
                    G4Color::Blue(),
                    colSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck,
                    nestVerbosity
                   );

      // shielding material outside of collimator
      double shieldScale[3] = {3.0, 3.0, 0.9};
      G4ThreeVector colShieldOuterOrigin = colOriginLocal;
      G4ThreeVector colShieldInnerOrigin(0.0, 0.0, 0.0);

      colShieldInfo[iCol].name = det.col(iCol)->name(extmonBaseName)+"_Shield";
      G4Box *colShieldOuter = new G4Box(colShieldInfo[iCol].name+"Outer", 
                                        colOuterParams[0]*shieldScale[0],
                                        colOuterParams[1]*shieldScale[1],
                                        colOuterParams[2]*shieldScale[2]);
      G4Box *colShieldInner = new G4Box(colShieldInfo[iCol].name+"Inner",
                                        colOuterParams[0],
                                        colOuterParams[1],
                                        colOuterParams[2]);
      colShieldInfo[iCol].solid = new G4SubtractionSolid(colShieldInfo[iCol].name,
                                                         colShieldOuter,
                                                         colShieldInner, 
                                                         colRot,
                                                         colShieldInnerOrigin);
      finishNesting(colShieldInfo[iCol],
                    colMaterial,
                    0,
                    colShieldOuterOrigin,
                    envelopeInfo.logical,
                    0,
                    colVisible,
                    G4Color::Cyan(),
                    colSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck,
                    nestVerbosity
                   );
    }

    YZYDEBUG("end");
  }
}
