//
// Free function to construct the stopping targets.
//
//
// Original author Peter Shanahan
//
// Notes:


// C++ includes
#include <iostream>
#include <string>

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/constructStoppingTarget.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eG4/inc/StrawSD.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestCons.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"
#include "Offline/MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Colour.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/G4VisAttributes.hh"
#include "Geant4/G4LogicalVolumeStore.hh"

using namespace std;

namespace mu2e {

    VolumeInfo constructStoppingTarget( VolumeInfo   const& parent,
                                      SimpleConfig const& config ){

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "stoppingTarget", "stoppingTarget");

    const bool stoppingTargetIsVisible = geomOptions->isVisible("stoppingTarget");
    const bool stoppingTargetIsSolid   = geomOptions->isSolid("stoppingTarget");
    const bool forceAuxEdgeVisible     = geomOptions->forceAuxEdgeVisible("stoppingTarget");
    const bool doSurfaceCheck          = geomOptions->doSurfaceCheck("stoppingTarget");
    const bool placePV                 = geomOptions->placePV("stoppingTarget");

    bool const inGaragePosition = config.getBool("inGaragePosition",false);
    bool const OPA_IPA_ST_Extracted = (inGaragePosition) ? config.getBool("garage.extractOPA_IPA_ST") : false;
    double zOffGarage = (inGaragePosition && OPA_IPA_ST_Extracted) ? config.getDouble("garage.zOffset") : 0.;
    CLHEP::Hep3Vector relPosFake(0.,0., zOffGarage); //for offsetting target in garage position

    int verbosity(config.getInt("stoppingTarget.verbosity",0));

    if ( verbosity > 1 ) std::cout << "In constructStoppingTarget" << std::endl;
    // Master geometry for the Target assembly
    GeomHandle<StoppingTarget> target;

    Mu2eG4Helper    & _helper = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();

    //get the proton absorber elements to prevent overlaps
    bool overlapPabs = false;
    double rAtZClosest, opaZCenter, opaR1, opaR2, opaHL, stZCenter, zclosest;
    if ( config.getBool("hasProtonAbsorber", true) ) {
      GeomHandle<MECOStyleProtonAbsorber> pabs;
      double cylinderRadius = target->cylinderRadius();
      if(pabs->isAvailable(2)) { //there is the OPA in DS2
        opaR1 = pabs->part(2).innerRadiusAtStart();
        opaR2 = pabs->part(2).innerRadiusAtEnd();
        if(!(cylinderRadius < opaR1 && cylinderRadius < opaR2)) { //could overlap
          if(cylinderRadius > opaR1 && cylinderRadius > opaR2) {
            throw cet::exception("GEOM") << "constructStoppingTarget::" << __func__
                                         << ": Stopping target mother overlaps with OPA!\n";
          }
          //could overlap, so check if it does in this z range
          opaZCenter = CLHEP::mm * pabs->part(2).center().z();
          opaHL = pabs->part(2).halfLength();
          stZCenter = target->centerInMu2e().z();
          int side = (opaR1 <= opaR2) ? 1 : -1; //check which way the cone opens
          zclosest = stZCenter - side*target->cylinderLength()/2.;
          //make sure closest point is within OPA region
          if(zclosest > opaZCenter + opaHL) zclosest = opaZCenter+opaHL;
          else if(zclosest < opaZCenter - opaHL) zclosest = opaZCenter-opaHL;
          rAtZClosest = opaR1 + (opaR2-opaR1)*(zclosest - (opaZCenter-opaHL))/(2.*opaHL); //linear radius change in z
          if(rAtZClosest < cylinderRadius + 0.001) overlapPabs = true; //require a small buffer
        } //end possible radius overlap check
      } //end OPA is available check
    }

    TubsParams targetMotherParams(0., target->cylinderRadius(), target->cylinderLength()/2.);

    VolumeInfo targetInfo;
    std::string targetMotherName = "StoppingTargetMother";
    if(!overlapPabs) {
      targetInfo = nestTubs(targetMotherName,
                            targetMotherParams,
                            findMaterialOrThrow(target->fillMaterial()),
                            0,
                            target->centerInMu2e() - parent.centerInMu2e() + relPosFake,
                            parent,
                            0,
                            false/*visible*/,
                            G4Colour::Black(),
                            false/*solid*/,
                            forceAuxEdgeVisible,
                            placePV,
                            doSurfaceCheck
                            );
    } else { //if going to overlap, make a conic section mother volume
      //simply follow the OPA with a small gap
      double motherParams[7] = {0., rAtZClosest-0.001,
                                0., rAtZClosest + abs(opaR2 - opaR1)/(2.*opaHL)*target->cylinderLength() - 0.001,
                                target->cylinderLength()/2.,
                                0.0, 360.0*CLHEP::degree};
      targetInfo = nestCons(targetMotherName,
                            motherParams,
                            findMaterialOrThrow(target->fillMaterial()),
                            0,
                            target->centerInMu2e() - parent.centerInMu2e() + relPosFake,
                            parent,
                            0,
                            false/*visible*/,
                            G4Colour::Black(),
                            false/*solid*/,
                            forceAuxEdgeVisible,
                            placePV,
                            doSurfaceCheck
                            );

    }

    // now create the individual target foils

    G4VPhysicalVolume* pv;

    for (int itf=0; itf<target->nFoils(); ++itf) {

        TargetFoil foil=target->foil(itf);

        VolumeInfo foilInfo;
        G4Material* foilMaterial = findMaterialOrThrow(foil.material());

        std::ostringstream os;
        os << std::setfill('0') << std::setw(2) << itf;
        foilInfo.name = "Foil_" + os.str();

        if ( verbosity > 0 )  std::cout << __func__ << " " << foilInfo.name << std::endl;

        foilInfo.solid = new G4Tubs(foilInfo.name
                                    ,foil.rIn()
                                    ,foil.rOut()
                                    ,foil.halfThickness()
                                    ,0.
                                    ,CLHEP::twopi
                                    );

        foilInfo.logical = new G4LogicalVolume( foilInfo.solid
                                                , foilMaterial
                                                , foilInfo.name
                                                );


        // rotation matrix...
        G4RotationMatrix* rot = 0; //... will have to wait

        G4ThreeVector foilOffset(foil.centerInMu2e() - targetInfo.centerInMu2e() + relPosFake);
        if ( verbosity > 1 ) std::cout << "foil "
                                  << itf
                                  << " centerInMu2e="
                                  << foil.centerInMu2e()
                                  << ", offset="<< foilOffset<< std::endl;

        // G4 manages the lifetime of this object.
        pv = new G4PVPlacement( rot,
                                foilOffset,
                                foilInfo.logical,
                                "Target"+foilInfo.name,
                                targetInfo.logical,
                                0,
                                itf,
                                false);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

        if (!stoppingTargetIsVisible) {
          foilInfo.logical->SetVisAttributes(G4VisAttributes::GetInvisible());
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(stoppingTargetIsSolid);
          foilInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils

    for (int itf=0; itf<target->nSupportStructures(); ++itf) {
        TargetFoilSupportStructure supportStructure=target->supportStructure(itf);
        //TargetFoil foil=target->foil(itf);

        //VolumeInfo foilInfo;
        VolumeInfo supportStructureInfo;

        G4Material* supportStructureMaterial = findMaterialOrThrow(supportStructure.material());

        std::ostringstream os;
        os << std::setfill('0') << std::setw(2) << itf;
        supportStructureInfo.name = "FoilSupportStructure_" + os.str();

        if ( verbosity > 0 )  std::cout << __func__ << " " << supportStructureInfo.name
                                          << std::endl;

        supportStructureInfo.solid = new G4Tubs(supportStructureInfo.name
                                    ,0
                                    ,supportStructure.radius()
                                    ,supportStructure.length()/2. // G4Tubs useses half-lengths as this parameter
                                    ,0.
                                    ,CLHEP::twopi
                                    );

        supportStructureInfo.logical = new G4LogicalVolume( supportStructureInfo.solid
                                                , supportStructureMaterial
                                                , supportStructureInfo.name
                                                );


        if ( verbosity > 1 ) std::cout << "supportStructure.support_id() = "
                                       << supportStructure.support_id()
                                       << "    target->nSupportStructures() = "
                                       << target->nSupportStructures()
                                       << "     target->nFoils() = "
                                       << target->nFoils()
                                       << "     supportStructure.length() = "
                                       << supportStructure.length()
                                       << std::endl;

        // rotation matrices to rotate the orientation of the
        // supporting wires. First rotate into xy-plane by 90deg
        // rotation around y-axis, then rotate within xy-plane by
        // appropiate rotation around z-axis

        CLHEP::HepRotationY secRy(-M_PI/2.);
        CLHEP::HepRotationZ secRz( -supportStructure.support_id() * 360.*CLHEP::deg /
                                   (target->nSupportStructures()/target->nFoils())
                                   - 90.*CLHEP::deg - supportStructure.angleOffset()*CLHEP::deg);
        G4RotationMatrix* supportStructure_rotMatrix = reg.add(G4RotationMatrix(secRy*secRz));

        if ( verbosity > 1 ) std::cout << "supportStructure_rotMatrix = "
                                       << *supportStructure_rotMatrix << std::endl;

        // vector where to place to support tube
        // first find target center
        G4ThreeVector supportStructureOffset(supportStructure.centerInMu2e() - targetInfo.centerInMu2e() + relPosFake);

        if ( verbosity > 1 ) std::cout << supportStructureInfo.name << " "
                                  << itf
                                  << " centerInMu2e="
                                  << supportStructure.centerInMu2e()
                                  << ", offset="
                                  << supportStructureOffset
                                  << std::endl;

        if ( verbosity > 1 ) std::cout << __func__ << " "
                                       << supportStructureInfo.name
                                       <<  std::endl;

        G4ThreeVector
          vector_supportStructure_Orientation( (supportStructure.length()/2.+supportStructure.foil_outer_radius()) *
                                               std::cos(supportStructure.support_id() * 360.*CLHEP::deg /
                                                        (target->nSupportStructures()/target->nFoils()) +
                                                        90.*CLHEP::deg +
                                                        supportStructure.angleOffset()*CLHEP::deg),
                                               (supportStructure.length()/2.+supportStructure.foil_outer_radius()) *
                                               std::sin(supportStructure.support_id() * 360.*CLHEP::deg /
                                                        (target->nSupportStructures()/target->nFoils()) +
                                                        90.*CLHEP::deg + supportStructure.angleOffset()*CLHEP::deg), 0);

        if ( verbosity > 1 ) std::cout << "vector_supportStructure_Orientation = "
                                       << vector_supportStructure_Orientation << std::endl;

        supportStructureOffset += vector_supportStructure_Orientation; // second add vector to support wire tube center

        if ( verbosity > 1 ) std::cout << "supportStructureOffset += vector_supportStructure_Orientation = "
                                       << supportStructureOffset << std::endl;

        pv = new G4PVPlacement( supportStructure_rotMatrix,
                                supportStructureOffset,
                                supportStructureInfo.logical,
                                supportStructureInfo.name,
                                targetInfo.logical,
                                0,
                                itf,
                                false);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

        if (!stoppingTargetIsVisible) {
          supportStructureInfo.logical->SetVisAttributes(G4VisAttributes::GetInvisible());
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Blue()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(stoppingTargetIsSolid);
          supportStructureInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils support structures

    return targetInfo;
  }

} // end namespace mu2e
