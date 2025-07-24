 // Andrei Gaponenko, 2012

#include "Offline/Mu2eG4/inc/constructPSEnclosure.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "Offline/ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/Cone.hh"
#include "Offline/GeomPrimitives/inc/Polycone.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestCons.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"

#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4UnionSolid.hh"

namespace mu2e {

  //================================================================

  void constructPSEnclosure(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<PSEnclosure> pse;
    GeomHandle<PSVacuum> psv;

    CLHEP::Hep3Vector extraOffset(0.,0.,pse->getExtraOffset());

    // Note to self (DNB):  It would be nice to decouple SimpleConfig
    // from Mu2eG4 construct functions.  Pipe dream?

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "psEnclosure", "psEnclosure");
    geomOptions->loadEntry( config, "psEnclosureVacuum", "psEnclosure.vacuum");

    const bool PSIsVisible         = geomOptions->isVisible("psEnclosure");
    const bool PSVacuumIsVisible   = geomOptions->isVisible("psEnclosureVacuum");
    const bool PSIsSolid           = geomOptions->isSolid("psEnclosure");
    const bool PSVacuumIsSolid     = geomOptions->isSolid("psEnclosureVacuum");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("psEnclosure");
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("psEnclosure");
    const bool placePV             = geomOptions->placePV("psEnclosure");
    const int  verbosityLevel      = config.getInt("PSEnclosure.verbosityLevel",0);



    //----------------------------------------------------------------
    std::string sName = "PSEnclosureShell";
    if ( pse->version() > 1 ) {
      double consParams[7] = { pse->shellCone().innerRadius1(),
                               pse->shellCone().outerRadius1(),
                               pse->shellCone().innerRadius2(),
                               pse->shellCone().outerRadius2(),
                               pse->shellCone().halfLength(),
                               pse->shellCone().phi0(),
                               pse->shellCone().deltaPhi() };
      nestCons( sName,
                consParams,
                findMaterialOrThrow(pse->shellCone().materialName()),
                0,
                pse->shellCone().originInMu2e() - parent.centerInMu2e()
                + extraOffset,
                parent,
                0,
               PSIsVisible,
               G4Colour::Blue(),
               PSIsSolid,
               forceAuxEdgeVisible,
               placePV,
               doSurfaceCheck
               );

    } else {  // Old cylindrical version (version 1)
      nestTubs(sName,
               pse->shell().getTubsParams(),
               findMaterialOrThrow(pse->shell().materialName()),
               0,
               pse->shell().originInMu2e() - parent.centerInMu2e(),
               parent,
               0,
               PSIsVisible,
               G4Colour::Blue(),
               PSIsSolid,
               forceAuxEdgeVisible,
               placePV,
               doSurfaceCheck
               );
    } // end adding shell, either conical or cylindrical

    // get the mass of the Shell
    Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>());
    verbosityLevel
      && std::cout << __func__ << " " << sName << " Mass in kg: "
                   << _helper->locateVolInfo(sName).logical->GetMass()/CLHEP::kg
                   << std::endl;


    //==========================================================
    // Append the flange at the end of the shell (versions 3+)
    //==========================================================

    if ( pse->version() > 2 ) {
    sName = "PSEnclosureFlange";
    const VolumeInfo flange =
      nestTubs(sName,
               pse->flange().getTubsParams(),
               findMaterialOrThrow(pse->flange().materialName()),
               0,
               pse->flange().originInMu2e() - parent.centerInMu2e() + extraOffset,
               parent,
               0,
               PSIsVisible,
               G4Colour::Blue(),
               PSIsSolid,
               forceAuxEdgeVisible,
               placePV,
               doSurfaceCheck
               );


    verbosityLevel
      && std::cout << __func__ << " " << sName << " Mass in kg: "
                   << _helper->locateVolInfo(sName).logical->GetMass()/CLHEP::kg
                   << std::endl;

    } // end adding flange


    //========================================================
    // Add the old, flat, endplate
    //========================================================
    if(pse->version() < 3) {
      sName = "PSEnclosureEndPlate";
      const VolumeInfo endPlate = nestTubs(sName,
                                           pse->endPlate().getTubsParams(),
                                           findMaterialOrThrow(pse->endPlate().materialName()),
                                           nullptr,
                                           pse->endPlate().originInMu2e() - parent.centerInMu2e() + extraOffset,
                                           parent,
                                           0,
                                           PSIsVisible,
                                           G4Colour::Blue(),
                                           PSIsSolid,
                                           forceAuxEdgeVisible,
                                           placePV,
                                           doSurfaceCheck
                                           );



      //----------------------------------------------------------------
      // Install the windows
      // MM version ====================================================

      for(unsigned i=0; i<pse->nWindows(); ++i) {

        CLHEP::Hep3Vector windOff(0,0,pse->windows()[i].halfLength());

        // Hole in the endPlate for the window
        const CLHEP::Hep3Vector vacCenter =
          pse->windows()[i].originInMu2e() +
          CLHEP::Hep3Vector(0,0, pse->endPlate().halfLength())   + extraOffset
          + windOff;


        const TubsParams vacTubs(pse->windows()[i].innerRadius(),
                                 pse->windows()[i].outerRadius(),
                                 pse->endPlate().halfLength()
                                 );

        std::ostringstream vname;
        vname<<"PSEnclosureWindowVac"<<i;
        nestTubs(vname.str(),
                 vacTubs,
                 findMaterialOrThrow(psv->vacuum().materialName()),
                 0,
                 vacCenter - endPlate.centerInMu2e(),
                 endPlate,
                 0,
                 PSVacuumIsVisible,
                 G4Colour::Black(),
                 PSVacuumIsSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck
                 );

        // The window itself
        std::ostringstream wname;
        wname<<"PSEnclosureWindow"<<i;

        nestTubs(wname.str(),
                 pse->windows()[i].getTubsParams(),
                 findMaterialOrThrow(pse->windows()[i].materialName()),
                 0,
                 pse->windows()[i].originInMu2e() - parent.centerInMu2e() + extraOffset,
                 parent,
                 0,
                 PSIsVisible,
                 G4Colour::Grey(),
                 PSIsSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck
                 );


        verbosityLevel
          && std::cout << __func__ << " " << wname.str() << " center in Mu2e: "
                       <<pse->windows()[i].originInMu2e()
                       << std::endl;

        //      verbosityLevel
        //        && std::cout << __func__ << " " << sName << " Mass in kg: "
        //                     << _helper->locateVolInfo(sName).logical->GetMass()/CLHEP::kg
        //                     << std::endl;


      } //end window for loop

    } //end version < 3 endplate construction

    else { //build version >= 3 endcap
      auto polycone = pse->endPlatePolycone();
      verbosityLevel
        && std::cout << __func__ << ": PSEndPlate polycone: "
                     << polycone
                     << std::endl;
      G4Polycone* basecone = new G4Polycone("PSEnclosureEndplateBaseCone",
                                            polycone.phi0(),
                                            polycone.phiTotal(),
                                            polycone.numZPlanes(),
                                            &polycone.zPlanes()[0],
                                            &polycone.rInner()[0],
                                            &polycone.rOuter()[0]
                                            );
      G4VSolid* solid = basecone;
      // double zlength = 0.;
      // if(polycone.zPlanes().size() > 0 ) {
      //   zlength = abs(polycone.zPlanes()[0] - polycone.zPlanes()[polycone.zPlanes().size() - 1]); //planes should be ordered
      // }
      CLHEP::Hep3Vector coneOrigin = polycone.originInMu2e();
      // coneOrigin.setZ(coneOrigin.z());
      //punch a hole in the endcap for each window
      for(unsigned iwindow = 0; iwindow<pse->nWindows(); ++iwindow) {
        Tube pipe = pse->windowPipes()[iwindow];
        std::ostringstream pipename, solidname, windowname, boxname;
        pipename << "PSEnclosurePipeHole_" << iwindow+1;
        solidname << "PSEnclosurePlateTmp_" << iwindow+1;
        boxname   << "PSEnclosurePlateBox_" << iwindow+1;
        windowname   << "PSEnclosureWindow_" << iwindow+1;
        //make a solid pipe that is the hole in the pipe, to make a hole in the endplate
        G4Tubs* pipetube = new G4Tubs(pipename.str(), 0., pipe.innerRadius(), pipe.halfLength(), 0., CLHEP::twopi);
        AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
        CLHEP::HepRotation* matrix = reg.add(CLHEP::HepRotation(pipe.rotation()));
        CLHEP::Hep3Vector tubeOrigin = pipe.originInMu2e() - coneOrigin;
        verbosityLevel
          && std::cout << __func__ << ": PSEndPlate tube " << iwindow << ": "
                       << pipe << " origin from end plate: " << tubeOrigin
                       << std::endl;
        solid = new G4SubtractionSolid(solidname.str(), solid, pipetube,
                                       matrix, tubeOrigin );
        //next fuse the pipe to the endplate
        pipetube = new G4Tubs(windowname.str()+"_2", pipe.innerRadius(), pipe.outerRadius(), pipe.halfLength(), 0., CLHEP::twopi);
        solid = new G4UnionSolid(solidname.str()+"_2", solid, pipetube,
                                 matrix, tubeOrigin);
        //next cut off a z-plane face for the window
        CLHEP::Hep3Vector pipeaxis(0.,.0,1.);
        pipeaxis = pipe.rotation()*pipeaxis;
        const double boxside = abs(pipe.outerRadius()/pipeaxis.z())*10.; //cut off everything in this plane essentially
        const double boxlength = pipe.halfLength()*10.;
        G4Box* box = new G4Box(boxname.str(), boxside, boxside, boxlength);
        CLHEP::Hep3Vector windowOrigin = pse->windows()[iwindow].originInMu2e() - coneOrigin;
        CLHEP::Hep3Vector boxOrigin = windowOrigin - CLHEP::Hep3Vector(0.,0.,boxlength-pse->windows()[iwindow].getTubsParams().zHalfLength()-0.01); // gap added
        solid = new G4SubtractionSolid(solidname.str()+"_3", solid, box, nullptr, boxOrigin);
        //now add the window
        windowOrigin = pse->windows()[iwindow].originInMu2e() - parent.centerInMu2e() + extraOffset;
        if(pse->hasWindowFrames()[iwindow]) { //offset the window to be at the outside edge of the window frame if it exists
          windowOrigin -= CLHEP::Hep3Vector(0., 0., (2.*pse->wFramesIn()[iwindow].getTubsParams().zHalfLength()
                                                     + pse->windows()[iwindow].getTubsParams().zHalfLength()-0.01)); // gap added
        }
        nestTubs(windowname.str(),
                 pse->windows()[iwindow].getTubsParams(),
                 findMaterialOrThrow(pse->windows()[iwindow].materialName()),
                 nullptr,
                 windowOrigin,
                 parent,
                 0,
                 PSIsVisible,
                 G4Colour::Grey(),
                 PSIsSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck
                 );

      // The window FRAME if version 3+
      if ( pse->version() >  2 && pse->hasWindowFrames()[iwindow] ) {
        std::ostringstream frname_in, frname_out;
        frname_in <<"PSEnclosureWindowFrameInside" <<iwindow;

        nestTubs(frname_in.str(),
                 pse->wFramesIn()[iwindow].getTubsParams(),
                 findMaterialOrThrow(pse->wFramesIn()[iwindow].materialName()),
                 0,
                 pse->wFramesIn()[iwindow].originInMu2e() - parent.centerInMu2e() + extraOffset,
                 parent,
                 0,
                 PSIsVisible,
                 G4Colour::Grey(),
                 PSIsSolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck
                 );

        if(pse->hasWindowFramesOut()[iwindow]) {
          frname_out<<"PSEnclosureWindowFrameOutside"<<iwindow;
          nestTubs(frname_out.str(),
                   pse->wFramesOut()[iwindow].getTubsParams(),
                   findMaterialOrThrow(pse->wFramesOut()[iwindow].materialName()),
                   0,
                   pse->wFramesOut()[iwindow].originInMu2e() - parent.centerInMu2e() + extraOffset,
                   parent,
                   0,
                   PSIsVisible,
                   G4Colour::Grey(),
                   PSIsSolid,
                   forceAuxEdgeVisible,
                   placePV,
                   doSurfaceCheck
                   );
        }

      } // end adding window FRAMES


        verbosityLevel
          && std::cout << __func__ << ": PSEndPlate window " << iwindow << ": "
                       << pse->windows()[iwindow]
                       << std::endl;
      }
      //ensure no pipes extrude into the PS
      std::vector<double> zplanes(polycone.zPlanes()), rins(polycone.rOuter()), routs;
      for(unsigned iplane = 0; iplane <= polycone.numZPlanes(); ++iplane) {
        routs.push_back(5.e3); //large outer radius to capture everything
      }
      zplanes.push_back(zplanes.back()+5.e3); //add another plane far into the detector
      rins.push_back(0.); //clean everything within this last step
      G4Polycone* outsidecone = new G4Polycone("PSEnclosureEndplateOutline",
                                               polycone.phi0(),
                                               polycone.phiTotal(),
                                               zplanes.size(),
                                               &zplanes[0],
                                               &rins[0],
                                               &routs[0]);
      solid = new G4SubtractionSolid("PSEnclosureEndPlateCleaned", solid, outsidecone);

      VolumeInfo endcap;
      endcap.name = "PSEnclosureEndPlate";
      sName = endcap.name;
      verbosityLevel && std::cout << __func__ << ": Created PS endcap volume, now nesting it...\n";
      endcap.solid = solid;
      finishNesting(endcap,
                    findMaterialOrThrow(pse->endPlatePolycone().materialName()),
                    nullptr,
                    polycone.originInMu2e() - parent.centerInMu2e() + extraOffset,
                    parent.logical,
                    0,
                    PSIsVisible,
                    G4Colour::Blue(),
                    PSIsSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);

    }
    verbosityLevel
      && std::cout << __func__ << " " << sName << " with windows Mass in kg: "
                   << _helper->locateVolInfo(sName).logical->GetMass()/CLHEP::kg
                   << std::endl;

    //----------------------------------------------------------------
  }

}
