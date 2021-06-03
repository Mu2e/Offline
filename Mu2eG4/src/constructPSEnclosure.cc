 // Andrei Gaponenko, 2012

#include "Mu2eG4/inc/constructPSEnclosure.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Cone.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"

#include "Geant4/G4LogicalVolume.hh"

namespace mu2e {

  //================================================================

  void constructPSEnclosure(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<PSEnclosure> pse;
    GeomHandle<PSVacuum> psv;

    CLHEP::Hep3Vector extraOffset(0.,0.,pse->getExtraOffset());

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
		

    } else {
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
    }

    // get the mass of the Shell
    Mu2eG4Helper* _helper = &(*art::ServiceHandle<Mu2eG4Helper>());
    verbosityLevel 
      && std::cout << __func__ << " " << sName << " Mass in kg: " 
                   << _helper->locateVolInfo(sName).logical->GetMass()/CLHEP::kg 
                   << std::endl;

    sName = "PSEnclosureEndPlate";
    const VolumeInfo endPlate = nestTubs(sName,
                                         pse->endPlate().getTubsParams(),
                                         findMaterialOrThrow(pse->endPlate().materialName()),
                                         0,
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


    verbosityLevel 
      && std::cout << __func__ << " " << sName << " Mass in kg: " 
                   << _helper->locateVolInfo(sName).logical->GetMass()/CLHEP::kg 
                   << std::endl;
    //----------------------------------------------------------------
    // Install the windows

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
    }

    verbosityLevel 
      && std::cout << __func__ << " " << sName << " with windows Mass in kg: " 
                   << _helper->locateVolInfo(sName).logical->GetMass()/CLHEP::kg 
                   << std::endl;

    //----------------------------------------------------------------
  }

}
