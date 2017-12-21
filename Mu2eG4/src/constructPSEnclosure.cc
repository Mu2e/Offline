// Andrei Gaponenko, 2012

#include "Mu2eG4/inc/constructPSEnclosure.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Cone.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestCons.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "G4Helper/inc/G4Helper.hh"

#include "G4LogicalVolume.hh"

namespace mu2e {

  //================================================================

  void constructPSEnclosure(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<PSEnclosure> pse;
    GeomHandle<PSVacuum> psv;

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false) 
      || config.getBool("ps.doSurfaceCheck",false);
    const bool placePV             = true;
    const int  verbosityLevel      = config.getInt("PSEnclosure.verbosityLevel",0);
    CLHEP::Hep3Vector extraOffset(0.,0.,pse->getExtraOffset());

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
	       config.getBool("PSEnclosure.visible"),
	       G4Colour::Blue(),
	       config.getBool("PSEnclosure.solid"),
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
	       config.getBool("PSEnclosure.visible"),
	       G4Colour::Blue(),
	       config.getBool("PSEnclosure.solid"),
	       forceAuxEdgeVisible,
	       placePV,
	       doSurfaceCheck
	       );
    }

    // get the mass of the Shell
    G4Helper* _helper = &(*art::ServiceHandle<G4Helper>());
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
                                         config.getBool("PSEnclosure.visible"),
                                         G4Colour::Blue(),
                                         config.getBool("PSEnclosure.solid"),
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

      // Hole in the endPlate for the window
      const CLHEP::Hep3Vector vacCenter =
        pse->windows()[i].originInMu2e() +
        CLHEP::Hep3Vector(0,0, pse->endPlate().halfLength() + pse->windows()[i].halfLength()) + extraOffset
        ;

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
               config.getBool("PSEnclosure.vacuum.visible"),
               G4Colour::Black(),
               config.getBool("PSEnclosure.vacuum.solid"),
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
               config.getBool("PSEnclosure.visible"),
               G4Colour::Grey(),
               config.getBool("PSEnclosure.solid"),
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
