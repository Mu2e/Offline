#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"


//#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

#define AGDEBUG(stuff) std::cerr<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;

namespace mu2e {
  void constructExtMonFNAL(const VolumeInfo& parent, const SimpleConfig& config) {
    AGDEBUG("start");

    ExtMonFNAL::ExtMon const & det = *(GeomHandle<ExtMonFNAL::ExtMon>());    

    MaterialFinder materialFinder(config);
    
    VolumeInfo logicalEnclosure = nestBox("ExtMonFNAL",
					  det.logicalEnclosureHalfDim(), 
					  materialFinder.get("extmon_fnal.logicalEnclosureMaterialName"),
					  &det.rotationInParent(),
					  det.offsetInParent(),
					  parent, 0, config.getBool("extmon_fnal.logicalEnclosureVisible"),
					  G4Colour::Red(),
					  false, true, true, true
					  );
    
    for(unsigned iplane = 0; iplane < det.nplanes(); ++iplane) {
      std::ostringstream oss;
      oss<<"EMFSensor"<<iplane;

      std::cout<<"AG: adding sensor at "<<det.sensorOffsetInParent(iplane)<<std::endl;

      VolumeInfo vplane = nestBox(oss.str(),
				  det.sensorHalfSize(iplane),
				  findMaterialOrThrow("G4_Si"),
				  0,
				  det.sensorOffsetInParent(iplane), 
				  logicalEnclosure,
				  0,
				  true, G4Colour::Gray(), true, true, true, true);
    }

    AGDEBUG("end");
  }
}
