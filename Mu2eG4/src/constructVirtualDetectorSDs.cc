//
// Free function to create the make the virtual detectors sensitive
//
// constructVirtualDetectorSDs.cc
// Author: Lisa Goodenough
// 2018/03/12
//

#include "Offline/Mu2eG4/inc/constructVirtualDetectorSDs.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"
#include "Geant4/G4LogicalVolumeStore.hh"

using namespace std;

namespace mu2e {

  // Construct the virtual detectors

  void constructVirtualDetectorSDs(SimpleConfig const & _config,
                                   Mu2eG4SensitiveDetector* vdSD){

      G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();

      for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++){

          G4String LVname = (*pos)->GetName();

            if (LVname.find("VirtualDetector") != std::string::npos) {
                //std::cout << "Setting this VirtualDetector LV: " << LVname << std::endl;
                store->GetVolume(LVname)->SetSensitiveDetector(vdSD);
            }

      }//for G4LogicalVolumeStore::iterator


  } // constructVirtualDetectorSDs()
}
