#ifndef ITGasLayerSD_ExtWireData_h
#define ITGasLayerSD_ExtWireData_h 1

// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD.hh"

namespace mu2e {

  class ITGasLayerSD_ExtWireData : public ITGasLayerSD {

  public:
    ITGasLayerSD_ExtWireData(G4String);
    ~ITGasLayerSD_ExtWireData(){}
    
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    
  };

} // namespace mu2e

#endif /*ITGasLayerSD_ExtWireData_h*/
