#ifndef ITGasLayerSD_Hexagonal_h
#define ITGasLayerSD_Hexagonal_h 1

//// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD.hh"

namespace mu2e {

  class ITGasLayerSD_Hexagonal : public ITGasLayerSD {

  public:
    ITGasLayerSD_Hexagonal(G4String, const SimpleConfig& config);
    ~ITGasLayerSD_Hexagonal();

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  };

} // namespace mu2e

#endif /*ITGasLayerSD_Hexagonal_h*/
