#ifndef Mu2eG4_ITGasLayerSD_Square_hh
#define Mu2eG4_ITGasLayerSD_Square_hh

//// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD.hh"

namespace mu2e {

  class ITGasLayerSD_Square : public ITGasLayerSD {

  public:
    ITGasLayerSD_Square(G4String, const SimpleConfig& config);
    ~ITGasLayerSD_Square();

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  };

} // namespace mu2e

#endif /* Mu2eG4_ITGasLayerSD_Square_hh */
