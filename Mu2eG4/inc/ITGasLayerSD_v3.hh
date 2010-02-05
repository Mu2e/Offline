#ifndef ITGasLayerSD_v3_h
#define ITGasLayerSD_v3_h 1

// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD.hh"

namespace mu2e {

class ITGasLayerSD_v3: public ITGasLayerSD {

public:
	ITGasLayerSD_v3(G4String);
	~ITGasLayerSD_v3() {
	}

	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

};

} // namespace mu2e

#endif /*ITGasLayerSD_v3_h*/
