#ifndef Mu2eG4_constructProductionTargetMon_hh
#define Mu2eG4_constructProductionTargetMon_hh

//
// Free function. Approach borrowed from constructPS
// Constructs the downstream production target scanning monitor.
// Parent volume is the air in the target hall. Probably?
//

//Mu2e includes
#include "Mu2eG4Helper/inc/VolumeInfo.hh"

namespace mu2e {

    class VolumeInfo;
    class SimpleConfig;

    void constructProductionTargetMon(VolumeInfo const & parent, SimpleConfig const & _config);

    void constructTargetHallPWC(VolumeInfo const & parent, SimpleConfig const & _config, std::string const & nameSuffix, G4ThreeVector const & position, int const wireNumStart);

} // namespace mu2e


#endif /* Mu2eG4_constructProductionTargetMon_hh */
