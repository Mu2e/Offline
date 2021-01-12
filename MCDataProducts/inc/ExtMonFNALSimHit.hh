#ifndef MCDataProducts_ExtMonFNALSimHit_hh
#define MCDataProducts_ExtMonFNALSimHit_hh
//
// A persistable class that keeps sufficient information about a G4
// tracking step in an ExtMonFNAL silicon sensor to produce digitized
// detector hits at a later stage.
//
//
// Original author Andrei Gaponenko
//

#include <ostream>

#include "canvas/Persistency/Common/Ptr.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "DataProducts/inc/ExtMonFNALModuleId.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {

  class ExtMonFNALSimHit{

  public:

    ExtMonFNALSimHit(const ExtMonFNALModuleId&    mid,
                     const art::Ptr<SimParticle>& particle,
                     double                       totalEnergyDeposit,
                     double                       nonIonizingEnergyDeposit,
                     const CLHEP::Hep3Vector&     startPosition,
                     double                       startTime,
                     const CLHEP::Hep3Vector&     endPosition,
                     double                       endTime
                     )
      : moduleId_(mid)
      , particle_(particle)
      , totalEnergyDeposit_(totalEnergyDeposit)
      , nonIonizingEnergyDeposit_(nonIonizingEnergyDeposit)
      , startPosition_(startPosition)
      , startTime_(startTime)
      , endPosition_(endPosition)
      , endTime_(endTime)
    {}

    // Defautl ctr required by genreflex persistency.
    // It should be not used by mu2e code.
    ExtMonFNALSimHit()
      : moduleId_(-1,  -1), totalEnergyDeposit_(), nonIonizingEnergyDeposit_(), startTime_(), endTime_()
    {}

    const ExtMonFNALModuleId& moduleId() const { return moduleId_; }

    const art::Ptr<SimParticle>& simParticle() const { return particle_; }
    art::Ptr<SimParticle>&       simParticle()       { return particle_; }

    double totalEnergyDeposit()        const { return totalEnergyDeposit_; }
    double nonIonizingEnergyDeposit()  const { return nonIonizingEnergyDeposit_; }
    double ionizingEnergyDeposit()     const { return totalEnergyDeposit_ - nonIonizingEnergyDeposit_; }

    // Start and end positions are in the coordinates of the module volume.
    // Start and end times are global.
    const CLHEP::Hep3Vector& localStartPosition() const { return startPosition_;  }
    double                   startTime() const { return startTime_;  }

    const CLHEP::Hep3Vector& localEndPosition() const { return endPosition_;  }
    double                   endTime() const { return endTime_;  }

    // needed for event mixing
    void setSimParticle(const art::Ptr<SimParticle>& p);

  private:
    ExtMonFNALModuleId    moduleId_;

    art::Ptr<SimParticle> particle_;

    double                totalEnergyDeposit_;
    double                nonIonizingEnergyDeposit_;

    CLHEP::Hep3Vector     startPosition_;
    double                startTime_;

    CLHEP::Hep3Vector     endPosition_;
    double                endTime_;
  };

  std::ostream& operator<<(std::ostream& os, const ExtMonFNALSimHit& hit);

} // namespace mu2e

#endif /* MCDataProducts_ExtMonFNALSimHit_hh */
