#ifndef Mu2eG4_UserTrackInformation_hh
#define Mu2eG4_UserTrackInformation_hh
//
// Mu2e specific information about one G4 track.
//
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "MCDataProducts/inc/ProcessCode.hh"

// Geant4 includes
#include "G4VUserTrackInformation.hh"
#include "G4ThreeVector.hh"

namespace mu2e{

  class UserTrackInformation : public G4VUserTrackInformation{

  public:
    UserTrackInformation();
    virtual ~UserTrackInformation();

    void setProcessCode ( ProcessCode code){
      _forcedStop = true;
      _code = code;
    }

    void setMuCapCode ( ProcessCode code){
      _muCapCode = code;
    }

    bool         isForced() const { return _forcedStop; }
    ProcessCode  code()    const { return _code; }
    ProcessCode  muCapCode() const { return _muCapCode; }

    //  Returns the normalized direction of the momentum
    const G4ThreeVector& GetMomentumDirection() const {return _momDirection;}

    //  Sets the normalized direction of the momentum (no check)
    void SetMomentumDirection(const G4ThreeVector &aDirection) {
      _momDirection = aDirection;
    }

    //  Returns the kinetic energy
    G4double GetKineticEnergy() const {return _kinEnergy;};

    //  Sets the kinetic energy
    void SetKineticEnergy(G4double kEnergy) {
      _kinEnergy = kEnergy;
    }

    virtual void Print() const;

  private:

    // Did Mu2e user stepping action force a stop?
    bool _forcedStop;

    // If it did, then this is the reason why.
    ProcessCode _code;

    // Label of muMinusCaptureAtRest daugter particles (if any)

    ProcessCode _muCapCode;

    // quantities recorded by a mu2e special process before geant4
    // poststepdoits acted

    G4ThreeVector _momDirection;
    G4double _kinEnergy;

  };

} // end namespace mu2e

#endif /* Mu2eG4_UserTrackInformation_hh */
