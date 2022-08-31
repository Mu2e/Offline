#ifndef Mu2eG4_Mu2eG4UserTrackInformation_hh
#define Mu2eG4_Mu2eG4UserTrackInformation_hh
//
// Mu2e specific information about one G4 track.
//
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Offline/MCDataProducts/inc/ProcessCode.hh"

// Geant4 includes
#include "Geant4/G4VUserTrackInformation.hh"
#include "Geant4/G4ThreeVector.hh"

namespace mu2e{

  class Mu2eG4UserTrackInformation : public G4VUserTrackInformation{

  public:
    Mu2eG4UserTrackInformation();
    virtual ~Mu2eG4UserTrackInformation();

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

    //  See more comments next to the data members
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

    //  Returns the global time
    G4double GetGlobalTime() const {return _globalTime;};

    //  Sets the global time
    void SetGlobalTime(G4double gTime) {
      _globalTime = gTime;
    }

    //  Returns the proper time
    G4double GetProperTime() const {return _properTime;};

    //  Sets the proper time
    void SetProperTime(G4double pTime) {
      _properTime = pTime;
    }

    // Returns the position
    const G4ThreeVector& GetPosition() const {return _position;}

    // Sets the position
    void SetPosition(const G4ThreeVector &aPosition) {
      _position = aPosition;
    }

    virtual void Print() const;

  private:

    // Did Mu2e user stepping action force a stop?
    bool _forcedStop;

    // If it did, then this is the reason why.
    ProcessCode _code;

    // Label of muMinusCaptureAtRest daugter particles (if any)

    ProcessCode _muCapCode;

    // quantities recorded by a Mu2e special recorder process before Geant4
    // post step doits acted

    G4ThreeVector _momDirection;
    G4double _kinEnergy;
    G4double _globalTime;
    G4double _properTime;
    G4ThreeVector _position;
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4UserTrackInformation_hh */
