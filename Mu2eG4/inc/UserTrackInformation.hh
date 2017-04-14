#ifndef Mu2eG4_UserTrackInformation_hh
#define Mu2eG4_UserTrackInformation_hh
//
// Mu2e specific information about one G4 track.
//
// $Id: UserTrackInformation.hh,v 1.6 2013/12/02 20:12:58 genser Exp $
// $Author: genser $
// $Date: 2013/12/02 20:12:58 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "MCDataProducts/inc/ProcessCode.hh"

// Geant4 includes
#include "G4VUserTrackInformation.hh"

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

    void setStepInfo (double kineticEnergy, int nSteps ) {
      _preLastStepKE = kineticEnergy;
      _nSteps = nSteps;
    }

    bool         isForced() const { return _forcedStop; }
    ProcessCode  code()    const { return _code; }
    double       preLastStepKE() const { return _preLastStepKE; }
    int          nSteps() const { return _nSteps; }
    ProcessCode  muCapCode() const { return _muCapCode; }

    virtual void Print() const;

  private:

    // Did Mu2e user stepping action force a stop?
    bool _forcedStop;

    // If it did, then this is the reason why.
    ProcessCode _code;

    // Kinetic energy of the particle at the beginning of the last step
    double _preLastStepKE;

    // Number of G4 steps the track if made of
    int _nSteps;

    // Label of muMinusCaptureAtRest daugter particles (if any)

    ProcessCode _muCapCode;

  };

} // end namespace mu2e

#endif /* Mu2eG4_UserTrackInformation_hh */
