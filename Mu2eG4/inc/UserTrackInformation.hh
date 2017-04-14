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
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"

// Geant4 includes
#include "G4VUserTrackInformation.hh"

#include "CLHEP/Vector/ThreeVector.h"

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

    void setPreLastStepMomentum(CLHEP::Hep3Vector preLastStepMom){
    _preLastStepMom = preLastStepMom;
    }

    void setLastStepPosition(CLHEP::Hep3Vector lastStepPosition){
    _lastStepPosition = lastStepPosition;
    }
    void setEndCode(ProcessCode endCode){
    _endCode = endCode;
    }

    bool               isForced() const { return _forcedStop; }
    ProcessCode        code()    const { return _code; }
    double             preLastStepKE() const { return _preLastStepKE; }
    int                nSteps() const { return _nSteps; }
    ProcessCode        muCapCode() const { return _muCapCode; }
    ProcessCode        endCode() const { return _endCode; }
    CLHEP::Hep3Vector  preLastStepMom() const {return _preLastStepMom;}
    CLHEP::Hep3Vector  lastStepPosition() const {return _lastStepPosition;}
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

    // Label of muMinusCaptureAtRest daughter particles (if any)

    ProcessCode _muCapCode;

    ProcessCode _endCode;

    // momentum just before last step.  This is needed, for example, in
    // modeling antiproton production with non-G4 model.
    CLHEP::Hep3Vector _preLastStepMom;
    //
    // position at last step. When this and the above are incorporated into the SimParticle class these can go away
    CLHEP::Hep3Vector _lastStepPosition;
    //
    // place to write the output
    std::string _fileSave;

  };

} // end namespace mu2e

#endif /* Mu2eG4_UserTrackInformation_hh */
