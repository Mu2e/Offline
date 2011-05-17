#ifndef Mu2eG4_UserTrackInformation_hh
#define Mu2eG4_UserTrackInformation_hh
//
// Mu2e specific information about one G4 track.
//
// $Id: UserTrackInformation.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "ToyDP/inc/ProcessCode.hh"

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

    bool       isForced() const { return _forcedStop; }
    ProcessCode code()    const { return _code; }

    virtual void Print() const;

  private:

    // Did Mu2e user stepping action force a stop?
    bool _forcedStop;

    // If it did, then this is the reason why.
    ProcessCode _code;

  };

} // end namespace mu2e

#endif /* Mu2eG4_UserTrackInformation_hh */
