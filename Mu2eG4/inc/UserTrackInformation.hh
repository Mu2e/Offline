#ifndef UserTrackInformation_H
#define UserTrackInformation_H
//
// Mu2e specific information about one G4 track.
//
// $Id: UserTrackInformation.hh,v 1.1 2010/12/17 22:07:56 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/12/17 22:07:56 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "ToyDP/inc/StoppingCode.hh"

// Geant4 includes
#include "G4VUserTrackInformation.hh"

namespace mu2e{

  class UserTrackInformation : public G4VUserTrackInformation{

  public:
    UserTrackInformation();
    virtual ~UserTrackInformation();

    void setStoppingCode ( StoppingCode code){
      _forcedStop =  true;
      _code = code;
    }

    bool         isForced() const { return _forcedStop; }
    StoppingCode code()     const { return _code; }

    virtual void Print() const;

  private:

    // Did Mu2e stepping action force a stop?
    bool _forcedStop;

    // If it did, then this is the reason why.
    StoppingCode _code;

  };

} // end namespace mu2e

#endif
