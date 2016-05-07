#ifndef Mu2eUtilities_KalRepInstanceNameDecoder_hh
#define Mu2eUtilities_KalRepInstanceNameDecoder_hh
//
// Decode the instance name of a KalRepCollection
// or a KalRepPtrCollection.
//
// Contact person, Rob Kutschke
//
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "Mu2eUtilities/inc/TrkSpecies.hh"

#include <string>

namespace mu2e {

  class KalRepInstanceNameDecoder {
  public:
    KalRepInstanceNameDecoder( std::string const& instanceName );

    TrkFitDirection               direction()    const { return direction_;    }
    TrkSpecies                    particleType() const { return particleType_; }
    int                           charge()       const { return charge_;       }
    std::string const&            instanceName() const { return instanceName_; }

  private:
    std::string                   instanceName_;
    TrkFitDirection               direction_;
    TrkSpecies                    particleType_;
    int                           charge_ = 0;

  };

}

#endif
