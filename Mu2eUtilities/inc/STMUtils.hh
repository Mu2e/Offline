#ifndef Mu2eUtilities_STMUtils_hh
#define Mu2eUtilities_STMUtils_hh
//
// Helper class to do work on STM objects.
//

#include "Offline/DataProducts/inc/STMChannel.hh"
#include "canvas/Utilities/InputTag.h"

namespace mu2e {

  namespace STMUtils {
    // Want to take an input tag and return an STMChannel
    STMChannel getChannel(art::InputTag const& tag);
  }
}

#endif
