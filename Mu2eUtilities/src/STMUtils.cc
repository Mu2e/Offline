// Functions for STMUtils

#include "Offline/Mu2eUtilities/inc/STMUtils.hh"

namespace mu2e {

  namespace STMUtils {
    // Function to get the STMChannel from the art::InputTag
    // (we will keep data from HPGe and LaBr in separate collections)
    STMChannel getChannel(art::InputTag const& tag) {
      if (tag.instance() != "") {
        // If we use instance name, it will only contain the channel name
        return STMChannel(STMChannel::findByName(tag.instance()));
      }
      else {
        std::string label = tag.label();
        // Look at last four characeters of module label to decide channel
        return STMChannel(STMChannel::findByName(label.substr(label.length()-4,4)));
      }
    }
  }
}
