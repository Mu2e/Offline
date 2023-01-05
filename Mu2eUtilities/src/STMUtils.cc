// Functions for STMUtils

#include "Offline/Mu2eUtilities/inc/STMUtils.hh"

#include "cetlib_except/exception.h"

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

    // A function to get the histogram binning for a certain STMWaveformDigi
    Binning getBinning(const STMWaveformDigi& waveform, const std::string& xAxis, const double nsPerCt) {

      int n_bins = waveform.adcs().size();
      double x_min = 0;
      double x_max = n_bins;
      double t_min = waveform.trigTimeOffset()*nsPerCt;
      double t_max = n_bins*nsPerCt;

      if (xAxis == "sample_number") {
        return Binning(n_bins, x_min, x_max);
      }
      else if (xAxis == "waveform_time") {
        return Binning(n_bins, x_min, t_max);
      }
      else if (xAxis == "event_time") {
        return Binning(n_bins, t_min, t_min+t_max);
      }
      else {
        throw cet::exception("STMUtils::getBinning") << "Invalid xAxis option: \"" << xAxis << "\"" << std::endl;
      }
    }
  }
}
