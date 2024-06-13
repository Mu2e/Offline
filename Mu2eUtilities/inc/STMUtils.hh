#ifndef Mu2eUtilities_STMUtils_hh
#define Mu2eUtilities_STMUtils_hh
//
// Helper class to do work on STM objects.
//
#include "Offline/GeneralUtilities/inc/Binning.hh"
#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

#include "canvas/Utilities/InputTag.h"

namespace mu2e {

  namespace STMUtils {
    // Want to take an input tag and return an STMChannel
    // TODO - this gets the channel from the last 4 characters of the module name, make this find the channel name by seaching the full str.
    STMChannel getChannel(art::InputTag const& tag);

    // To get the binning for a specific STMWaveformDigi
    Binning getBinning(const STMWaveformDigi& waveform, const std::string& xAxis, const double nsPerCt);

    // To convert from time (in nanoseconds) to clock ticks and vice versa
    unsigned int convertToClockTicks(double time, const STMChannel channel, const STMEnergyCalib& stmEnergyCalib);
    double convertToTime(unsigned int ct, const STMChannel channel, const STMEnergyCalib& stmEnergyCalib);
  }
}

#endif
