#ifndef RecoDataProducts_STMWaveformCollection_hh
#define RecoDataProducts_STMWaveformCollection_hh
//
// The STMWaveformCollection isn't just a vector of STMWaveforms.
// It contains a channel number since we plan to keep the waveforms from each detector separate
//

#include "Offline/RecoDataProducts/inc/STMWaveform.hh"

namespace mu2e {

  class STMWaveformCollection {
  public:
    STMWaveformCollection() : _channel(STMChannel::unknown) {}
    STMWaveformCollection(STMChannel ch) : _channel(ch) {}

    STMChannel channel() const { return _channel; }

    const size_t size() const { return _waveforms.size(); }
    void push_back(const mu2e::STMWaveform& waveform) { _waveforms.push_back(waveform); }
    const std::vector<mu2e::STMWaveform>& waveforms() const { return _waveforms; }

  private:
    mu2e::STMChannel _channel;
    std::vector<mu2e::STMWaveform> _waveforms;
  };
}
#endif
