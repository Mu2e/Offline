#ifndef RecoDataProducts_STMWaveformDigi_hh
#define RecoDataProducts_STMWaveformDigi_hh
//
// Data product that represents the digitized waveforms coming from the STM detectors
// This is used for both unsuppressed and zero-suppressed waveforms
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>
#include <Rtypes.h>
#include <map>

#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/RecoDataProducts/inc/STMEventHeader.hh"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {
  class STMWaveformDigi {

  public:
 // Initialise all variables
    STMWaveformDigi() : _trigTimeOffset(0), _adcs(std::vector<int16_t>()) {};
    // Basic constructor
    STMWaveformDigi(uint32_t trigTimeOffset, std::vector<int16_t> &adcs) : _trigTimeOffset(trigTimeOffset), _adcs(adcs) {};

    uint32_t                    trigTimeOffset() const { return _trigTimeOffset; }
    const std::vector<int16_t>& adcs   () const { return _adcs; }
    void set_data ( size_t n_data, int16_t const* data ) { _adcs.resize(n_data); std::copy(data, data+n_data, _adcs.begin()); }
    art::Ptr<STMWaveformDigi> const& parent() const { return _parent; }
    bool hasParent() const { return _parent.isNonnull(); }
    void setParent(art::Ptr<STMWaveformDigi> const& parent) {_parent = parent; }
  private:
    uint32_t _trigTimeOffset; // time offset from EWT? to first ADC value [ct]
    std::vector<int16_t> _adcs; // vector of ADC values for the waveform
    art::Ptr<STMWaveformDigi> _parent; // get parent raw waveform for zs
  };
  typedef std::vector<STMWaveformDigi> STMWaveformDigiCollection;
  typedef std::map<STMEventHeader,STMWaveformDigiCollection> STMWaveformDigiCollectionMap;
}
#endif
