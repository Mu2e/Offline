#ifndef RecoDataProducts_STMEventHeader_hh_
#define RecoDataProducts_STMEventHeader_hh_

#include "artdaq-core-mu2e/Data/Mu2eEventHeader.hh"

namespace mu2e
{
  class STMEventHeader
  {
  public:
    explicit STMEventHeader() {}
    STMEventHeader(uint64_t ewt, uint8_t mode, uint64_t adcClock, uint64_t dtcClock) :
      _mu2eEventHeader(ewt, mode), _adcClock(adcClock), _dtcClock(dtcClock) { }

    uint64_t eventWindowTag() const { return _mu2eEventHeader.EventWindowTag; }
    uint8_t eventMode() const { return _mu2eEventHeader.EventMode; }

    uint64_t adcClock() const { return _adcClock; }
    uint64_t dtcClock() const { return _dtcClock; }

    // Ordered by EWT, then use remaining header fields to further classify
    bool operator<(STMEventHeader const& other) const {
      if (eventWindowTag() != other.eventWindowTag()) {
        return eventWindowTag() < other.eventWindowTag();
      }
      if (eventMode() != other.eventMode()) {
        return eventMode() < other.eventMode();
      }
      if (adcClock() != other.adcClock()) {
        return adcClock() < other.adcClock();
      }
      return dtcClock() < other.dtcClock();
    }

  private:
    Mu2eEventHeader _mu2eEventHeader;
    uint64_t _adcClock;
    uint64_t _dtcClock;
  };
  typedef std::vector<STMEventHeader> STMEventHeaderCollection;
}  // namespace mu2e


#endif
