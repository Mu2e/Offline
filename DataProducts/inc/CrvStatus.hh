#ifndef RecoDataProducts_CrvStatus_hh
#define RecoDataProducts_CrvStatus_hh
//
//
// Contact person Simon Corrodi
//

#include <array>
#include <vector>

namespace mu2e
{
  class CrvStatus
  {
    public:

    CrvStatus() {}

    CrvStatus(uint8_t ControllerID, uint16_t WordCount, uint32_t ActiveFEBFlags, uint16_t TriggerCount, uint16_t Status, uint32_t EventWindowTag)
    : _ControllerID(ControllerID)
    , _WordCount(WordCount)
    , _ActiveFEBs(ActiveFEBFlags)
    , _TriggerCount(TriggerCount)
    , _Status(Status)
    , _EventWindowTag(EventWindowTag) {}

    uint8_t  GetDTCId() const {return _ControllerID;}
    uint16_t GetWordCount() const {return _WordCount;}
    uint32_t GetActiveFEBs() const {return _ActiveFEBs;}
    uint16_t GetTriggerCount() const {return _TriggerCount;}
    uint16_t GetStatus() const {return _Status;}
    uint32_t GetEventWindowTag() const {return _EventWindowTag;}

    private:
    uint8_t  _ControllerID;
    uint16_t _WordCount;
    uint32_t _ActiveFEBs;
    uint16_t _TriggerCount;
    uint16_t _Status;
    uint32_t _EventWindowTag;

  };
  typedef std::vector<mu2e::CrvStatus> CrvStatusCollection;
}

#endif /* RecoDataProducts_CrvStatus_hh */
