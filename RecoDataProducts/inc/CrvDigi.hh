#ifndef RecoDataProducts_CrvDigi_hh
#define RecoDataProducts_CrvDigi_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include <vector>
#include <cstdint>

namespace mu2e
{
  class CrvDigi
  {
    public:

    CrvDigi() {}

    CrvDigi(const std::vector<int16_t> &ADCs, uint16_t startTDC, bool NZS, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, uint8_t SiPMNumber) :
            _ADCs(ADCs), _startTDC(startTDC), _NZS(NZS), _scintillatorBarIndex(scintillatorBarIndex), _SiPMNumber(SiPMNumber) {}

    const std::vector<int16_t>           &GetADCs() const     {return _ADCs;}
    uint16_t                              GetStartTDC() const {return _startTDC;}
    bool                                  IsNZS() const       {return _NZS;}

    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    uint8_t                       GetSiPMNumber() const           {return _SiPMNumber;}

    private:

    std::vector<int16_t>           _ADCs{0};
    uint16_t                       _startTDC{0};
    bool                           _NZS{false};

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    uint8_t                        _SiPMNumber{0};
  };
  typedef std::vector<mu2e::CrvDigi> CrvDigiCollection;
}

#endif /* RecoDataProducts_CrvDigi_hh */
