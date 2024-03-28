#ifndef RecoDataProducts_CrvDigi_hh
#define RecoDataProducts_CrvDigi_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include <array>
#include <vector>
#include <cstdint>

namespace mu2e
{
  class CrvDigi
  {
    public:

    static constexpr size_t NSamples = 8; //FIXME: this is also a parameter in CrvDigiMC

    CrvDigi() {}

    CrvDigi(const std::array<int16_t, NSamples> &ADCs, uint16_t startTDC, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, uint8_t SiPMNumber) :
            _ADCs(ADCs), _startTDC(startTDC), _scintillatorBarIndex(scintillatorBarIndex), _SiPMNumber(SiPMNumber) {}

    const std::array<int16_t, NSamples>  &GetADCs() const     {return _ADCs;}
    uint16_t                              GetStartTDC() const {return _startTDC;}

    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    uint8_t                       GetSiPMNumber() const           {return _SiPMNumber;}

    private:

    std::array<int16_t, NSamples>  _ADCs{0};
    uint16_t                       _startTDC{0};

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    uint8_t                        _SiPMNumber{0};
  };
  typedef std::vector<mu2e::CrvDigi> CrvDigiCollection;
}

#endif /* RecoDataProducts_CrvDigi_hh */
