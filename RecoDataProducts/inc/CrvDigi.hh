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

    CrvDigi(const std::vector<int16_t> &ADCs, uint16_t startTDC, bool NZS, bool oddTimestamp, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, uint8_t SiPMNumber,
            uint8_t ROC, uint8_t FEB, uint8_t FEBchannel) :
            _ADCs(ADCs), _startTDC(startTDC), _NZS(NZS), _oddTimestamp(oddTimestamp), _scintillatorBarIndex(scintillatorBarIndex), _SiPMNumber(SiPMNumber),
            _ROC(ROC), _FEB(FEB), _FEBchannel(FEBchannel) {}

    CrvDigi(std::vector<int16_t> &&ADCs, uint16_t startTDC, bool NZS, bool oddTimestamp, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, uint8_t SiPMNumber,
            uint8_t ROC, uint8_t FEB, uint8_t FEBchannel) :
            _ADCs(std::move(ADCs)), _startTDC(startTDC), _NZS(NZS), _oddTimestamp(oddTimestamp), _scintillatorBarIndex(scintillatorBarIndex), _SiPMNumber(SiPMNumber),
            _ROC(ROC), _FEB(FEB), _FEBchannel(FEBchannel) {}

    const std::vector<int16_t>           &GetADCs() const     {return _ADCs;}
    uint16_t                              GetStartTDC() const {return _startTDC;}
    bool                                  IsNZS() const       {return _NZS;}
    bool                                  HasOddTimestamp() const {return _oddTimestamp;}

    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    uint8_t                       GetSiPMNumber() const           {return _SiPMNumber;}
    uint8_t                       GetROC() const                  {return _ROC;}
    uint8_t                       GetFEB() const                  {return _FEB;}
    uint8_t                       GetFEBchannel() const           {return _FEBchannel;}

    private:

    std::vector<int16_t>           _ADCs{0};
    uint16_t                       _startTDC{0};
    bool                           _NZS{false};
    bool                           _oddTimestamp{false};

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    uint8_t                        _SiPMNumber{0};
    uint8_t                        _ROC{0};
    uint8_t                        _FEB{0};
    uint8_t                        _FEBchannel{0};
  };
  typedef std::vector<mu2e::CrvDigi> CrvDigiCollection;
}

#endif /* RecoDataProducts_CrvDigi_hh */
