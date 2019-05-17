#ifndef RecoDataProducts_CrvDigi_hh
#define RecoDataProducts_CrvDigi_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include <array>
#include <vector>

namespace mu2e 
{
  class CrvDigi
  {
    public:

    static constexpr size_t NSamples = 8; //FIXME: this is also a parameter in CrvDigiMC

    CrvDigi() {}

    CrvDigi(const std::array<unsigned int, NSamples> &ADCs, unsigned int startTDC, mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) :
            _ADCs(ADCs), _startTDC(startTDC), _scintillatorBarIndex(scintillatorBarIndex), _SiPMNumber(SiPMNumber) {}

    const std::array<unsigned int, NSamples> &GetADCs() const     {return _ADCs;}
    unsigned int                              GetStartTDC() const {return _startTDC;}

    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}

    private:

    std::array<unsigned int, NSamples> _ADCs;
    unsigned int                       _startTDC;

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber; 
  };
  typedef std::vector<mu2e::CrvDigi> CrvDigiCollection;
}

#endif /* RecoDataProducts_CrvDigi_hh */
