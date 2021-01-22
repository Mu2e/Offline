#ifndef MCDataProducts_CrvPhotons_hh
#define MCDataProducts_CrvPhotons_hh
//
//
// Contact person Ralf Ehrlich
//

#include "MCDataProducts/inc/CrvStep.hh"
#include <vector>

namespace mu2e 
{
  class CrvPhotons
  {
    public:

    struct SinglePhoton
    {
      double _time;
      art::Ptr<CrvStep> _step;
      SinglePhoton(double time, art::Ptr<CrvStep> step) : _time(time), _step(step) {}
      SinglePhoton() : _time(NAN) {}  //to make Root happy
    };

    CrvPhotons() : _SiPMNumber(-1) {}  //automatically sets the barindex to invalid

    CrvPhotons(mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber,
               const std::vector<SinglePhoton> &photons) :
                                         _scintillatorBarIndex(scintillatorBarIndex),
                                         _SiPMNumber(SiPMNumber),
                                         _photons(photons) {}

    mu2e::CRSScintillatorBarIndex    GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                              GetSiPMNumber() const           {return _SiPMNumber;}
    std::vector<SinglePhoton>       &GetPhotons()                    {return _photons;} 
    const std::vector<SinglePhoton> &GetPhotons() const              {return _photons;}

    private:

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber; 
    std::vector<SinglePhoton>      _photons;
  };

  typedef std::vector<mu2e::CrvPhotons> CrvPhotonsCollection;
}

#endif /* MCDataProducts_CrvPhotons_hh */
