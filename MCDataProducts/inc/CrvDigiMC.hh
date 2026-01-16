#ifndef MCDataProducts_CrvDigiMC_hh
#define MCDataProducts_CrvDigiMC_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include <vector>

namespace mu2e
{
  class CrvDigiMC
  {
    public:

    CrvDigiMC() {}
    CrvDigiMC(const std::vector<double> &voltages, const std::vector<art::Ptr<CrvStep> > &steps,
              art::Ptr<SimParticle> simParticle, double startTime, double TDC0Time, bool NZS,
              mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) :
                          _voltages(voltages),
                          _steps(steps),
                          _simParticle(simParticle),
                          _startTime(startTime),
                          _TDC0Time(TDC0Time),
                          _NZS(NZS),
                          _scintillatorBarIndex(scintillatorBarIndex),
                          _SiPMNumber(SiPMNumber) {}

    const std::vector<double>                 &GetVoltages() const        {return _voltages;}
    const std::vector<art::Ptr<CrvStep> >     &GetCrvSteps() const        {return _steps;}
    std::vector<art::Ptr<CrvStep> >           &GetCrvSteps()              {return _steps;}  //needed for mixer
    const art::Ptr<SimParticle>               &GetSimParticle() const     {return _simParticle;}
    double                                     GetStartTime() const       {return _startTime;}
    double                                     GetTDC0Time() const        {return _TDC0Time;}
    bool                                       IsNZS() const              {return _NZS;}

    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}

    void setSimParticle(const art::Ptr<SimParticle>& sim) {_simParticle = sim;}
    void setCrvSteps(const std::vector<art::Ptr<CrvStep> >& steps) {_steps = steps;}

    private:

    std::vector<double>                 _voltages;
    std::vector<art::Ptr<CrvStep> >     _steps;          //crv steps responsible for this waveform
    art::Ptr<SimParticle>               _simParticle;    //most likely sim particle responsible for this waveform
    double                              _startTime{0};
    double                              _TDC0Time{0};
    bool                                _NZS{false};     //non-zero suppressed

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber{0};
  };
  typedef std::vector<mu2e::CrvDigiMC> CrvDigiMCCollection;
}
#endif /* MCDataProducts_CrvDigiMC_hh */
