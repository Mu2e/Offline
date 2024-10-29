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

    static constexpr size_t NSamples = 8; //FIXME: this is also a parameter in CrvDigi

    CrvDigiMC() {}
    CrvDigiMC(const std::array<double,NSamples> &voltages, const std::vector<art::Ptr<CrvStep> > &steps,
              art::Ptr<SimParticle> simParticle, double startTime, double TDC0Time,
              mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) :
                          _voltages(voltages),
                          _steps(steps),
                          _simParticle(simParticle),
                          _startTime(startTime),
                          _TDC0Time(TDC0Time),
                          _scintillatorBarIndex(scintillatorBarIndex),
                          _SiPMNumber(SiPMNumber) {}

    const std::array<double,NSamples>         &GetVoltages() const        {return _voltages;}
    const std::vector<art::Ptr<CrvStep> >     &GetCrvSteps() const        {return _steps;}
    const art::Ptr<SimParticle>               &GetSimParticle() const     {return _simParticle;}
    const double                              &GetStartTime() const       {return _startTime;}
    const double                              &GetTDC0Time() const        {return _TDC0Time;}

    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}

    void setSimParticle(const art::Ptr<SimParticle>& sim) {_simParticle = sim;}
    void setCrvSteps(const std::vector<art::Ptr<CrvStep> >& steps) {_steps = steps;}

    private:

    std::array<double,NSamples>         _voltages{0};
    std::vector<art::Ptr<CrvStep> >     _steps  ;        //crv steps responsible for this waveform
    art::Ptr<SimParticle>               _simParticle;    //most likely sim particle responsible for this waveform
    double                              _startTime{0};
    double                              _TDC0Time{0};

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber{0};
  };
  typedef std::vector<mu2e::CrvDigiMC> CrvDigiMCCollection;
}
#endif /* MCDataProducts_CrvDigiMC_hh */
