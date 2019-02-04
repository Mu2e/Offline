#ifndef MCDataProducts_CrvDigiMC_hh
#define MCDataProducts_CrvDigiMC_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include <vector>

namespace mu2e 
{
  class CrvDigiMC
  {
    public:

    static constexpr size_t NSamples = 8; //FIXME: this is also a parameter in CrvDigi

    CrvDigiMC() {}
    CrvDigiMC(const std::array<double,NSamples> &voltages, const std::vector<art::Ptr<StepPointMC> > &steps, 
              art::Ptr<SimParticle> simParticle, double startTime, 
              mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) :
                          _voltages(voltages), 
                          _steps(steps),
                          _simParticle(simParticle),
                          _startTime(startTime),
                          _scintillatorBarIndex(scintillatorBarIndex),
                          _SiPMNumber(SiPMNumber) {}

    const std::array<double,NSamples>         &GetVoltages() const        {return _voltages;}
    const std::vector<art::Ptr<StepPointMC> > &GetStepPoints() const      {return _steps;}
    const art::Ptr<SimParticle>               &GetSimParticle() const     {return _simParticle;}
    const double                              &GetStartTime() const       {return _startTime;}

    mu2e::CRSScintillatorBarIndex GetScintillatorBarIndex() const {return _scintillatorBarIndex;}
    int                           GetSiPMNumber() const           {return _SiPMNumber;}

    void setSimParticle(const art::Ptr<SimParticle>& sim) {_simParticle = sim;}
    void setStepPoints(const std::vector<art::Ptr<StepPointMC> >& steps) {_steps = steps;}

    private:

    std::array<double,NSamples>         _voltages;
    std::vector<art::Ptr<StepPointMC> > _steps;        //step points responsible for this waveform
    art::Ptr<SimParticle>               _simParticle;  //most likely sim particle responsible for this waveform
    double                              _startTime;

    mu2e::CRSScintillatorBarIndex  _scintillatorBarIndex;
    int                            _SiPMNumber; 
  };
  typedef std::vector<mu2e::CrvDigiMC> CrvDigiMCCollection;
}
#endif /* MCDataProducts_CrvDigiMC_hh */
