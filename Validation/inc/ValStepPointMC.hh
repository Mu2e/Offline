
#ifndef ValStepPointMC_HH_
#define ValStepPointMC_HH_

#include "art/Framework/Principal/Event.h"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Validation/inc/ValId.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValStepPointMC {

  public:
    ValStepPointMC(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const StepPointMCCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    ValId _id;
    TH1D* _hp;
    TH1D* _het1;
    TH1D* _het2;
    TH1D* _ht;
    TH1D* _hx;
    TH1D* _hy;
    TH1D* _hz;
    TH1D* _hl1;
    TH1D* _hl2;
    TH1D* _hl3;
    TH1D* _hxDS;
    TH1D* _hyDS;
    TH1D* _hxTrk;
    TH1D* _hzTrk;
    TH1D* _hzCal;
    TH1D* _hxCRV;
    TH1D* _hyCRV;
    TH1D* _hzCRV;

  };
}


#endif
