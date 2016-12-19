
#ifndef ValStepPointMC_HH_
#define ValStepPointMC_HH_

#include "art/Framework/Principal/Event.h"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Validation/inc/ValId.hh"
#include "Validation/inc/ValPosition.hh"
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
    ValPosition _pos;
  };
}


#endif
