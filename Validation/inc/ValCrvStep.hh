
#ifndef ValCrvStep_HH_
#define ValCrvStep_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "MCDataProducts/inc/CrvStep.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValCrvStep {

  public:
    ValCrvStep(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const CrvStepCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hb;
    TH1D* _ht;
    TH1D* _hE;
    TH1D* _hposx;
    TH1D* _hposy;
    TH1D* _hposz;
    TH1D* _hp;
    TH1D* _hp2;
    TH1D* _hpL;
  };
}


#endif
