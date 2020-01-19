
#ifndef ValStrawDigiMC_HH_
#define ValStrawDigiMC_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValStrawDigiMC {

  public:
    ValStrawDigiMC(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const StrawDigiMCCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    TH1D* _htime0;
    TH1D* _htime1;
    TH1D* _hener;
    TH1D* _henerT;
    TH1D* _hcross;
    TH1D* _hgStep;
    TH1D* _hSI;
  };
}


#endif
