
#ifndef ValCrvDigiMC_HH_
#define ValCrvDigiMC_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValCrvDigiMC {

  public:
    ValCrvDigiMC(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const CrvDigiMCCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    TH1D* _hI;
    TH1D* _hIS;
    TH1D* _hNS;
    TH1D* _ht;
    TH1D* _hV;
  };
}


#endif
