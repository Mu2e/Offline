
#ifndef ValTriggerResults_HH_
#define ValTriggerResults_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValTriggerResults {

  public:
    ValTriggerResults(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const art::TriggerResults & obj, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hNpath;
    TH1D* _hState;
    TH1D* _hIndex;

  };
}


#endif
