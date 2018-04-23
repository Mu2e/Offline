
#ifndef ValCrvCoincidenceCluster_HH_
#define ValCrvCoincidenceCluster_HH_

#include "art/Framework/Principal/Event.h"
#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "Validation/inc/ValPosition.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValCrvCoincidenceCluster {

  public:
    ValCrvCoincidenceCluster(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const CrvCoincidenceClusterCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hSec;
    TH1D* _hPE;
    TH1D* _ht;
    ValPosition _pos;
  };
}


#endif
