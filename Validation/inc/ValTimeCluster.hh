
#ifndef ValTimeCluster_HH_
#define ValTimeCluster_HH_

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "Validation/inc/ValPosition.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValTimeCluster {

  public:
    ValTimeCluster(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const TimeClusterCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hNhit;
    ValPosition _pos;
    TH1D* _ht;
  };
}


#endif
