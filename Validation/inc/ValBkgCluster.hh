
#ifndef ValBkgCluster_HH_
#define ValBkgCluster_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/BkgClusterHit.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValBkgCluster {

  public:
    ValBkgCluster(std::string name):_name(name){}
    int declare( art::TFileDirectory tfs);
    int fill(const BkgClusterCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;
    
    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hr;
    TH1D* _ht;
    TH1D* _hd;
    TH1D* _hBits;
  };
}


#endif
