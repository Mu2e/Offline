
#ifndef ValTimeCluster_HH_
#define ValTimeCluster_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "RecoDataProducts/inc/TimeCluster.hh"
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
    TH1D* _hx;
    TH1D* _hy;
    TH1D* _hz;
    TH1D* _ht;
    TH1D* _hterr;
    TH1D* _nc;
  };
}


#endif
