
#ifndef ValBkgCluster_HH_
#define ValBkgCluster_HH_

#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValBkgCluster {
 public:
  ValBkgCluster(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const BkgClusterCollection& coll, art::Event const& event);
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
}  // namespace mu2e

#endif
