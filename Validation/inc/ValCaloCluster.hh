
#ifndef ValCaloCluster_HH_
#define ValCaloCluster_HH_

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValCaloCluster {
 public:
  ValCaloCluster(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const CaloClusterCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _ht;
  TH1D* _hE;
  TH1D* _hR;
};
}  // namespace mu2e

#endif
