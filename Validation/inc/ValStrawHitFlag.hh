
#ifndef ValStrawHitFlag_HH_
#define ValStrawHitFlag_HH_

#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValStrawHitFlag {
 public:
  ValStrawHitFlag(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const StrawHitFlagCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hN2;
  TH1D* _hBits;
};
}  // namespace mu2e

#endif
