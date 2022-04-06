
#ifndef ValBkgQual_HH_
#define ValBkgQual_HH_

#include "Offline/RecoDataProducts/inc/BkgQual.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValBkgQual {
 public:
  ValBkgQual(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const BkgQualCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hmva;
  TH1D* _hstat;
};
}  // namespace mu2e

#endif
