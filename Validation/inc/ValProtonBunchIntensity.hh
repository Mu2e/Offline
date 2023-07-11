
#ifndef ValProtonBunchIntensity_HH_
#define ValProtonBunchIntensity_HH_

#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValProtonBunchIntensity {
 public:
  ValProtonBunchIntensity(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const ProtonBunchIntensity& obj, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hint;
};
}  // namespace mu2e

#endif
