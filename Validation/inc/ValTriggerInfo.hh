
#ifndef ValTriggerInfo_HH_
#define ValTriggerInfo_HH_

#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValTriggerInfo {
 public:
  ValTriggerInfo(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const mu2e::TriggerInfo& obj, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hccs;
  TH1D* _htrk;
  TH1D* _hhel;
  TH1D* _hhit;
  TH1D* _hcts;
  TH1D* _hcos;
};
}  // namespace mu2e

#endif
