
#ifndef ValCaloShowerStep_HH_
#define ValCaloShowerStep_HH_

#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValCaloShowerStep {
 public:
  ValCaloShowerStep(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const CaloShowerStepCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hN2;
  TH1D* _ht;
  TH1D* _ht2;
  TH1D* _hE;
  TH1D* _hE2;
  TH1D* _hposx;
  TH1D* _hposy;
  TH1D* _hposz;
};
}  // namespace mu2e

#endif
