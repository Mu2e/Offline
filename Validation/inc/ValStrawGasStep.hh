
#ifndef ValStrawGasStep_HH_
#define ValStrawGasStep_HH_

#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValStrawGasStep {
 public:
  ValStrawGasStep(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const StrawGasStepCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hN2;
  TH1D* _ht;
  TH1D* _ht2;
  TH1D* _hE;
  TH1D* _hlE;
  TH1D* _hlen;
  TH1D* _pmom;
  TH1D* _hz;
  TH1D* _hSI;
  TH1D* _hpla;
  TH1D* _hstr;
};
}  // namespace mu2e

#endif
