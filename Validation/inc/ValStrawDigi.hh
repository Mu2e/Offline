
#ifndef ValStrawDigi_HH_
#define ValStrawDigi_HH_

#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValStrawDigi {
 public:
  ValStrawDigi(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const StrawDigiCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hN2;
  TH1D* _htdc;
  TH1D* _htdc2;
  TH1D* _htot;
  TH1D* _hpmp;
  TH1D* _hSI;
};
}  // namespace mu2e

#endif
