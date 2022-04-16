
#ifndef ValCaloDigi_HH_
#define ValCaloDigi_HH_

#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValCaloDigi {
 public:
  ValCaloDigi(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const CaloDigiCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hN2;
  TH1D* _hI;
  TH1D* _ht;
  TH1D* _ht2;
  TH1D* _hm;
  TH1D* _hE;
};
}  // namespace mu2e

#endif
