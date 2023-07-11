
#ifndef ValCrvRecoPulse_HH_
#define ValCrvRecoPulse_HH_

#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValCrvRecoPulse {
 public:
  ValCrvRecoPulse(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const CrvRecoPulseCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hN;
  TH1D* _hN2;
  TH1D* _hI;
  TH1D* _hIS;
  TH1D* _hPE;
  TH1D* _hPH;
  TH1D* _ht;
  TH1D* _ht2;
  TH1D* _hChi2;
  TH1D* _hLChi2;
  TH1D* _hLE;
  TH1D* _hLE2;
};
}  // namespace mu2e

#endif
