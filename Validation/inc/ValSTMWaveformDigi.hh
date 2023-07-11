#ifndef ValSTMWaveformDigi_HH_
#define ValSTMWaveformDigi_HH_

#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "TH1D.h"
#include <string>

namespace mu2e {

class ValSTMWaveformDigi {
 public:
  ValSTMWaveformDigi(std::string name) : _name(name) {}
  int declare(const art::TFileDirectory& tfs);
  int fill(const STMWaveformDigiCollection& coll, art::Event const& event);
  std::string& name() { return _name; }

 private:
  std::string _name;

  TH1D* _hVer;
  TH1D* _hNwf;
  TH1D* _hlen;
  TH1D* _hadc;
  TH1D* _hamax;
};
}  // namespace mu2e

#endif
