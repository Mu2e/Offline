
#ifndef ValStrawDigiADCWaveform_HH_
#define ValStrawDigiADCWaveform_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValStrawDigiADCWaveform {

  public:
    ValStrawDigiADCWaveform(std::string name):_name(name){}
    int declare( const art::TFileDirectory& tfs);
    int fill(const StrawDigiADCWaveformCollection & coll, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;

    TH1D* _hVer;
    TH1D* _hN;
    TH1D* _hN2;
    TH1D* _hlen;
    TH1D* _hadc;
    TH1D* _hpmp;

  };
}


#endif
