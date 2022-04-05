
#ifndef ValProtonBunchTime_HH_
#define ValProtonBunchTime_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValProtonBunchTime {

  public:
    ValProtonBunchTime(std::string name):_name(name){}
    int declare( const art::TFileDirectory& tfs);
    int fill(const ProtonBunchTime & obj, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;

    TH1D* _hVer;
    TH1D* _htime;
    TH1D* _hterr;

  };
}


#endif
