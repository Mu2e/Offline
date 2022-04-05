
#ifndef ValProtonBunchTimeMC_HH_
#define ValProtonBunchTimeMC_HH_

#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileDirectory.h"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "TH1D.h"
#include <string>

namespace mu2e {

  class ValProtonBunchTimeMC {

  public:
    ValProtonBunchTimeMC(std::string name):_name(name){}
    int declare( const art::TFileDirectory& tfs);
    int fill(const ProtonBunchTimeMC & obj, art::Event const& event);
    std::string& name() { return _name; }

  private:
    std::string _name;

    TH1D* _hVer;
    TH1D* _htime;

  };
}


#endif
